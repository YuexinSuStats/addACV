# =============================================================================
# Compare EM, WEE, and Oracle-WEE for a right-censored additive-risk
# mixture-cure model.
#
# Methods
#   WEE       : Beran-estimated posterior susceptible weights.
#   OracleWEE : true posterior susceptible weights under the data-generating
#               model (not the realized latent cure labels).
#   EM        : logistic-incidence Bernstein-sieve EM from armc_em_right.R.
#
# For every setting this script writes:
#   *_raw.csv                 replication-level estimates and diagnostics
#   *_summary.csv             requested parameter-performance summary
#   *_survival.csv            mean/quantile marginal population-survival curves
#   *_survival_metrics.csv    integrated squared error of the survival curves
#   *_incidence.csv           gridwise incidence means, quantiles, bias, RMSE
#   *_incidence_surface.png   incidence-probability surface comparison plot
#
# The 12 settings are S1/S2/S3 x L/H x n300/n600. Example:
#   SETTING_ID=S1-L-n300 OUTPUT_DIR=compare/S1-L-n300 \
#   REP_START=1 REP_END=500 PARALLEL_BACKEND=mpi \
#   mpirun Rscript compare.R
#
# A quick local pilot can use:
#   SETTING_ID=S1-L-n300 REP_START=1 REP_END=10 BOOT_B=100 \
#   PARALLEL_BACKEND=sequential Rscript compare.R
#
# Keep armc_em_right.R in the working directory, or set EM_CODE_PATH.
# =============================================================================

suppressPackageStartupMessages(library(survival))

PARALLEL_BACKEND <- tolower(Sys.getenv("PARALLEL_BACKEND", "mpi"))
if (!PARALLEL_BACKEND %in% c("mpi", "sequential")) {
  stop("PARALLEL_BACKEND must be 'mpi' or 'sequential'.")
}
if (PARALLEL_BACKEND == "mpi") {
  suppressPackageStartupMessages({
    library(iterators)
    library(Rmpi)
    library(foreach)
    library(doMPI)
  })
}

EM_CODE_PATH <- Sys.getenv("EM_CODE_PATH", "armc_em_right.R")
if (!file.exists(EM_CODE_PATH)) {
  stop(
    "Cannot find armc_em_right.R. Put it in the working directory or set ",
    "EM_CODE_PATH to its location."
  )
}
source(EM_CODE_PATH, local = .GlobalEnv)

timestart <- Sys.time()

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
SETTING_ID_INPUT <- Sys.getenv("SETTING_ID", "S1-L-n300")
# Accept both S1/S2/S3 and pai1/pai2/pai3 prefixes.
SETTING_ID <- sub(
  "^pai([123])",
  "S\\1",
  SETTING_ID_INPUT,
  ignore.case = TRUE
)
REP_START <- as.integer(Sys.getenv("REP_START", "1"))
REP_END <- as.integer(Sys.getenv("REP_END", "500"))
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", file.path("compare_EM_WEE", SETTING_ID))
PARALLEL_SEED_BASE <- as.integer(Sys.getenv("PARALLEL_SEED_BASE", "20260922"))
SEED_INTERNAL_DATA <- as.integer(Sys.getenv("SEED_INTERNAL_DATA", "3"))

valid.settings <- as.vector(outer(
  c("S1", "S2", "S3"),
  c("L-n300", "L-n600", "H-n300", "H-n600"),
  paste,
  sep = "-"
))
if (!SETTING_ID %in% valid.settings) {
  stop("Invalid SETTING_ID. Use S1/S2/S3-L/H-n300/n600.")
}
if (is.na(REP_START) || is.na(REP_END) || REP_START < 1L || REP_END < REP_START) {
  stop("REP_START and REP_END must satisfy 1 <= REP_START <= REP_END.")
}

setting.parts <- strsplit(SETTING_ID, "-", fixed = TRUE)[[1]]
pai.setting <- setting.parts[1]
censor.setting <- setting.parts[2]
n <- as.integer(sub("^n", "", setting.parts[3]))
rep.ids <- seq.int(REP_START, REP_END)
repetition <- length(rep.ids)

beta.true <- c(1, 1)
if (censor.setting == "L") {
  lambdaC <- 0.50
  gamma.true <- c(-1.5, 1, 1)
} else {
  lambdaC <- 1.00
  gamma.true <- c(-2.5, 1, 1)
}

pai <- switch(
  pai.setting,
  S1 = function(u) plogis(u),
  S2 = function(u) plogis(0.4 * (u^3 - u^2 - 0.9 * u + 1)),
  S3 = function(u) exp(-exp(-u))
)

tau0 <- 2
kappa <- 1

BOOT_B_ENV <- Sys.getenv("BOOT_B", "")
BOOT_B <- if (nzchar(BOOT_B_ENV)) as.integer(BOOT_B_ENV) else 4L * n
MIN_VALID_BOOT_ENV <- Sys.getenv("MIN_VALID_BOOT", "")
MIN_VALID_BOOT <- if (nzchar(MIN_VALID_BOOT_ENV)) {
  as.integer(MIN_VALID_BOOT_ENV)
} else {
  max(30L, ceiling(0.80 * BOOT_B))
}
if (is.na(BOOT_B) || BOOT_B < 2L) stop("BOOT_B must be at least 2.")
if (is.na(MIN_VALID_BOOT) || MIN_VALID_BOOT < 2L || MIN_VALID_BOOT > BOOT_B) {
  stop("MIN_VALID_BOOT must lie between 2 and BOOT_B.")
}

EM_DEGREE_ENV <- Sys.getenv("EM_DEGREE", "")
EM_DEGREE <- if (nzchar(EM_DEGREE_ENV)) {
  as.integer(EM_DEGREE_ENV)
} else {
  max(3L, floor(n^(1 / 4)))
}
EM_MAXIT <- as.integer(Sys.getenv("EM_MAXIT", "2000"))
EM_TOL <- as.numeric(Sys.getenv("EM_TOL", "1e-6"))

TIME_GRID_N <- as.integer(Sys.getenv("TIME_GRID_N", "101"))
Z_GRID_N <- as.integer(Sys.getenv("Z_GRID_N", "41"))
time.grid <- seq(0, tau0, length.out = TIME_GRID_N)
if (is.na(Z_GRID_N) || Z_GRID_N < 2L) stop("Z_GRID_N must be at least 2.")

# Midpoint grid approximates expectation over Z1~U(0,2), Z2~U(0,4), while
# avoiding kernel-boundary evaluation exactly at 0, 2, or 4.
z1.grid <- (seq_len(Z_GRID_N) - 0.5) * 2 / Z_GRID_N
z2.grid <- (seq_len(Z_GRID_N) - 0.5) * 4 / Z_GRID_N
Z.eval <- as.matrix(expand.grid(Z1 = z1.grid, Z2 = z2.grid))

# -----------------------------------------------------------------------------
# Data generation
# -----------------------------------------------------------------------------
H0.true <- function(t) {
  t <- as.numeric(t)
  ans <- rep(Inf, length(t))
  id <- is.finite(t) & t >= 0 & t < tau0
  ans[id] <- 2 * log1p(t[id]) - kappa * log1p(-t[id] / tau0)
  ans
}

generate_data <- function(nn, seed) {
  set.seed(seed)
  Z <- cbind(
    Z1 = runif(nn, 0, 2),
    Z2 = runif(nn, 0, 4)
  )
  u <- runif(nn)
  T0 <- numeric(nn)
  for (i in seq_len(nn)) {
    zb <- sum(Z[i, ] * beta.true)
    T0[i] <- uniroot(
      function(t) 1 - exp(-H0.true(t) - zb * t) - u[i],
      interval = c(0, tau0),
      tol = .Machine$double.eps^0.5
    )$root
  }

  pi0 <- as.numeric(pai(drop(cbind(1, Z) %*% gamma.true)))
  eta <- rbinom(nn, 1, pi0)
  C <- rexp(nn, lambdaC)
  TT <- ifelse(eta == 1, T0, Inf)
  Y <- pmin(TT, C)
  delta <- as.integer(eta == 1 & T0 <= C)

  dat <- data.frame(Y = Y, delta = delta, Z1 = Z[, 1], Z2 = Z[, 2])
  dat <- dat[order(dat$Y), , drop = FALSE]
  rownames(dat) <- NULL
  list(
    data = dat,
    censor.rate = mean(delta == 0),
    cure.rate = mean(eta == 0)
  )
}

# -----------------------------------------------------------------------------
# Beran weights copied from addICV.R
# -----------------------------------------------------------------------------
KK <- function(x, hn) {
  3 / 4 * (1 - (x / hn)^2) * (abs(x / hn) <= 1)
}

NWweight2 <- function(Z, Z.sample, hn, lam = rep(1, nrow(Z.sample))) {
  Z <- as.matrix(Z)
  Z.sample <- as.matrix(Z.sample)
  m <- nrow(Z)
  p <- ncol(Z)
  nn <- nrow(Z.sample)
  if (ncol(Z.sample) != p || length(hn) != p) stop("Kernel dimensions disagree.")

  Kprod <- matrix(1, nrow = m, ncol = nn)
  for (j in seq_len(p)) {
    distance <- outer(Z[, j], Z.sample[, j], "-")
    Kprod <- Kprod * KK(distance, hn[j])
  }
  weighted <- sweep(Kprod, 2, lam, `*`)
  rs <- rowSums(weighted)
  zero <- !is.finite(rs) | rs <= 0
  rs[zero] <- 1
  ans <- weighted / rs
  ans[zero, ] <- 0
  ans
}

S.Ti <- function(data, w) {
  m <- nrow(w)
  nn <- ncol(w)
  k <- cbind(rep(0, m), w[, seq_len(nn - 1L), drop = FALSE])
  cumsum.k <- t(apply(k, 1, cumsum))
  cum.w <- 1 - cumsum.k
  cum.w <- cum.w * (cum.w > w + 1e-5) + (w + 1e-4) * (cum.w <= w + 1e-5)
  q <- 1 - w / cum.w
  delta.mat <- matrix(data$delta, nrow = m, ncol = nn, byrow = TRUE)
  q[delta.mat == 0] <- 1
  q <- t(apply(q, 1, cumprod))
  r <- rle(data$Y)
  last.index <- rep(cumsum(r$lengths), r$lengths)
  q[cbind(seq_len(m), last.index[seq_len(m)])]
}

P.T <- function(data, w) {
  m <- nrow(w)
  nn <- ncol(w)
  k <- cbind(rep(0, m), w[, seq_len(nn - 1L), drop = FALSE])
  cumsum.k <- t(apply(k, 1, cumsum))
  cum.w <- 1 - cumsum.k
  w[cum.w == 0] <- 0
  cum.w[cum.w == 0] <- 1
  q <- 1 - w / cum.w
  delta.mat <- matrix(data$delta, nrow = m, ncol = nn, byrow = TRUE)
  q[delta.mat == 0] <- 1
  q <- apply(q, 1, prod)
  as.numeric(pmax(q, 0))
}

bandwidth_for <- function(Z) {
  h <- 3 * nrow(Z)^(-2 / 7) * apply(Z, 2, stats::sd)
  if (any(!is.finite(h)) || any(h <= 0)) {
    stop("All bandwidths must be finite and strictly positive.")
  }
  h
}

beran_posterior_weights <- function(data, h) {
  Z <- as.matrix(data[, c("Z1", "Z2")])
  w.kernel <- NWweight2(Z, Z, h)
  S.beran <- S.Ti(data, w.kernel)
  P.inf <- P.T(data, w.kernel)
  S.safe <- pmax(S.beran, 1e-8)
  last.failure <- max(data$Y[data$delta == 1])
  posterior <- data$delta +
    (1 - data$delta) * (S.beran - P.inf) / S.safe *
    (data$Y <= last.failure)
  posterior <- pmin(pmax(posterior, 0), 1)
  list(weight = posterior, P.inf = P.inf, h = h, last.failure = last.failure)
}

true_posterior_weights <- function(data) {
  Z <- as.matrix(data[, c("Z1", "Z2")])
  pi0 <- as.numeric(pai(drop(cbind(1, Z) %*% gamma.true)))
  Su0 <- exp(-H0.true(data$Y) - data$Y * drop(Z %*% beta.true))
  Su0[data$Y >= tau0] <- 0
  denominator <- 1 - pi0 + pi0 * Su0
  posterior.censored <- ifelse(denominator > 0, pi0 * Su0 / denominator, 0)
  out <- ifelse(data$delta == 1, 1, posterior.censored)
  pmin(pmax(out, 0), 1)
}

# -----------------------------------------------------------------------------
# Weighted estimating equation given susceptible weights
# -----------------------------------------------------------------------------
fit_weighted_ee <- function(data, case.weight, eval.times = NULL) {
  ord <- order(data$Y)
  data <- data[ord, , drop = FALSE]
  case.weight <- as.numeric(case.weight)[ord]
  Y <- data$Y
  delta <- data$delta
  Z <- as.matrix(data[, c("Z1", "Z2")])
  nn <- nrow(data)
  p <- ncol(Z)
  if (length(case.weight) != nn || any(!is.finite(case.weight))) {
    stop("Invalid susceptible weights.")
  }

  dif <- diff(c(0, Y))
  rle.Y <- rle(Y)
  start.index <- c(1L, 1L + head(cumsum(rle.Y$lengths), -1L))
  first.index <- rep(start.index, rle.Y$lengths)

  Zbar <- matrix(0, nrow = nn, ncol = p)
  for (j in seq_len(nn)) {
    idx <- first.index[j]:nn
    ww <- case.weight[idx]
    sw <- sum(ww)
    if (sw > 0) Zbar[j, ] <- colSums(Z[idx, , drop = FALSE] * ww) / sw
  }

  numerator <- colSums(case.weight * delta * (Z - Zbar))
  denominator <- matrix(0, nrow = p, ncol = p)
  for (j in seq_len(nn)) {
    if (dif[j] <= 1e-8) next
    idx <- first.index[j]:nn
    ww <- case.weight[idx]
    Zc <- sweep(Z[idx, , drop = FALSE], 2, Zbar[j, ], "-")
    denominator <- denominator +
      (crossprod(Zc, Zc * ww)) * dif[j]
  }
  if (any(!is.finite(denominator)) || rcond(denominator) < 1e-12) {
    stop("Singular WEE estimating-equation matrix.")
  }
  beta.hat <- as.numeric(solve(denominator, numerator))
  if (any(!is.finite(beta.hat))) stop("Non-finite WEE beta estimate.")

  risk.w <- rev(cumsum(rev(case.weight)))
  zbeta <- drop(Z %*% beta.hat)
  risk.wz <- rev(cumsum(rev(case.weight * zbeta)))
  risk.event <- risk.w[first.index]

  H0.at <- function(ts) {
    if (ts <= 0) return(0)
    # The susceptible distribution has the known finite endpoint tau0.  Do not
    # truncate at the random largest observed failure, which would force the
    # estimated susceptible survival to zero too early.
    if (ts >= tau0) return(Inf)
    event.num <- delta * case.weight * (Y <= ts)
    if (any(event.num != 0 & risk.event <= 1e-12)) {
      stop("Zero risk denominator at an event time.")
    }
    term1 <- numeric(nn)
    ok.event <- risk.event > 1e-12
    term1[ok.event] <- event.num[ok.event] / risk.event[ok.event]
    r1 <- sum(term1)

    point <- findInterval(ts, Y)
    if (point == 0L) {
      if (risk.w[1] <= 1e-12) stop("Zero initial risk denominator.")
      r2 <- ts * risk.wz[1] / risk.w[1]
    } else {
      id <- seq_len(point)
      valid <- dif[id] > 0 & risk.w[id] > 1e-12
      if (any(dif[id] > 0 & !valid)) stop("Zero risk denominator in H0 integral.")
      r2 <- sum(dif[id][valid] * risk.wz[id][valid] / risk.w[id][valid])
      if (point < nn && ts > Y[point]) {
        if (risk.w[point + 1L] <= 1e-12) stop("Zero tail risk denominator.")
        r2 <- r2 +
          (ts - Y[point]) * risk.wz[point + 1L] / risk.w[point + 1L]
      }
    }
    ans <- r1 - r2
    if (!is.finite(ans)) stop("Non-finite baseline cumulative hazard.")
    ans
  }

  H0.grid <- NULL
  if (!is.null(eval.times)) {
    H0.grid <- vapply(eval.times, H0.at, numeric(1))
    finite <- is.finite(H0.grid)
    # Shape projection keeps the exported survival-function summaries valid;
    # beta estimates and their inference remain completely unchanged.
    H0.grid[finite] <- cummax(pmax(H0.grid[finite], 0))
  }

  list(
    beta = beta.hat,
    H0 = H0.grid,
    H0.fun = H0.at,
    last.failure = max(Y[delta == 1]),
    weight = case.weight
  )
}

fit_wee <- function(data, mode = c("WEE", "OracleWEE"), eval.times = NULL) {
  mode <- match.arg(mode)
  data <- data[order(data$Y), , drop = FALSE]
  Z <- as.matrix(data[, c("Z1", "Z2")])
  h <- bandwidth_for(Z)
  if (mode == "WEE") {
    weight.info <- beran_posterior_weights(data, h)
    case.weight <- weight.info$weight
    P.inf <- weight.info$P.inf
  } else {
    case.weight <- true_posterior_weights(data)
    P.inf <- NULL
  }
  fit <- fit_weighted_ee(data, case.weight, eval.times)
  fit$data <- data
  fit$h <- h
  fit$P.inf <- P.inf
  fit$mode <- mode
  fit
}

evaluate_weighted_H0 <- function(fit, eval.times) {
  H0.grid <- vapply(eval.times, fit$H0.fun, numeric(1))
  finite <- is.finite(H0.grid)
  H0.grid[finite] <- cummax(pmax(H0.grid[finite], 0))
  H0.grid
}

# -----------------------------------------------------------------------------
# Population survival curves
# -----------------------------------------------------------------------------
predict_pi_wee <- function(fit, Z.new) {
  Z.sample <- as.matrix(fit$data[, c("Z1", "Z2")])
  w.new <- NWweight2(Z.new, Z.sample, fit$h)
  pmin(pmax(1 - P.T(fit$data, w.new), 0), 1)
}

predict_pi_oracle <- function(fit, Z.new) {
  Z.sample <- as.matrix(fit$data[, c("Z1", "Z2")])
  w.new <- NWweight2(Z.new, Z.sample, fit$h)
  pmin(pmax(drop(w.new %*% fit$weight), 0), 1)
}

weighted_population_curve <- function(fit, Z.new, eval.times) {
  pi.hat <- if (fit$mode == "WEE") {
    predict_pi_wee(fit, Z.new)
  } else {
    predict_pi_oracle(fit, Z.new)
  }
  zb <- drop(Z.new %*% fit$beta)
  raw.cumhaz <- outer(fit$H0, rep(1, length(zb))) + outer(eval.times, zb)
  raw.invalid <- any(raw.cumhaz[is.finite(raw.cumhaz)] < -1e-8, na.rm = TRUE) ||
    any(vapply(seq_len(ncol(raw.cumhaz)), function(j) {
      x <- raw.cumhaz[, j]
      x <- x[is.finite(x)]
      length(x) >= 2L && any(diff(x) < -1e-8)
    }, logical(1)))
  projected.cumhaz <- raw.cumhaz
  for (j in seq_len(ncol(projected.cumhaz))) {
    finite <- is.finite(projected.cumhaz[, j])
    projected.cumhaz[finite, j] <- cummax(pmax(projected.cumhaz[finite, j], 0))
  }
  Su <- exp(-projected.cumhaz)
  Su[eval.times >= tau0, ] <- 0
  Spop <- sweep(Su, 2, pi.hat, `*`)
  Spop <- sweep(Spop, 2, 1 - pi.hat, `+`)
  list(curve = rowMeans(Spop), projected = raw.invalid)
}

em_population_curve <- function(fit, Z.new, eval.times) {
  pi.hat <- plogis(drop(cbind(1, Z.new) %*% fit$gamma))
  basis.times <- pmin(eval.times, fit$tau)
  basis <- bernstein_hazard_basis(basis.times, fit$degree, fit$tau)
  H0.hat <- drop(basis$cumulative %*% fit$theta)
  zb <- drop(Z.new %*% fit$beta)
  cumhaz <- outer(H0.hat, rep(1, length(zb))) + outer(eval.times, zb)
  Su <- exp(-cumhaz)
  Su[eval.times >= tau0, ] <- 0
  Spop <- sweep(Su, 2, pi.hat, `*`)
  Spop <- sweep(Spop, 2, 1 - pi.hat, `+`)
  rowMeans(Spop)
}

true_population_curve <- function(Z.new, eval.times) {
  pi0 <- as.numeric(pai(drop(cbind(1, Z.new) %*% gamma.true)))
  zb <- drop(Z.new %*% beta.true)
  H0 <- H0.true(eval.times)
  cumhaz <- outer(H0, rep(1, length(zb))) + outer(eval.times, zb)
  Su <- exp(-cumhaz)
  Su[eval.times >= tau0, ] <- 0
  Spop <- sweep(Su, 2, pi0, `*`)
  Spop <- sweep(Spop, 2, 1 - pi0, `+`)
  rowMeans(Spop)
}

truth.curve <- true_population_curve(Z.eval, time.grid)

# -----------------------------------------------------------------------------
# Standard errors: nonparametric subject bootstrap for all three methods
# -----------------------------------------------------------------------------
fit_em <- function(data) {
  em_armc_right(
    Y = data$Y,
    delta = data$delta,
    Z = as.matrix(data[, c("Z1", "Z2")]),
    degree = EM_DEGREE,
    support.end = tau0,
    control = list(maxit = EM_MAXIT, tol = EM_TOL, trace = FALSE)
  )
}

joint_bootstrap <- function(data, B, min.valid, run.EM = TRUE) {
  nn <- nrow(data)
  beta.WEE <- matrix(NA_real_, nrow = B, ncol = 2)
  beta.Oracle <- matrix(NA_real_, nrow = B, ncol = 2)
  beta.EM <- matrix(NA_real_, nrow = B, ncol = 2)
  time.WEE <- 0
  time.Oracle <- 0
  time.EM <- 0

  for (b in seq_len(B)) {
    # All methods use the same resampled subject indices in each bootstrap
    # draw, which makes their Monte Carlo comparison directly comparable.
    bs <- data[sample.int(nn, nn, replace = TRUE), , drop = FALSE]

    start <- proc.time()[["elapsed"]]
    fw <- tryCatch(fit_wee(bs, "WEE")$beta, error = function(e) NULL)
    time.WEE <- time.WEE + proc.time()[["elapsed"]] - start
    if (!is.null(fw) && all(is.finite(fw))) beta.WEE[b, ] <- fw

    start <- proc.time()[["elapsed"]]
    fo <- tryCatch(fit_wee(bs, "OracleWEE")$beta, error = function(e) NULL)
    time.Oracle <- time.Oracle + proc.time()[["elapsed"]] - start
    if (!is.null(fo) && all(is.finite(fo))) beta.Oracle[b, ] <- fo

    if (isTRUE(run.EM)) {
      start <- proc.time()[["elapsed"]]
      fe <- tryCatch(fit_em(bs), error = function(e) NULL)
      time.EM <- time.EM + proc.time()[["elapsed"]] - start
      if (!is.null(fe) && isTRUE(fe$converged) && all(is.finite(fe$beta))) {
        beta.EM[b, ] <- fe$beta
      }
    }
  }

  summarize <- function(mat) {
    ok <- complete.cases(mat)
    list(
      se = if (sum(ok) >= min.valid) apply(mat[ok, , drop = FALSE], 2, sd) else rep(NA_real_, 2),
      nvalid = sum(ok),
      ok = sum(ok) >= min.valid
    )
  }
  list(
    WEE = c(summarize(beta.WEE), list(time = time.WEE)),
    OracleWEE = c(summarize(beta.Oracle), list(time = time.Oracle)),
    EM = c(summarize(beta.EM), list(time = time.EM))
  )
}

safe_timed <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- tryCatch(
    list(ok = TRUE, value = force(expr), error = NA_character_),
    error = function(e) list(ok = FALSE, value = NULL, error = conditionMessage(e))
  )
  value$time <- proc.time()[["elapsed"]] - start
  value
}

# -----------------------------------------------------------------------------
# One Monte Carlo replication
# -----------------------------------------------------------------------------
one_replication <- function(rep.id) {
  sim <- generate_data(n, SEED_INTERNAL_DATA + rep.id)
  data <- sim$data

  WEE.fit <- safe_timed(fit_wee(data, "WEE"))
  Oracle.fit <- safe_timed(fit_wee(data, "OracleWEE"))

  EM.fit <- safe_timed(fit_em(data))
  EM.converged <- isTRUE(EM.fit$ok) && isTRUE(EM.fit$value$converged)
  EM.point.ok <- EM.converged && all(is.finite(EM.fit$value$beta))

  # The same bootstrap resamples are used for WEE, Oracle-WEE, and EM. Their
  # method-specific elapsed times are accumulated separately.
  set.seed(PARALLEL_SEED_BASE + rep.id)
  boot <- joint_bootstrap(
    data,
    BOOT_B,
    MIN_VALID_BOOT,
    run.EM = EM.point.ok
  )

  EM.se.beta <- boot$EM$se
  EM.se.ok <- EM.point.ok && isTRUE(boot$EM$ok) && all(is.finite(EM.se.beta))

  WEE.point.ok <- isTRUE(WEE.fit$ok) && all(is.finite(WEE.fit$value$beta))
  Oracle.point.ok <- isTRUE(Oracle.fit$ok) && all(is.finite(Oracle.fit$value$beta))
  WEE.se.ok <- WEE.point.ok && isTRUE(boot$WEE$ok) && all(is.finite(boot$WEE$se))
  Oracle.se.ok <- Oracle.point.ok && isTRUE(boot$OracleWEE$ok) &&
    all(is.finite(boot$OracleWEE$se))

  curve.WEE <- curve.Oracle <- curve.EM <- rep(NA_real_, length(time.grid))
  incidence.WEE <- incidence.Oracle <- incidence.EM <-
    rep(NA_real_, nrow(Z.eval))
  projection.WEE <- projection.Oracle <- NA
  if (WEE.point.ok) {
    incidence.WEE <- tryCatch(
      predict_pi_wee(WEE.fit$value, Z.eval),
      error = function(e) rep(NA_real_, nrow(Z.eval))
    )
    WEE.fit$value$H0 <- tryCatch(
      evaluate_weighted_H0(WEE.fit$value, time.grid),
      error = function(e) rep(NA_real_, length(time.grid))
    )
    cw <- tryCatch(
      weighted_population_curve(WEE.fit$value, Z.eval, time.grid),
      error = function(e) NULL
    )
    if (!is.null(cw)) {
      curve.WEE <- cw$curve
      projection.WEE <- cw$projected
    }
  }
  if (Oracle.point.ok) {
    incidence.Oracle <- tryCatch(
      predict_pi_oracle(Oracle.fit$value, Z.eval),
      error = function(e) rep(NA_real_, nrow(Z.eval))
    )
    Oracle.fit$value$H0 <- tryCatch(
      evaluate_weighted_H0(Oracle.fit$value, time.grid),
      error = function(e) rep(NA_real_, length(time.grid))
    )
    co <- tryCatch(
      weighted_population_curve(Oracle.fit$value, Z.eval, time.grid),
      error = function(e) NULL
    )
    if (!is.null(co)) {
      curve.Oracle <- co$curve
      projection.Oracle <- co$projected
    }
  }
  if (EM.point.ok) {
    incidence.EM <- tryCatch(
      plogis(drop(cbind(1, Z.eval) %*% EM.fit$value$gamma)),
      error = function(e) rep(NA_real_, nrow(Z.eval))
    )
    curve.EM <- tryCatch(
      em_population_curve(EM.fit$value, Z.eval, time.grid),
      error = function(e) rep(NA_real_, length(time.grid))
    )
  }

  list(
    rep = rep.id,
    fatal.error = NA_character_,
    censor.rate = sim$censor.rate,
    cure.rate = sim$cure.rate,

    beta1.WEE = if (WEE.point.ok) WEE.fit$value$beta[1] else NA_real_,
    beta2.WEE = if (WEE.point.ok) WEE.fit$value$beta[2] else NA_real_,
    se1.WEE = if (WEE.se.ok) boot$WEE$se[1] else NA_real_,
    se2.WEE = if (WEE.se.ok) boot$WEE$se[2] else NA_real_,
    point.ok.WEE = WEE.point.ok,
    inference.ok.WEE = WEE.se.ok,
    point.time.WEE = WEE.fit$time,
    se.time.WEE = boot$WEE$time,
    nvalid.boot.WEE = boot$WEE$nvalid,
    projection.WEE = projection.WEE,
    error.WEE = WEE.fit$error,

    beta1.OracleWEE = if (Oracle.point.ok) Oracle.fit$value$beta[1] else NA_real_,
    beta2.OracleWEE = if (Oracle.point.ok) Oracle.fit$value$beta[2] else NA_real_,
    se1.OracleWEE = if (Oracle.se.ok) boot$OracleWEE$se[1] else NA_real_,
    se2.OracleWEE = if (Oracle.se.ok) boot$OracleWEE$se[2] else NA_real_,
    point.ok.OracleWEE = Oracle.point.ok,
    inference.ok.OracleWEE = Oracle.se.ok,
    point.time.OracleWEE = Oracle.fit$time,
    se.time.OracleWEE = boot$OracleWEE$time,
    nvalid.boot.OracleWEE = boot$OracleWEE$nvalid,
    projection.OracleWEE = projection.Oracle,
    error.OracleWEE = Oracle.fit$error,

    beta1.EM = if (EM.point.ok) EM.fit$value$beta[1] else NA_real_,
    beta2.EM = if (EM.point.ok) EM.fit$value$beta[2] else NA_real_,
    se1.EM = if (EM.se.ok) EM.se.beta[1] else NA_real_,
    se2.EM = if (EM.se.ok) EM.se.beta[2] else NA_real_,
    point.ok.EM = EM.point.ok,
    inference.ok.EM = EM.se.ok,
    converged.EM = EM.converged,
    iterations.EM = if (isTRUE(EM.fit$ok)) EM.fit$value$iterations else NA_integer_,
    point.time.EM = EM.fit$time,
    se.time.EM = boot$EM$time,
    nvalid.boot.EM = boot$EM$nvalid,
    projection.EM = FALSE,
    error.EM = EM.fit$error,

    curve.WEE = curve.WEE,
    curve.OracleWEE = curve.Oracle,
    curve.EM = curve.EM,
    incidence.WEE = incidence.WEE,
    incidence.OracleWEE = incidence.Oracle,
    incidence.EM = incidence.EM
  )
}

failed_replication_result <- function(rep.id, error.message) {
  out <- list(
    rep = rep.id,
    fatal.error = as.character(error.message),
    censor.rate = NA_real_,
    cure.rate = NA_real_
  )
  for (method in methods) {
    out[[paste0("beta1.", method)]] <- NA_real_
    out[[paste0("beta2.", method)]] <- NA_real_
    out[[paste0("se1.", method)]] <- NA_real_
    out[[paste0("se2.", method)]] <- NA_real_
    out[[paste0("point.ok.", method)]] <- FALSE
    out[[paste0("inference.ok.", method)]] <- FALSE
    out[[paste0("point.time.", method)]] <- NA_real_
    out[[paste0("se.time.", method)]] <- NA_real_
    out[[paste0("nvalid.boot.", method)]] <- NA_integer_
    out[[paste0("projection.", method)]] <- NA
    out[[paste0("error.", method)]] <- as.character(error.message)
    out[[paste0("curve.", method)]] <- rep(NA_real_, length(time.grid))
    out[[paste0("incidence.", method)]] <- rep(NA_real_, nrow(Z.eval))
  }
  out$converged.EM <- FALSE
  out$iterations.EM <- NA_integer_
  out
}

# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------
methods <- c("WEE", "OracleWEE", "EM")

scalar_result_frame <- function(results) {
  scalar.names <- names(results[[1]])[
    !vapply(results[[1]], function(x) length(x) > 1L, logical(1))
  ]
  rows <- lapply(results, function(x) {
    as.data.frame(x[scalar.names], stringsAsFactors = FALSE, check.names = FALSE)
  })
  do.call(rbind, rows)
}

parameter_summary <- function(raw) {
  z <- qnorm(0.975)
  ans <- list()
  index <- 1L
  for (method in methods) {
    point.ok <- as.logical(raw[[paste0("point.ok.", method)]])
    inference.ok <- as.logical(raw[[paste0("inference.ok.", method)]])
    point.time <- as.numeric(raw[[paste0("point.time.", method)]])
    for (j in 1:2) {
      est <- as.numeric(raw[[paste0("beta", j, ".", method)]])
      se <- as.numeric(raw[[paste0("se", j, ".", method)]])
      ok.point <- point.ok & is.finite(est)
      ok.inf <- inference.ok & is.finite(est) & is.finite(se)
      truth <- beta.true[j]
      ans[[index]] <- data.frame(
        Setting = SETTING_ID,
        Method = method,
        Parameter = paste0("beta", j),
        Truth = truth,
        Bias = if (any(ok.point)) mean(est[ok.point]) - truth else NA_real_,
        SD = if (sum(ok.point) >= 2L) sd(est[ok.point]) else NA_real_,
        SE = if (any(ok.inf)) mean(se[ok.inf]) else NA_real_,
        CP = if (any(ok.inf)) {
          mean(est[ok.inf] - z * se[ok.inf] <= truth & truth <= est[ok.inf] + z * se[ok.inf])
        } else NA_real_,
        MSE = if (any(ok.point)) mean((est[ok.point] - truth)^2) else NA_real_,
        MeanPointTimeSec = if (any(is.finite(point.time))) {
          mean(point.time[is.finite(point.time)])
        } else NA_real_,
        `PointValidN/n` = sprintf("%d/%d", sum(ok.point), repetition),
        `InferenceValidN/n` = sprintf("%d/%d", sum(ok.inf), repetition),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      index <- index + 1L
    }
  }
  do.call(rbind, ans)
}

curve_matrix <- function(results, method) {
  do.call(rbind, lapply(results, function(x) x[[paste0("curve.", method)]]))
}

incidence.truth <- as.numeric(pai(drop(cbind(1, Z.eval) %*% gamma.true)))

incidence_matrix <- function(results, method) {
  do.call(rbind, lapply(results, function(x) x[[paste0("incidence.", method)]]))
}

incidence_summary <- function(results) {
  truth <- data.frame(
    Setting = SETTING_ID,
    Method = "Truth",
    Z1 = Z.eval[, 1],
    Z2 = Z.eval[, 2],
    Truth = incidence.truth,
    Mean = incidence.truth,
    Q025 = incidence.truth,
    Q975 = incidence.truth,
    Bias = 0,
    RMSE = 0,
    ValidN = length(results),
    stringsAsFactors = FALSE
  )
  out <- list(truth)
  for (method in methods) {
    mat <- incidence_matrix(results, method)
    mean.hat <- apply(mat, 2, function(x) {
      x <- x[is.finite(x)]
      if (length(x)) mean(x) else NA_real_
    })
    out[[length(out) + 1L]] <- data.frame(
      Setting = SETTING_ID,
      Method = method,
      Z1 = Z.eval[, 1],
      Z2 = Z.eval[, 2],
      Truth = incidence.truth,
      Mean = mean.hat,
      Q025 = apply(mat, 2, function(x) {
        x <- x[is.finite(x)]
        if (length(x)) unname(quantile(x, 0.025)) else NA_real_
      }),
      Q975 = apply(mat, 2, function(x) {
        x <- x[is.finite(x)]
        if (length(x)) unname(quantile(x, 0.975)) else NA_real_
      }),
      Bias = mean.hat - incidence.truth,
      RMSE = vapply(seq_len(ncol(mat)), function(j) {
        x <- mat[, j]
        ok <- is.finite(x)
        if (any(ok)) sqrt(mean((x[ok] - incidence.truth[j])^2)) else NA_real_
      }, numeric(1)),
      ValidN = apply(mat, 2, function(x) sum(is.finite(x))),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

curve_summary <- function(results) {
  out <- list(data.frame(
    Setting = SETTING_ID,
    Method = "Truth",
    Time = time.grid,
    Mean = truth.curve,
    Q025 = truth.curve,
    Q975 = truth.curve,
    ValidN = length(results),
    stringsAsFactors = FALSE
  ))
  for (method in methods) {
    mat <- curve_matrix(results, method)
    out[[length(out) + 1L]] <- data.frame(
      Setting = SETTING_ID,
      Method = method,
      Time = time.grid,
      Mean = apply(mat, 2, function(x) {
        x <- x[is.finite(x)]
        if (length(x)) mean(x) else NA_real_
      }),
      Q025 = apply(mat, 2, function(x) {
        x <- x[is.finite(x)]
        if (length(x)) unname(quantile(x, 0.025)) else NA_real_
      }),
      Q975 = apply(mat, 2, function(x) {
        x <- x[is.finite(x)]
        if (length(x)) unname(quantile(x, 0.975)) else NA_real_
      }),
      ValidN = apply(mat, 2, function(x) sum(is.finite(x))),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

survival_metrics <- function(results) {
  do.call(rbind, lapply(methods, function(method) {
    mat <- curve_matrix(results, method)
    ise <- apply(mat, 1, function(curve) {
      if (any(!is.finite(curve))) return(NA_real_)
      trapz(time.grid, (curve - truth.curve)^2)
    })
    data.frame(
      Setting = SETTING_ID,
      Method = method,
      MeanISE = mean(ise, na.rm = TRUE),
      MedianISE = median(ise, na.rm = TRUE),
      SD_ISE = sd(ise, na.rm = TRUE),
      ValidN = sum(is.finite(ise)),
      stringsAsFactors = FALSE
    )
  }))
}

plot_incidence_surface <- function(d, filename) {
  required <- c("Setting", "Method", "Z1", "Z2", "Mean")
  missing.columns <- setdiff(required, names(d))
  if (length(missing.columns)) {
    stop("Missing column(s): ", paste(missing.columns, collapse = ", "))
  }

  surface.methods <- c("Truth", "WEE", "OracleWEE", "EM")
  missing.methods <- setdiff(surface.methods, unique(d$Method))
  if (length(missing.methods)) {
    stop("Missing method(s): ", paste(missing.methods, collapse = ", "))
  }

  z1 <- sort(unique(d$Z1))
  z2 <- sort(unique(d$Z2))
  if (length(z1) < 2L || length(z2) < 2L) {
    stop("The incidence results must contain a two-dimensional Z grid.")
  }

  surface_matrix <- function(method, value) {
    dm <- d[d$Method == method, c("Z1", "Z2", value), drop = FALSE]
    ans <- matrix(NA_real_, nrow = length(z1), ncol = length(z2))
    ans[cbind(match(dm$Z1, z1), match(dm$Z2, z2))] <- dm[[value]]
    ans
  }

  method.label <- c(
    Truth = "Truth",
    WEE = "WEE",
    OracleWEE = "Oracle WEE",
    EM = "EM"
  )
  method.col <- c(
    Truth = "#BDBDBD",
    WEE = "#56B4E9",
    OracleWEE = "#009E73",
    EM = "#E69F00"
  )

  png(filename, width = 2200, height = 1800, res = 220)
  on.exit(dev.off(), add = TRUE)
  par(
    mfrow = c(2, 2),
    mar = c(4.8, 4.8, 3.2, 2.2),
    oma = c(0, 0, 0, 0),
    mgp = c(4.2, 1.0, 0)
  )
  for (method in surface.methods) {
    zz <- surface_matrix(method, "Mean")
    if (!any(is.finite(zz))) {
      plot.new()
      title(main = method.label[method])
      text(0.5, 0.5, "No valid estimates")
    } else {
      persp(
        x = z1,
        y = z2,
        z = zz,
        zlim = c(0, 1),
        theta = -35,
        phi = 25,
        expand = 0.72,
        ticktype = "detailed",
        col = adjustcolor(method.col[method], alpha.f = 0.82),
        border = adjustcolor("grey25", alpha.f = 0.35),
        xlab = expression(Z[1]),
        ylab = expression(Z[2]),
        zlab = expression(pi(z[1], z[2])),
        cex.lab = 1.15,
        cex.axis = 0.90,
        main = method.label[method]
      )
    }
  }
}

# -----------------------------------------------------------------------------
# Run simulation
# -----------------------------------------------------------------------------
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("SETTING_ID:", SETTING_ID, "\n")
cat("Replications:", REP_START, "to", REP_END, "\n")
cat("Bootstrap B:", BOOT_B, "; minimum valid:", MIN_VALID_BOOT, "\n")
cat("EM degree:", EM_DEGREE, "; EM tolerance:", EM_TOL, "\n")
cat("Parallel backend:", PARALLEL_BACKEND, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")

cl <- NULL
results <- NULL

if (PARALLEL_BACKEND == "mpi") {
  cl <- startMPIcluster()
  registerDoMPI(cl)
  cat("MPI workers:", foreach::getDoParWorkers(), "\n")
  export.names <- setdiff(ls(envir = .GlobalEnv), c("cl", "results"))
  results <- foreach(
    rep.id = rep.ids,
    .packages = "survival",
    .export = export.names,
    .errorhandling = "pass",
    .inorder = TRUE
  ) %dopar% {
    tryCatch(
      one_replication(rep.id),
      error = function(e) failed_replication_result(rep.id, conditionMessage(e))
    )
  }
} else {
  results <- lapply(rep.ids, function(rep.id) {
    tryCatch(
      one_replication(rep.id),
      error = function(e) failed_replication_result(rep.id, conditionMessage(e))
    )
  })
}

tryCatch({
  for (i in seq_along(results)) {
    if (inherits(results[[i]], "error")) {
      results[[i]] <- failed_replication_result(rep.ids[i], conditionMessage(results[[i]]))
    }
  }
  if (!length(results)) stop("No replication was returned.")

  raw <- scalar_result_frame(results)
  summary.table <- parameter_summary(raw)
  curve.table <- curve_summary(results)
  survival.table <- survival_metrics(results)
  incidence.table <- incidence_summary(results)

  prefix <- file.path(
    OUTPUT_DIR,
    sprintf("%s_rep%d-%d", SETTING_ID, REP_START, REP_END)
  )
  write.csv(raw, paste0(prefix, "_raw.csv"), row.names = FALSE)
  write.csv(summary.table, paste0(prefix, "_summary.csv"), row.names = FALSE)
  write.csv(curve.table, paste0(prefix, "_survival.csv"), row.names = FALSE)
  write.csv(survival.table, paste0(prefix, "_survival_metrics.csv"), row.names = FALSE)
  write.csv(incidence.table, paste0(prefix, "_incidence.csv"), row.names = FALSE)
  plot_incidence_surface(
    incidence.table,
    paste0(prefix, "_incidence_surface.png")
  )

  cat("\n================ Parameter performance ================\n")
  print(transform(
    summary.table,
    Bias = round(Bias, 4), SD = round(SD, 4), SE = round(SE, 4),
    CP = round(CP, 3), MSE = round(MSE, 4),
    MeanPointTimeSec = round(MeanPointTimeSec, 4)
  ), row.names = FALSE)
  cat("\n================ Survival-function ISE ================\n")
  print(survival.table, row.names = FALSE)
  cat("\nTotal wall-clock time (minutes):", round(as.numeric(
    difftime(Sys.time(), timestart, units = "mins")
  ), 3), "\n")
}, finally = {
  if (!is.null(cl)) {
    closeCluster(cl)
    Rmpi::mpi.quit()
  }
})
