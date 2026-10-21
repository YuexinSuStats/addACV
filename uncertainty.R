library(iterators)
library(foreach)
library(survival)
library(ggplot2)

# =============================================================================
# Comparison of uncertainty-adjusted and naive ICV estimators
#
# Design inherited from addICV.R:
#   * pai1: logistic cure-incidence link;
#   * n = 300;
#   * low/high censoring selected by SETTING_ID;
#   * two continuous covariates;
#   * bandwidth h_j = 3 * n^(-2/7) * sd(Z_j);
#   * the same subtf() rule as addICV.R;
#   * four subgroups defined by Z1 = 1 and Z2 = 2 cutoffs;
#   * two auxiliary time points.
#
# Methods:
#   WEE       : weighted estimating-equation estimator;
#   ICV       : uncertainty-adjusted ICV, using
#               V22.U = V22 + (n/N) D.pi Sigma.phi D.pi;
#   ICV_naive : ICV that treats the external KM estimates as fixed.
#
# A fresh external sample is generated in every Monte Carlo replication. This
# is essential because the Monte Carlo SD must contain external-sample
# variation when assessing whether the uncertainty-adjusted SE is calibrated.
# Within a replication, the five external samples are nested subsets of one
# randomly permuted sample of size max(N), providing common random numbers for
# a smoother and more efficient comparison across N/n.
# =============================================================================

timestart <- Sys.time()

# -----------------------------------------------------------------------------
# Runtime settings
# -----------------------------------------------------------------------------
SETTING_ID <- Sys.getenv("SETTING_ID", "S1-L-n300")
REP_START <- as.integer(Sys.getenv("REP_START", "1"))
REP_END <- as.integer(Sys.getenv("REP_END", "500"))
OUTPUT_DIR <- Sys.getenv(
  "OUTPUT_DIR",
  file.path("uncertainty", SETTING_ID)
)
PARALLEL_SEED_BASE <- as.integer(
  Sys.getenv("PARALLEL_SEED_BASE", "20260825")
)
PARALLEL_BACKEND <- tolower(Sys.getenv("PARALLEL_BACKEND", "mpi"))

valid.settings <- c("S1-L-n300", "S1-H-n300")
if (!SETTING_ID %in% valid.settings) {
  stop("SETTING_ID must be S1-L-n300 or S1-H-n300.")
}
if (is.na(REP_START) || is.na(REP_END) ||
    REP_START < 1L || REP_END < REP_START) {
  stop("REP_START and REP_END must satisfy 1 <= REP_START <= REP_END.")
}
if (is.na(PARALLEL_SEED_BASE)) {
  stop("PARALLEL_SEED_BASE must be an integer.")
}
if (!PARALLEL_BACKEND %in% c("mpi", "serial")) {
  stop("PARALLEL_BACKEND must be 'mpi' or 'serial'.")
}

setting.parts <- strsplit(SETTING_ID, "-", fixed = TRUE)[[1]]
pai.setting <- setting.parts[1]
censor.setting <- setting.parts[2]
n <- as.integer(sub("^n", "", setting.parts[3]))
rep.ids <- seq.int(REP_START, REP_END)
repetition <- length(rep.ids)

# The default grid puts more resolution where the external uncertainty changes
# most rapidly. It can be overridden, but exactly five ratios are required.
ratio.text <- Sys.getenv("N_RATIOS", "3,5,10,20,50")
N.ratios.requested <- suppressWarnings(as.numeric(
  strsplit(gsub("[[:space:]]", "", ratio.text), ",", fixed = TRUE)[[1]]
))
if (length(N.ratios.requested) != 5L ||
    any(!is.finite(N.ratios.requested)) ||
    any(N.ratios.requested <= 0)) {
  stop("N_RATIOS must contain exactly five positive numbers, e.g. 3,5,10,20,50.")
}
N.grid <- as.integer(round(n * N.ratios.requested))
if (anyDuplicated(N.grid)) {
  stop("N_RATIOS gives duplicated external sample sizes after rounding.")
}
N.ratios <- N.grid / n
N.max <- max(N.grid)

N.ref <- as.integer(Sys.getenv("N_REFERENCE", "20000"))
B <- as.integer(Sys.getenv("NBOOT", as.character(n * 4L)))
MIN_VALID_BOOT <- as.integer(Sys.getenv("MIN_VALID_BOOT", "30"))
COV_EIGEN_FLOOR <- as.numeric(Sys.getenv("COV_EIGEN_FLOOR", "1e-8"))

if (is.na(N.ref) || N.ref < 100L) {
  stop("N_REFERENCE must be an integer of at least 100.")
}
if (is.na(B) || B < 2L) {
  stop("NBOOT must be an integer greater than one.")
}
if (is.na(MIN_VALID_BOOT) || MIN_VALID_BOOT < 2L ||
    MIN_VALID_BOOT > B) {
  stop("MIN_VALID_BOOT must satisfy 2 <= MIN_VALID_BOOT <= NBOOT.")
}
if (!is.finite(COV_EIGEN_FLOOR) || COV_EIGEN_FLOOR <= 0) {
  stop("COV_EIGEN_FLOOR must be positive.")
}

# Independent random-number streams.
SEED_TSTAR_REFERENCE <- 100000L
SEED_INTERNAL_DATA <- 200000L
SEED_EXTERNAL_DATA <- 300000L

# -----------------------------------------------------------------------------
# Simulation model: identical to pai1 in addICV.R
# -----------------------------------------------------------------------------
beta.true <- c(1, 1)
k.num <- 4L
m.num <- 2L
AUX_NAME <- "AUX2"
method.order <- c("WEE", "ICV", "ICV_naive")

if (censor.setting == "L") {
  lambdaC <- 0.50
  gamma.true <- c(-1.5, 1, 1)
} else {
  lambdaC <- 1.00
  gamma.true <- c(-2.5, 1, 1)
}

pai <- function(u) stats::plogis(u)

subgf <- function(Z) {
  rbind(
    Z[, 1] <  1 & Z[, 2] <  2,
    Z[, 1] >= 1 & Z[, 2] <  2,
    Z[, 1] <  1 & Z[, 2] >= 2,
    Z[, 1] >= 1 & Z[, 2] >= 2
  )
}

tau0 <- 2
kappa <- 1

H_0 <- function(t) {
  ans <- rep(Inf, length(t))
  id <- t < tau0
  ans[id] <- 2 * log1p(t[id]) - kappa * log1p(-t[id] / tau0)
  ans
}

# -----------------------------------------------------------------------------
# Data generation
# -----------------------------------------------------------------------------
Get.random <- function(sample.size) {
  Z <- matrix(1, nrow = sample.size, ncol = length(beta.true))
  Z[, 1] <- stats::runif(sample.size, 0, 2)
  Z[, 2] <- stats::runif(sample.size, 0, 4)

  T0 <- numeric(sample.size)
  u <- stats::runif(sample.size, 0, 1)
  for (i in seq_len(sample.size)) {
    S0 <- function(t) {
      if (t >= tau0) return(0)
      exp(-H_0(t) - as.numeric(Z[i, ] %*% beta.true) * t)
    }
    T0[i] <- stats::uniroot(
      function(t) 1 - S0(t) - u[i],
      interval = c(0, tau0)
    )$root
  }

  Z.multi <- cbind(1, Z)
  peta <- pai(Z.multi %*% gamma.true)
  eta <- stats::rbinom(sample.size, 1, peta)
  TT <- ifelse(eta == 1, T0, Inf)

  C <- stats::rexp(sample.size, lambdaC)
  Y <- pmin(TT, C)
  delta <- as.numeric(eta == 1 & T0 <= C)

  ord <- order(Y)
  data <- data.frame(Y = Y, delta = delta, Z = Z)[ord, , drop = FALSE]
  eta.sorted <- eta[ord]

  list(
    data = data,
    eta = eta.sorted,
    censor.rate = mean(delta == 0),
    cure.rate = mean(eta == 0)
  )
}

Get.Data <- function(rep.ids) {
  data <- vector("list", length(rep.ids))
  censor.rate <- cure.rate <- numeric(length(rep.ids))

  for (i in seq_along(rep.ids)) {
    set.seed(SEED_INTERNAL_DATA + rep.ids[i])
    tmp <- Get.random(n)
    data[[i]] <- tmp$data
    censor.rate[i] <- tmp$censor.rate
    cure.rate[i] <- tmp$cure.rate
  }

  list(data = data, censor.rate = censor.rate, cure.rate = cure.rate)
}

# -----------------------------------------------------------------------------
# Auxiliary time points and subgroup KM estimates
# -----------------------------------------------------------------------------
# Kept identical to addICV.R, apart from boundary checks.
subtf <- function(data, m.num) {
  data <- data[order(data[, 1]), , drop = FALSE]
  Y <- data$Y
  failure.Y <- data$Y[data$delta == 1]
  if (length(failure.Y) == 0L) {
    stop("No observed failure time is available for selecting tstar.")
  }

  tau <- which.max(failure.Y)
  index <- if (m.num == 1L) {
    round(tau * 0.5)
  } else {
    round(seq_len(m.num) / (m.num + 1) * tau)
  }
  index <- pmin(length(Y), pmax(1L, index))
  as.numeric(Y[index])
}

getprob <- function(data, tstar, m.num) {
  pro <- matrix(NA_real_, k.num, m.num)
  group <- subgf(as.matrix(data[, -(1:2), drop = FALSE]))

  for (i in seq_len(k.num)) {
    gi <- which(group[i, ])
    if (length(gi) == 0L) next

    gd <- data[gi, , drop = FALSE]
    fit <- survival::survfit(
      survival::Surv(gd[, 1], gd[, 2]) ~ 1
    )
    pro[i, ] <- summary(
      fit,
      times = tstar,
      extend = TRUE
    )$surv
  }

  pro
}

# -----------------------------------------------------------------------------
# Nonparametric weights
# -----------------------------------------------------------------------------
KK <- function(x, hn) {
  3 / 4 * (1 - (x / hn)^2) * (abs(x / hn) <= 1)
}

NWweight2 <- function(Z, Z.sample, hn, lam = rep(1, nrow(Z.sample))) {
  m <- nrow(Z)
  p1 <- ncol(Z)
  nn <- nrow(Z.sample)
  Z <- as.data.frame(Z)
  Z.sample <- as.data.frame(Z.sample)

  MZ <- array(apply(Z, 2, function(x) rep(x, nn)), dim = c(m, nn, p1))
  tMZ <- aperm(
    array(apply(Z.sample, 2, function(x) rep(x, m)), dim = c(nn, m, p1)),
    c(2, 1, 3)
  )
  dMZ <- MZ - tMZ

  Kp <- array(NA_real_, dim = c(m, nn, p1))
  Kp[, , 1] <- KK(dMZ[, , 1], hn[1])
  Kp[, , 2] <- KK(dMZ[, , 2], hn[2])
  prodKp <- apply(Kp, c(1, 2), prod)

  Lam <- matrix(lam, nrow = m, ncol = nn, byrow = TRUE)
  lamKp <- prodKp * Lam
  rs <- matrix(rowSums(lamKp), nrow = m, ncol = nn, byrow = FALSE)
  rs[rs == 0] <- 1e-10
  lamKp / rs
}

S.Ti <- function(data, w) {
  m <- nrow(w)
  nn <- ncol(w)
  k <- cbind(rep(0, m), w[, seq_len(nn - 1L), drop = FALSE])
  cumsum.k <- t(apply(k, 1, cumsum))
  cum.w <- 1 - cumsum.k
  cum.w <- cum.w * (cum.w > w + 1e-5) +
    (w + 1e-4) * (cum.w <= w + 1e-5)
  q <- 1 - w / cum.w
  delta <- matrix(data$delta, nrow = m, ncol = nn, byrow = TRUE)
  q <- replace(q, delta == 0, 1)
  q <- t(apply(q, 1, cumprod))

  Y <- data$Y
  r <- rle(Y)
  ei <- cumsum(r$lengths)
  li <- rep(ei, r$lengths)
  q[cbind(seq_len(m), li[seq_len(m)])]
}

P.T <- function(data, w) {
  m <- nrow(w)
  nn <- ncol(w)
  k <- cbind(rep(0, m), w[, seq_len(nn - 1L), drop = FALSE])
  cumsum.k <- t(apply(k, 1, cumsum))
  cum.w <- 1 - cumsum.k
  w <- replace(w, cum.w == 0, 0)
  cum.w <- replace(cum.w, cum.w == 0, 1)
  q <- 1 - w / cum.w
  delta <- matrix(data$delta, nrow = m, ncol = nn, byrow = TRUE)
  q <- replace(q, delta == 0, 1)
  q <- apply(q, 1, prod)
  as.numeric(replace(q, q <= 0, 0))
}

# -----------------------------------------------------------------------------
# WEE and the internal component of the control variable
# -----------------------------------------------------------------------------
get.est <- function(data, aux.info) {
  data <- data[order(data[, 1]), , drop = FALSE]
  Y <- data$Y
  delta <- data$delta
  Z <- as.matrix(data[, -(1:2), drop = FALSE])
  p <- ncol(Z)
  nn <- nrow(data)

  if (!any(delta == 1)) {
    return(list(ok = FALSE, reason = "no observed failure"))
  }
  T1.max <- max(Y[delta == 1])

  # Bandwidth required by addICV.R.
  h <- 3 * nn^(-2 / 7) * apply(Z, 2L, stats::sd)
  if (any(!is.finite(h)) || any(h <= 0)) {
    return(list(ok = FALSE, reason = "invalid bandwidth"))
  }

  wght <- NWweight2(Z = Z, Z.sample = Z, hn = h)
  S.beran <- S.Ti(data, wght)
  P.T.inf <- P.T(data, wght)
  S.safe <- S.beran
  S.safe[S.safe == 0] <- 1e-8
  case.wght <- delta +
    (1 - delta) * (S.beran - P.T.inf) / S.safe * (Y <= T1.max)

  dif <- diff(c(0, Y))
  rle.Y <- rle(Y)
  start.index <- c(1L, 1L + head(cumsum(rle.Y$lengths), -1L))
  first.indices <- rep(start.index, rle.Y$lengths)

  Zbarr <- matrix(0, nn, p)
  for (j in seq_len(nn)) {
    index <- first.indices[j]:nn
    wj <- case.wght[index]
    sw <- sum(wj)
    if (sw > 0) {
      Zbarr[j, ] <- colSums(Z[index, , drop = FALSE] * wj) / sw
    }
  }

  fenzi <- colSums(case.wght * delta * (Z - Zbarr))
  fenmu <- matrix(0, p, p)
  for (j in seq_len(nn)) {
    if (dif[j] > 1e-8) {
      index <- first.indices[j]:nn
      wj <- case.wght[index]
      Zc <- sweep(Z[index, , drop = FALSE], 2, Zbarr[j, ], "-")
      fenmu <- fenmu +
        (t(Zc) %*% sweep(Zc, 1, wj, "*")) * dif[j]
    }
  }

  beta.WEE <- tryCatch(
    as.vector(fenzi %*% solve(fenmu)),
    error = function(e) rep(NA_real_, length(beta.true))
  )
  if (any(!is.finite(beta.WEE))) {
    return(list(ok = FALSE, reason = "non-finite beta"))
  }

  lambest <- function(ts) {
    eps <- 1e-12
    risk.w <- rev(cumsum(rev(case.wght)))
    zbeta <- as.vector(Z %*% beta.WEE)
    risk.wz <- rev(cumsum(rev(case.wght * zbeta)))

    event.num <- delta * case.wght * (Y <= ts)
    risk.event <- risk.w[first.indices]
    if (any(event.num != 0 & risk.event <= eps)) {
      stop("Zero risk-set denominator at an event time.")
    }

    term1 <- numeric(nn)
    valid.event <- risk.event > eps
    term1[valid.event] <-
      event.num[valid.event] / risk.event[valid.event]
    r1 <- sum(term1)
    point <- findInterval(ts, Y)

    if (point == 0L) {
      if (risk.w[1] <= eps) stop("Zero initial risk-set denominator.")
      r2 <- ts * risk.wz[1] / risk.w[1]
    } else {
      id <- seq_len(point)
      if (any(dif[id] > 0 & risk.w[id] <= eps)) {
        stop("Zero risk-set denominator in baseline-hazard integral.")
      }
      valid.interval <- dif[id] > 0 & risk.w[id] > eps
      r2 <- sum(
        dif[id][valid.interval] *
          risk.wz[id][valid.interval] /
          risk.w[id][valid.interval]
      )
      if (point < nn && ts > Y[point]) {
        if (risk.w[point + 1L] <= eps) {
          stop("Zero risk-set denominator after the last observed time.")
        }
        r2 <- r2 +
          (ts - Y[point]) *
          risk.wz[point + 1L] /
          risk.w[point + 1L]
      }
    }

    ans <- r1 - r2
    if (!is.finite(ans)) stop("Non-finite baseline cumulative hazard.")
    ans
  }

  subdata <- cbind(1 - P.T.inf, Z)
  group <- subgf(Z)
  pai.hat <- rowMeans(group)

  phi.by.method <- vector("list", length(aux.info))
  names(phi.by.method) <- names(aux.info)

  for (method in names(aux.info)) {
    tstar <- aux.info[[method]]$tstar
    tpro <- aux.info[[method]]$tpro
    current.m <- length(tstar)

    alpha <- tryCatch(
      vapply(tstar, lambest, numeric(1)),
      error = function(e) NULL
    )
    if (is.null(alpha) || any(!is.finite(alpha))) {
      phi.by.method[[method]] <- rep(NA_real_, k.num * current.m)
      next
    }

    phi <- matrix(0, k.num, current.m)
    for (m in seq_len(current.m)) {
      for (k in seq_len(k.num)) {
        gi <- which(group[k, ])
        if (length(gi) > 0L) {
          gd <- subdata[gi, , drop = FALSE]
          res <- numeric(nrow(gd))
          for (i in seq_len(nrow(gd))) {
            res[i] <- (1 - gd[i, 1]) +
              gd[i, 1] * exp(
                -alpha[m] -
                  as.numeric(beta.WEE %*% gd[i, -1]) * tstar[m]
              ) - tpro[k, m]
          }
          # Division by n, rather than subgroup size, is retained from addICV.
          phi[k, m] <- sum(res) / nn
        }
      }
    }

    phi.by.method[[method]] <- as.vector(phi)
  }

  list(
    ok = TRUE,
    est.beta = beta.WEE,
    est.phi = phi.by.method,
    pai = pai.hat
  )
}

# -----------------------------------------------------------------------------
# Bootstrap covariance components
# -----------------------------------------------------------------------------
# aux.zero uses the fixed tstar and tpro = 0. For any external KM vector p,
# phi(p) is recovered exactly as phi(0) - diag(replicated pai) p. Therefore the
# expensive WEE/bootstrap calculation is performed once per internal sample,
# not five times.
boot.cov <- function(data, aux.zero, B) {
  nn <- nrow(data)
  p <- length(beta.true)
  q <- k.num * m.num

  beta.boot <- matrix(NA_real_, B, p)
  phi0.boot <- matrix(NA_real_, B, q)
  pai.boot <- matrix(NA_real_, B, k.num)
  km.boot <- matrix(NA_real_, B, q)

  for (b in seq_len(B)) {
    bs <- data[sample.int(nn, nn, replace = TRUE), , drop = FALSE]
    fit <- tryCatch(
      get.est(bs, aux.zero),
      error = function(e) list(ok = FALSE)
    )
    if (!isTRUE(fit$ok)) next

    beta.boot[b, ] <- fit$est.beta
    phi0.boot[b, ] <- fit$est.phi[[AUX_NAME]]
    pai.boot[b, ] <- fit$pai
    km.boot[b, ] <- tryCatch(
      as.vector(getprob(bs, aux.zero[[AUX_NAME]]$tstar, m.num)),
      error = function(e) rep(NA_real_, q)
    )
  }

  ok.beta <- complete.cases(beta.boot)
  ok.joint <- ok.beta &
    complete.cases(phi0.boot) &
    complete.cases(pai.boot) &
    complete.cases(km.boot)

  list(
    beta = beta.boot,
    phi0 = phi0.boot,
    pai = pai.boot,
    km = km.boot,
    ok.beta = ok.beta,
    ok.joint = ok.joint,
    nvalid.beta = sum(ok.beta),
    nvalid.joint = sum(ok.joint)
  )
}

# -----------------------------------------------------------------------------
# Matrix utilities and ICV calculation
# -----------------------------------------------------------------------------
symmetrize <- function(A) (A + t(A)) / 2

regularized_inverse <- function(A, rel.floor = COV_EIGEN_FLOOR) {
  A <- symmetrize(A)
  ee <- eigen(A, symmetric = TRUE)
  scale <- max(abs(ee$values), 1)
  floor.value <- rel.floor * scale
  values.used <- pmax(ee$values, floor.value)

  A.used <- ee$vectors %*%
    diag(values.used, nrow = length(values.used)) %*%
    t(ee$vectors)
  A.inv <- ee$vectors %*%
    diag(1 / values.used, nrow = length(values.used)) %*%
    t(ee$vectors)

  list(
    matrix = symmetrize(A.used),
    inverse = symmetrize(A.inv),
    min.eigen = min(ee$values),
    floor = floor.value,
    regularized = any(ee$values < floor.value)
  )
}

cv.result <- function(beta.WEE, psi.est, V11, V12, V.aux, nn) {
  failed <- list(
    beta = rep(NA_real_, length(beta.true)),
    se = rep(NA_real_, length(beta.true)),
    valid = FALSE,
    min.eigen = NA_real_,
    regularized = NA
  )

  if (any(!is.finite(V11)) || any(!is.finite(V12)) ||
      any(!is.finite(V.aux)) || any(!is.finite(psi.est))) {
    return(failed)
  }

  inv <- tryCatch(
    regularized_inverse(V.aux),
    error = function(e) NULL
  )
  if (is.null(inv)) return(failed)

  beta.candidate <- as.vector(
    beta.WEE - t(V12) %*% inv$inverse %*% psi.est
  )
  var.candidate <- symmetrize(
    (V11 - t(V12) %*% inv$inverse %*% V12) / nn
  )

  valid <- all(is.finite(beta.candidate)) &&
    all(is.finite(var.candidate)) &&
    all(diag(var.candidate) >= -1e-10)
  if (!valid) {
    failed$min.eigen <- inv$min.eigen
    failed$regularized <- inv$regularized
    return(failed)
  }

  list(
    beta = beta.candidate,
    se = sqrt(pmax(diag(var.candidate), 0)),
    valid = TRUE,
    min.eigen = inv$min.eigen,
    regularized = inv$regularized
  )
}

method.rows <- function(
    rep.id, N.external, method, estimate, se, valid,
    cens.rate, cure.rate, external.cens.rate, external.cure.rate,
    external.min.group.size, nvalid.beta, nvalid.joint,
    regularized = NA, min.eigen = NA_real_,
    min.eigen.sigma = NA_real_) {
  data.frame(
    Rep = rep.id,
    Setting = SETTING_ID,
    InternalN = n,
    ExternalN = N.external,
    N_over_n = N.external / n,
    n_over_N = n / N.external,
    Method = method,
    Parameter = c("beta1", "beta2"),
    Estimate = estimate,
    EstSE = se,
    Valid = rep(valid, 2L),
    CensorRate = cens.rate,
    CureRate = cure.rate,
    ExternalCensorRate = external.cens.rate,
    ExternalCureRate = external.cure.rate,
    ExternalMinGroupSize = external.min.group.size,
    ValidBootWEE = nvalid.beta,
    ValidBootJoint = nvalid.joint,
    Regularized = rep(regularized, 2L),
    MinEigenAux = min.eigen,
    MinEigenSigmaPhi = min.eigen.sigma,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# One Monte Carlo replication
# -----------------------------------------------------------------------------
Repeat <- function(r, rep.id, aux.zero) {
  data <- Data$data[[r]]
  cens.rate <- Data$censor.rate[r]
  cure.rate <- Data$cure.rate[r]
  nn <- nrow(data)

  fit0 <- tryCatch(
    get.est(data, aux.zero),
    error = function(e) list(ok = FALSE)
  )
  if (!isTRUE(fit0$ok)) return(NULL)

  beta.WEE <- fit0$est.beta
  phi0.est <- fit0$est.phi[[AUX_NAME]]
  pai.vector <- rep(fit0$pai, times = m.num)

  set.seed(PARALLEL_SEED_BASE + rep.id)
  bc <- boot.cov(data, aux.zero, B)
  if (bc$nvalid.beta < MIN_VALID_BOOT) return(NULL)

  beta.boot.WEE <- bc$beta[bc$ok.beta, , drop = FALSE]
  V11.WEE <- stats::cov(beta.boot.WEE * sqrt(nn))
  if (any(!is.finite(V11.WEE))) return(NULL)
  se.WEE <- sqrt(pmax(diag(V11.WEE / nn), 0))

  joint.ok <- bc$nvalid.joint >= MIN_VALID_BOOT
  if (joint.ok) {
    ok <- bc$ok.joint
    beta.boot <- bc$beta[ok, , drop = FALSE]
    phi0.boot <- bc$phi0[ok, , drop = FALSE]
    pai.boot <- bc$pai[ok, , drop = FALSE]
    km.boot <- bc$km[ok, , drop = FALSE]

    V11 <- stats::cov(beta.boot * sqrt(nn))
    Sigma.phi <- symmetrize(stats::cov(km.boot * sqrt(nn)))
    min.eigen.sigma <- tryCatch(
      min(eigen(Sigma.phi, symmetric = TRUE, only.values = TRUE)$values),
      error = function(e) NA_real_
    )
    pai.boot.expanded <- pai.boot[
      , rep(seq_len(k.num), times = m.num), drop = FALSE
    ]
  } else {
    V11 <- Sigma.phi <- NULL
    min.eigen.sigma <- NA_real_
    pai.boot.expanded <- NULL
  }

  # One new external sample per MC replication; nested random subsets give the
  # five marginally correct external samples while pairing the N comparisons.
  set.seed(SEED_EXTERNAL_DATA + rep.id)
  external.full <- Get.random(N.max)
  random.order <- sample.int(N.max, N.max, replace = FALSE)
  external.data <- external.full$data[random.order, , drop = FALSE]
  external.eta <- external.full$eta[random.order]

  output <- vector("list", length(N.grid) * length(method.order))
  output.index <- 0L

  for (g in seq_along(N.grid)) {
    current.N <- N.grid[g]
    take <- seq_len(current.N)
    ext.data <- external.data[take, , drop = FALSE]
    ext.eta <- external.eta[take]
    ext.group <- subgf(as.matrix(ext.data[, -(1:2), drop = FALSE]))
    ext.censor <- mean(ext.data$delta == 0)
    ext.cure <- mean(ext.eta == 0)
    ext.min.group <- min(rowSums(ext.group))

    output.index <- output.index + 1L
    output[[output.index]] <- method.rows(
      rep.id = rep.id,
      N.external = current.N,
      method = "WEE",
      estimate = beta.WEE,
      se = se.WEE,
      valid = TRUE,
      cens.rate = cens.rate,
      cure.rate = cure.rate,
      external.cens.rate = ext.censor,
      external.cure.rate = ext.cure,
      external.min.group.size = ext.min.group,
      nvalid.beta = bc$nvalid.beta,
      nvalid.joint = bc$nvalid.joint
    )

    fit.adjusted <- fit.naive <- list(
      beta = rep(NA_real_, length(beta.true)),
      se = rep(NA_real_, length(beta.true)),
      valid = FALSE,
      min.eigen = NA_real_,
      regularized = NA
    )

    tpro <- tryCatch(
      getprob(ext.data, aux.zero[[AUX_NAME]]$tstar, m.num),
      error = function(e) matrix(NA_real_, k.num, m.num)
    )

    if (joint.ok && all(is.finite(tpro)) && all(is.finite(phi0.est))) {
      external.phi <- as.vector(tpro)
      psi.est <- phi0.est - pai.vector * external.phi
      psi.boot <- phi0.boot - sweep(
        pai.boot.expanded,
        2L,
        external.phi,
        "*"
      )

      V22 <- symmetrize(stats::cov(psi.boot * sqrt(nn)))
      V12 <- stats::cov(
        psi.boot * sqrt(nn),
        beta.boot * sqrt(nn)
      )
      D.pi <- diag(pai.vector, nrow = length(pai.vector))
      V22.U <- symmetrize(
        V22 + (nn / current.N) * D.pi %*% Sigma.phi %*% D.pi
      )

      # Adjusted ICV is the primary method and therefore receives the name ICV.
      fit.adjusted <- cv.result(
        beta.WEE, psi.est, V11, V12, V22.U, nn
      )
      fit.naive <- cv.result(
        beta.WEE, psi.est, V11, V12, V22, nn
      )
    }

    output.index <- output.index + 1L
    output[[output.index]] <- method.rows(
      rep.id = rep.id,
      N.external = current.N,
      method = "ICV",
      estimate = fit.adjusted$beta,
      se = fit.adjusted$se,
      valid = fit.adjusted$valid,
      cens.rate = cens.rate,
      cure.rate = cure.rate,
      external.cens.rate = ext.censor,
      external.cure.rate = ext.cure,
      external.min.group.size = ext.min.group,
      nvalid.beta = bc$nvalid.beta,
      nvalid.joint = bc$nvalid.joint,
      regularized = fit.adjusted$regularized,
      min.eigen = fit.adjusted$min.eigen,
      min.eigen.sigma = min.eigen.sigma
    )

    output.index <- output.index + 1L
    output[[output.index]] <- method.rows(
      rep.id = rep.id,
      N.external = current.N,
      method = "ICV_naive",
      estimate = fit.naive$beta,
      se = fit.naive$se,
      valid = fit.naive$valid,
      cens.rate = cens.rate,
      cure.rate = cure.rate,
      external.cens.rate = ext.censor,
      external.cure.rate = ext.cure,
      external.min.group.size = ext.min.group,
      nvalid.beta = bc$nvalid.beta,
      nvalid.joint = bc$nvalid.joint,
      regularized = fit.naive$regularized,
      min.eigen = fit.naive$min.eigen,
      min.eigen.sigma = min.eigen.sigma
    )
  }

  do.call(rbind, output)
}

# -----------------------------------------------------------------------------
# Summary and calibration plot
# -----------------------------------------------------------------------------
GetSummary <- function(raw.result, beta.true) {
  truth <- c(beta1 = beta.true[1], beta2 = beta.true[2])
  groups <- interaction(
    raw.result$ExternalN,
    factor(raw.result$Method, levels = method.order),
    factor(raw.result$Parameter, levels = names(truth)),
    drop = TRUE,
    lex.order = TRUE
  )

  split.data <- split(raw.result, groups)
  summary.rows <- lapply(split.data, function(d) {
    method <- as.character(d$Method[1])
    parameter <- as.character(d$Parameter[1])
    true.value <- truth[[parameter]]
    ok <- d$Valid & is.finite(d$Estimate) & is.finite(d$EstSE)
    est <- d$Estimate[ok]
    se <- d$EstSE[ok]

    if (length(est) >= 2L) {
      bias <- mean(est) - true.value
      mc.sd <- stats::sd(est)
      mean.se <- mean(se)
      z <- stats::qnorm(0.975)
      cp <- mean(est - z * se <= true.value & true.value <= est + z * se)
      mse <- mean((est - true.value)^2)
    } else {
      bias <- mc.sd <- mean.se <- cp <- mse <- NA_real_
    }

    reg <- d$Regularized[ok]
    reg.rate <- if (method == "WEE" || length(reg) == 0L || all(is.na(reg))) {
      NA_real_
    } else {
      mean(reg, na.rm = TRUE)
    }

    data.frame(
      Setting = SETTING_ID,
      InternalN = n,
      ExternalN = d$ExternalN[1],
      N_over_n = d$N_over_n[1],
      Method = method,
      Parameter = parameter,
      Bias = bias,
      SE = mean.se,
      SD = mc.sd,
      CP = cp,
      SE_SD = mean.se / mc.sd,
      MSE = mse,
      RMSE = sqrt(mse),
      ValidN = length(est),
      RegularizedRate = reg.rate,
      CensorRate = mean(d$CensorRate, na.rm = TRUE),
      CureRate = mean(d$CureRate, na.rm = TRUE),
      ExternalCensorRate = mean(d$ExternalCensorRate, na.rm = TRUE),
      ExternalCureRate = mean(d$ExternalCureRate, na.rm = TRUE),
      MeanValidBootWEE = mean(d$ValidBootWEE, na.rm = TRUE),
      MeanValidBootJoint = mean(d$ValidBootJoint, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  ans <- do.call(rbind, summary.rows)
  ans$Method <- factor(ans$Method, levels = method.order)
  ans$Parameter <- factor(ans$Parameter, levels = c("beta1", "beta2"))
  ans <- ans[order(ans$ExternalN, ans$Method, ans$Parameter), , drop = FALSE]
  ans$Method <- as.character(ans$Method)
  ans$Parameter <- as.character(ans$Parameter)
  rownames(ans) <- NULL
  ans
}

plot.calibration <- function(summary.result, file.base) {
  dat <- summary.result
  dat$SE_SD <- dat$SE / dat$SD
  dat <- subset(
    dat,
    Method %in% method.order &
      Parameter %in% c("beta1", "beta2") &
      is.finite(N_over_n) &
      is.finite(SE_SD)
  )
  if (nrow(dat) == 0L) {
    warning("No finite SE/SD values are available; calibration plot skipped.")
    return(invisible(NULL))
  }

  dat$Method <- factor(dat$Method, levels = method.order)
  dat$Parameter <- factor(
    dat$Parameter,
    levels = c("beta1", "beta2")
  )
  dat <- dat[
    order(dat$Parameter, dat$Method, dat$N_over_n),
    ,
    drop = FALSE
  ]

  x.breaks <- sort(unique(dat$N_over_n))
  method.labels <- c(
    WEE = "WEE",
    ICV = "UICV",
    ICV_naive = "Naive ICV"
  )
  method.colors <- c(
    WEE = "#000000",
    ICV = "#0072B2",
    ICV_naive = "#D55E00"
  )
  method.linetypes <- c(
    WEE = "dashed",
    ICV = "solid",
    ICV_naive = "dotted"
  )
  method.shapes <- c(
    WEE = 15,
    ICV = 16,
    ICV_naive = 17
  )

  p <- ggplot(
    dat,
    aes(
      x = N_over_n,
      y = SE_SD,
      color = Method,
      linetype = Method,
      shape = Method,
      group = Method
    )
  ) +
    geom_hline(
      yintercept = 1,
      color = "grey50",
      linewidth = 0.55,
      linetype = "dashed"
    ) +
    geom_line(linewidth = 0.85, na.rm = TRUE) +
    geom_point(size = 2.5, stroke = 0.7, na.rm = TRUE) +
    facet_wrap(~ Parameter, nrow = 1) +
    scale_x_log10(
      breaks = x.breaks,
      labels = format(x.breaks, trim = TRUE, scientific = FALSE)
    ) +
    scale_color_manual(
      values = method.colors,
      breaks = names(method.labels),
      labels = method.labels
    ) +
    scale_linetype_manual(
      values = method.linetypes,
      breaks = names(method.labels),
      labels = method.labels
    ) +
    scale_shape_manual(
      values = method.shapes,
      breaks = names(method.labels),
      labels = method.labels
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.12))) +
    labs(
      title = NULL,
      x = expression(N / n),
      y = "SE/SD",
      color = NULL,
      linetype = NULL,
      shape = NULL
    ) +
    guides(
      color = guide_legend(ncol = 1),
      linetype = guide_legend(ncol = 1),
      shape = guide_legend(ncol = 1)
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = c(0.985, 0.025),
      legend.justification = c(1, 0),
      legend.direction = "vertical",
      legend.background = element_rect(
        fill = "white",
        color = "grey40",
        linewidth = 0.3
      ),
      legend.key = element_blank(),
      legend.text = element_text(size = 9),
      plot.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.spacing = grid::unit(1, "lines")
    )

  ggsave(
    filename = paste0(file.base, ".pdf"),
    plot = p,
    width = 8.5,
    height = 4.8,
    units = "in"
  )
  ggsave(
    filename = paste0(file.base, ".png"),
    plot = p,
    width = 8.5,
    height = 4.8,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  invisible(p)
}

# -----------------------------------------------------------------------------
# Fixed auxiliary time points
# -----------------------------------------------------------------------------
# tstar is fixed before the Monte Carlo loop so that changing N affects only
# the precision of the external KM information, not the definition of the
# auxiliary estimand.
set.seed(SEED_TSTAR_REFERENCE)
reference.data <- Get.random(N.ref)$data
tstar <- subtf(reference.data, m.num)
rm(reference.data)

aux.zero <- setNames(
  list(list(tstar = tstar, tpro = matrix(0, k.num, m.num))),
  AUX_NAME
)

# -----------------------------------------------------------------------------
# Internal samples and parallel computation
# -----------------------------------------------------------------------------
Data <- Get.Data(rep.ids)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

tstar.file <- file.path(
  OUTPUT_DIR,
  sprintf("%s_tstar.csv", SETTING_ID)
)
utils::write.csv(
  data.frame(
    TimeIndex = seq_along(tstar),
    Tstar = tstar,
    ReferenceN = N.ref,
    stringsAsFactors = FALSE
  ),
  file = tstar.file,
  row.names = FALSE
)

grid.file <- file.path(
  OUTPUT_DIR,
  sprintf("%s_N_grid.csv", SETTING_ID)
)
utils::write.csv(
  data.frame(
    GridIndex = seq_along(N.grid),
    InternalN = n,
    ExternalN = N.grid,
    N_over_n = N.ratios,
    n_over_N = n / N.grid,
    stringsAsFactors = FALSE
  ),
  file = grid.file,
  row.names = FALSE
)

exports <- c(
  "Repeat", "get.est", "boot.cov", "getprob", "cv.result",
  "method.rows", "symmetrize", "regularized_inverse", "Get.random",
  "NWweight2", "S.Ti", "P.T", "KK", "subgf", "pai", "H_0",
  "Data", "aux.zero", "beta.true", "gamma.true", "lambdaC",
  "tau0", "kappa", "B", "n", "N.grid", "N.max", "m.num", "k.num",
  "AUX_NAME", "SETTING_ID", "method.order", "MIN_VALID_BOOT",
  "COV_EIGEN_FLOOR", "SEED_EXTERNAL_DATA", "PARALLEL_SEED_BASE",
  "rep.ids"
)

cluster.object <- NULL
if (PARALLEL_BACKEND == "mpi") {
  suppressPackageStartupMessages(library(Rmpi))
  suppressPackageStartupMessages(library(doMPI))
  cluster.object <- startMPIcluster()
  registerDoMPI(cluster.object)
} else {
  foreach::registerDoSEQ()
}

cat("Parallel backend:", PARALLEL_BACKEND, "\n")
cat("Workers:", foreach::getDoParWorkers(), "\n")
cat("SETTING_ID:", SETTING_ID, "\n")
cat("Replications:", REP_START, "to", REP_END, "\n")
cat("Internal n:", n, "\n")
cat("External N:", paste(N.grid, collapse = ", "), "\n")
cat("N/n ratios:", paste(N.ratios, collapse = ", "), "\n")
cat("Auxiliary time points:", paste(tstar, collapse = ", "), "\n")
cat("Bandwidth: 3 * n^(-2/7) * sample SD(Z_j)\n")
cat("Output directory:", OUTPUT_DIR, "\n")

tryCatch({
  MyResult <- foreach(
    r = seq_along(rep.ids),
    .packages = c("survival"),
    .export = exports,
    .errorhandling = "pass",
    .inorder = TRUE
  ) %dopar% {
    rep.id <- rep.ids[r]
    Repeat(r = r, rep.id = rep.id, aux.zero = aux.zero)
  }

  valid.output <- Filter(
    function(x) !is.null(x) && !inherits(x, "error") && is.data.frame(x),
    MyResult
  )
  if (length(valid.output) == 0L) {
    errors <- vapply(
      MyResult,
      function(x) {
        if (inherits(x, "error")) conditionMessage(x) else NA_character_
      },
      character(1)
    )
    stop(
      "All iterations failed. First worker error: ",
      errors[which(!is.na(errors))[1]]
    )
  }

  raw.result <- do.call(rbind, valid.output)
  rownames(raw.result) <- NULL
  summary.result <- GetSummary(raw.result, beta.true)

  file.prefix <- sprintf(
    "%s_rep%d-%d",
    SETTING_ID,
    REP_START,
    REP_END
  )
  raw.file <- file.path(OUTPUT_DIR, paste0(file.prefix, "_raw.csv"))
  summary.file <- file.path(
    OUTPUT_DIR,
    paste0(file.prefix, "_summary.csv")
  )
  plot.base <- file.path(
    OUTPUT_DIR,
    paste0(file.prefix, "_SE_SD_curve")
  )

  utils::write.csv(raw.result, raw.file, row.names = FALSE)
  utils::write.csv(summary.result, summary.file, row.names = FALSE)
  plot.calibration(summary.result, plot.base)

  runtime <- as.numeric(difftime(Sys.time(), timestart, units = "mins"))
  cat("================ Simulation Summary ================\n")
  cat("Run time (minutes):", round(runtime, 3), "\n")
  cat("Successful internal replications:",
      length(unique(raw.result$Rep)), "/", length(rep.ids), "\n")
  cat("Mean censoring rate:",
      round(mean(raw.result$CensorRate, na.rm = TRUE), 3), "\n")
  cat("Mean cure rate:",
      round(mean(raw.result$CureRate, na.rm = TRUE), 3), "\n")
  cat("Raw results:", raw.file, "\n")
  cat("Five-point summary:", summary.file, "\n")
  cat("Calibration plots:", paste0(plot.base, c(".pdf", ".png")), "\n")

  print.copy <- summary.result[
    , c(
      "ExternalN", "N_over_n", "Method", "Parameter",
      "Bias", "SE", "SD", "CP", "SE_SD", "ValidN"
    ),
    drop = FALSE
  ]
  numeric.columns <- vapply(print.copy, is.numeric, logical(1))
  print.copy[numeric.columns] <- lapply(
    print.copy[numeric.columns],
    round,
    digits = 3
  )
  print(print.copy, row.names = FALSE)
  cat("====================================================\n")
}, finally = {
  if (!is.null(cluster.object)) {
    closeCluster(cluster.object)
    Rmpi::mpi.quit()
  }
})
