
library(iterators)
library(Rmpi)
library(foreach)
library(doMPI)

library(survival)
library(MASS)

timestart <- Sys.time()

# =============================================================================
# Uncertainty-adjusted ICV simulation based on addICV011.R
#
# The external sample is generated once and then held fixed across all Monte
# Carlo replications.  The reported ICV1/ICV2/ICV4 estimators use the
# uncertainty-adjusted control-variate covariance
#
#   V22.U = V22 + (n / N_ex) * D.pi %*% Sigma.phi %*% D.pi,
#
# where Sigma.phi is estimated by applying the same subgroup KM procedure to
# the internal bootstrap samples, as in addICV_U01.R.  Both the ICV point
# estimator and its variance use V22.U.  Only WEE and these adjusted ICV
# estimators are reported.
# =============================================================================

# ---------------------------------------------------------------------------
# Setting selected through environment variables
# SETTING_ID: S1/S2/S3 x L/H x n300/n600, e.g. S1-L-n300
# ---------------------------------------------------------------------------
SETTING_ID <- Sys.getenv("SETTING_ID", "S1-L-n300")
REP_START <- as.integer(Sys.getenv("REP_START", "1"))
REP_END <- as.integer(Sys.getenv("REP_END", "500"))
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", file.path("addICV", SETTING_ID))
PARALLEL_SEED_BASE <- as.integer(Sys.getenv("PARALLEL_SEED_BASE", "20260825"))

valid.settings <- as.vector(outer(
  c("S1", "S2", "S3"),
  c("L-n300", "L-n600", "H-n300", "H-n600"),
  paste,
  sep = "-"
))

if (!SETTING_ID %in% valid.settings)
  stop("Invalid SETTING_ID. Use S1/S2/S3-L/H-n300/n600.")
if (is.na(REP_START) || is.na(REP_END) || REP_START < 1L || REP_END < REP_START)
  stop("REP_START and REP_END must satisfy 1 <= REP_START <= REP_END.")
if (is.na(PARALLEL_SEED_BASE))
  stop("PARALLEL_SEED_BASE must be an integer.")

setting.parts <- strsplit(SETTING_ID, "-", fixed = TRUE)[[1]]
pai.setting <- setting.parts[1]
censor.setting <- setting.parts[2]
n <- as.integer(sub("^n", "", setting.parts[3]))
rep.ids <- seq.int(REP_START, REP_END)
repetition <- length(rep.ids)

# Reproducible random-number settings
SEED_EXTERNAL_DATA <- 1L
SEED_INTERNAL_DATA <- 3L

# Simulation settings: two continuous covariates (original design)
beta.true <- c(1, 1)


if (censor.setting == "L") {
  lambdaC <- 0.50
  gamma.true <- c(-1.5, 1, 1)
} else {
  lambdaC <- 1.00
  gamma.true <- c(-2.5, 1, 1)
}

N_ex <- 20000L
B <- n * 4L
MIN_VALID_BOOT <- 30L
COV_EIGEN_FLOOR <- 1e-8

### WEE and ICV estimators using 1, 2, and 4 auxiliary time points ###
m.values <- c(1L, 2L, 4L)
icv.names <- paste0("ICV", m.values)
k.num <- 4L


# ---------------------------------------------------------------------------
# Link, subgroup, and baseline-hazard functions
# ---------------------------------------------------------------------------
pai <- switch(
  pai.setting,
  S1 = function(u) plogis(u),
  S2 = function(u) plogis(0.4 * (u^3 - u^2 - 0.9 * u + 1)),
  S3 = function(u) exp(-exp(-u))
)

subgf <- function(Z){
  rbind(
    Z[,1] <  1 & Z[,2] <  2,
    Z[,1] >= 1 & Z[,2] <  2,
    Z[,1] <  1 & Z[,2] >= 2,
    Z[,1] >= 1 & Z[,2] >= 2
  )
}

tau0 <- 2
kappa <- 1

H_0 <- function(t) {
  ans <- rep(Inf, length(t))
  id <- t < tau0

  ans[id] <- 2 * log1p(t[id]) -
    kappa * log1p(-t[id] / tau0)

  ans
}

# ---------------------------------------------------------------------------
# Data generation
# ---------------------------------------------------------------------------
Get.random <- function(n) {
  Z <- matrix(1, nrow = n, ncol = length(beta.true))
  Z[,1] <- runif(n, 0, 2)
  Z[,2] <- runif(n, 0, 4)
  
  T0 <- numeric(n); u <- runif(n, 0, 1)
  for (i in 1:n) {
    S0 <- function(t) {
      if (t >= tau0) return(0)
      exp(-H_0(t) - as.numeric(Z[i,] %*% beta.true) * t)
    }
    T0[i] <- uniroot(
      function(t) 1 - S0(t) - u[i],
      interval = c(0, tau0)
    )$root
  }
  
  Z.multi <- cbind(1, Z)
  peta <- pai(Z.multi %*% gamma.true)
  eta  <- rbinom(n, 1, peta)
  TT   <- ifelse(eta == 1, T0, Inf)
  
  C     <- rexp(n, lambdaC)
  Y     <- pmin(TT, C)
  delta <- as.numeric(eta == 1 & T0 <= C)
  
  data <- data.frame(Y = Y, delta = delta, Z = Z)
  data <- data[order(data[,1]), ]
  list(data = data, censor.rate = mean(delta == 0), cure.rate = mean(eta == 0))
}

Get.Data <- function(rep.ids){
  data <- vector('list', repetition)
  censor.rate <- cure.rate <- numeric(repetition)
  for (i in seq_along(rep.ids)) {
    set.seed(SEED_INTERNAL_DATA + rep.ids[i])
    tmp <- Get.random(n)
    data[[i]]      <- tmp$data
    censor.rate[i] <- tmp$censor.rate
    cure.rate[i]   <- tmp$cure.rate
  }
  list(data = data, censor.rate = censor.rate, cure.rate = cure.rate)
}

# ---------------------------------------------------------------------------
# Auxiliary time points and subgroup survival probabilities
# ---------------------------------------------------------------------------
subtf <- function(data, m.num) {
  data <- data[order(data[,1]), ]; Y <- data$Y
  tau <- which.max(data$Y[data$delta == 1])
  if (m.num == 1) Y[round(tau * 0.5)]
  else            Y[round((1:m.num) / (m.num + 1) * tau)]
}

getprob <- function(data, tstar, m.num){
  pro   <- matrix(0, k.num, m.num)
  group <- subgf(as.matrix(data[,-(1:2)]))
  for (i in 1:k.num) {
    gi  <- which(group[i, ]); gd <- data[gi, ]
    fit <- survival::survfit(survival::Surv(gd[,1], gd[,2]) ~ 1)   # KM
    pro[i, ] <- summary(fit, times = tstar, extend = TRUE)$surv
  }
  pro
}

protrue <- function(data, tstar, m.num){
  pro   <- matrix(0, k.num, m.num)
  Z     <- as.matrix(data[,-(1:2)])
  alpha <- sapply(tstar, H_0); group <- subgf(Z)
  for (m in 1:m.num) for (k in 1:k.num) {
    gi  <- which(group[k, ]); gd <- Z[gi, ]
    res <- numeric(nrow(gd))
    for (i in 1:nrow(gd)) {
      pu <- pai(gamma.true[1] + gd[i,1]*gamma.true[2] + gd[i,2]*gamma.true[3])
      res[i] <- 1 - pu + pu * exp(-alpha[m] - as.numeric(beta.true %*% gd[i,]) * tstar[m])
    }
    pro[k, m] <- mean(res)
  }
  pro
}

# ---------------------------------------------------------------------------
# Nonparametric weights
# ---------------------------------------------------------------------------
KK  <- function(x, hn) 3/4 * (1 - (x/hn)^2) * (abs(x/hn) <= 1)

NWweight2 <- function(Z, Z.sample, hn, lam = rep(1, nrow(Z.sample))){
  m  <- dim(Z)[1]; p1 <- dim(Z)[2]; nn <- dim(Z.sample)[1]
  Z <- as.data.frame(Z); Z.sample <- as.data.frame(Z.sample)
  MZ  <- array(apply(Z, 2, function(x) rep(x, nn)), dim = c(m, nn, p1))
  tMZ <- aperm(array(apply(Z.sample, 2, function(x) rep(x, m)), dim = c(nn, m, p1)), c(2,1,3))
  dMZ <- MZ - tMZ
  Kp <- array(NA, dim = c(m, nn, p1))
  Kp[,,1] <- KK(dMZ[,,1], hn[1])
  Kp[,,2] <- KK(dMZ[,,2], hn[2])
  prodKp <- apply(Kp, c(1,2), prod)
  Lam    <- matrix(lam, nrow = m, ncol = nn, byrow = TRUE)
  lamKp  <- prodKp * Lam
  rs <- matrix(apply(lamKp, 1, sum), nrow = m, ncol = nn, byrow = FALSE)
  rs[rs == 0] <- 1e-10
  lamKp / rs
}

S.Ti <- function(data, w){
  m <- nrow(w); nn <- ncol(w)
  k <- cbind(rep(0, m), w[, 1:(nn-1)])
  cumsum.k <- t(apply(k, 1, cumsum))
  cum.w <- 1 - cumsum.k
  cum.w <- cum.w * (cum.w > w + 1e-5) + (w + 1e-4) * (cum.w <= w + 1e-5)
  q <- 1 - w / cum.w
  delta <- matrix(data$delta, nrow = m, ncol = nn, byrow = TRUE)
  q <- replace(q, delta == 0, 1)
  q <- t(apply(q, 1, cumprod))
  Y <- data$Y; r <- rle(Y); ei <- cumsum(r$lengths); li <- rep(ei, r$lengths)
  q[cbind(1:m, li[1:m])]
}

P.T <- function(data, w){
  m <- nrow(w); nn <- ncol(w)
  k <- cbind(rep(0, m), w[, 1:(nn-1)])
  cumsum.k <- t(apply(k, 1, cumsum))
  cum.w <- 1 - cumsum.k
  w     <- replace(w, cum.w == 0, 0)
  cum.w <- replace(cum.w, cum.w == 0, 1)
  q <- 1 - w / cum.w
  delta <- matrix(data$delta, nrow = m, ncol = nn, byrow = TRUE)
  q <- replace(q, delta == 0, 1)
  q <- t(apply(q, 1, prod))
  as.numeric(replace(q, q <= 0, 0))
}

# ---------------------------------------------------------------------------
# Weighted estimating equation and control variables
# ---------------------------------------------------------------------------
get.est <- function(data, aux.info){
  data  <- data[order(data[,1]), ]
  Y     <- data$Y; delta <- data$delta
  Z     <- as.matrix(data[,-(1:2)]); p <- ncol(Z); nn <- nrow(data)
  T1.max <- max(Y[delta == 1])
  # Component-specific sample-SD scaling for the two continuous covariates.
  h <- 3 * nn^(-2/7) * apply(Z, 2L, stats::sd)
  if (any(!is.finite(h)) || any(h <= 0)) {
    stop("All bandwidths must be finite and strictly positive.")
  }
  
  wght      <- NWweight2(Z = Z, Z.sample = Z, hn = h)
  S.beran   <- S.Ti(data, wght)
  P.T_inf   <- P.T(data, wght)
  S.safe    <- S.beran; S.safe[S.safe == 0] <- 1e-8
  case.wght <- delta + (1 - delta) * (S.beran - P.T_inf) / S.safe * (Y <= T1.max)
  
  dif <- diff(c(0, Y))
  rle_Y <- rle(Y)
  start_idx     <- c(1, 1 + head(cumsum(rle_Y$lengths), -1))
  first_indices <- rep(start_idx, rle_Y$lengths)
  
  Zbarr <- matrix(0, nn, p)
  for (j in 1:nn) {
    idx <- first_indices[j]:nn
    wj  <- case.wght[idx]; sw <- sum(wj)
    if (sw > 0) Zbarr[j, ] <- colSums(Z[idx, , drop = FALSE] * wj) / sw
  }
  
  fenzi <- colSums(case.wght * delta * (Z - Zbarr))
  fenmu <- matrix(0, p, p)
  for (j in 1:nn) if (dif[j] > 1e-8) {
    idx <- first_indices[j]:nn
    wj  <- case.wght[idx]
    Zc  <- sweep(Z[idx, , drop = FALSE], 2, Zbarr[j, ], "-")
    fenmu <- fenmu + (t(Zc) %*% sweep(Zc, 1, wj, "*")) * dif[j]
  }
  
  beta.WEE <- as.vector(fenzi %*% solve(fenmu))
  if (any(!is.finite(beta.WEE))) return(list(ok = FALSE, reason = "non-finite beta"))
  
  lambest <- function(ts) {
    
    eps <- 1e-12
    
    risk.w <- rev(cumsum(rev(case.wght)))
    
    zbeta <- as.vector(Z %*% beta.WEE)
    risk.wz <- rev(cumsum(rev(case.wght * zbeta)))

    event.num <- delta * case.wght * (Y <= ts)
    risk.event <- risk.w[first_indices]
    
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
      
      if (risk.w[1] <= eps) {
        stop("Zero initial risk-set denominator.")
      }
      
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
    
    if (!is.finite(ans)) {
      stop("Non-finite baseline cumulative hazard.")
    }
    
    ans
  }
  
  subdata <- cbind(1 - P.T_inf, Z)
  group   <- subgf(Z)
  pai.hat <- apply(group, 1, mean)

  phi.by.method <- vector("list", length(aux.info))
  names(phi.by.method) <- names(aux.info)

  for (method in names(aux.info)) {
    tstar <- aux.info[[method]]$tstar
    tpro <- aux.info[[method]]$tpro
    m.num <- length(tstar)

    alpha <- tryCatch(
      vapply(tstar, lambest, numeric(1)),
      error = function(e) NULL
    )
    if (is.null(alpha) || any(!is.finite(alpha))) {
      phi.by.method[[method]] <- rep(NA_real_, k.num * m.num)
      next
    }

    phi <- matrix(0, k.num, m.num)
    for (m in seq_len(m.num)) for (k in seq_len(k.num)) {
      gi <- which(group[k, ])
      if (length(gi) > 0) {
        gd <- subdata[gi, , drop = FALSE]
        res <- 0
        for (i in seq_len(nrow(gd))) {
          term <- (1 - gd[i,1]) +
            gd[i,1] * exp(
              -alpha[m] -
                as.numeric(beta.WEE %*% gd[i,-1]) * tstar[m]
            ) - tpro[k, m]
          res <- res + term
        }
        phi[k, m] <- res / nn
      }
    }

    est.phi <- as.vector(phi)
    if (any(!is.finite(est.phi))) {
      phi.by.method[[method]] <- rep(NA_real_, k.num * m.num)
    } else {
      phi.by.method[[method]] <- est.phi
    }
  }

  list(
    ok = TRUE,
    est.beta = beta.WEE,
    est.phi = phi.by.method,
    pai = pai.hat
  )
}

# ---------------------------------------------------------------------------
# Bootstrap covariance estimation
# ---------------------------------------------------------------------------
boot.cov <- function(data, aux.info, B){
  nn <- nrow(data)
  p <- length(beta.true)
  beta.boot <- matrix(NA, B, p)

  # Bootstrap versions of the internal control variates.
  phi.boot <- lapply(
    aux.info,
    function(x) matrix(NA_real_, B, k.num * length(x$tstar))
  )

  # The same subgroup KM summaries, evaluated on the internal bootstrap
  # samples, estimate Sigma.phi on the sqrt(n) scale under homogeneity.
  km.boot <- lapply(
    aux.info,
    function(x) matrix(NA_real_, B, k.num * length(x$tstar))
  )

  for (b in seq_len(B)) {
    bs  <- data[sample(1:nn, nn, replace = TRUE), ]
    fit <- tryCatch(get.est(bs, aux.info), error = function(e) list(ok = FALSE))
    if (!isTRUE(fit$ok)) next
    beta.boot[b, ] <- fit$est.beta
    for (method in names(aux.info)) {
      phi.boot[[method]][b, ] <- fit$est.phi[[method]]
      km.boot[[method]][b, ] <- tryCatch(
        as.vector(
          getprob(
            bs,
            aux.info[[method]]$tstar,
            length(aux.info[[method]]$tstar)
          )
        ),
        error = function(e) rep(
          NA_real_,
          k.num * length(aux.info[[method]]$tstar)
        )
      )
    }
  }

  ok.beta <- complete.cases(beta.boot)
  ok.phi <- lapply(
    phi.boot,
    function(x) ok.beta & complete.cases(x)
  )

  # A common complete-bootstrap set is required because V11, V12, V22, and
  # Sigma.phi must describe the same resampling distribution.
  ok.km <- lapply(
    km.boot,
    function(x) complete.cases(x)
  )
  ok.icv <- Map(
    function(ok1, ok2) ok1 & ok2,
    ok.phi,
    ok.km
  )
  nvalid <- c(
    WEE = sum(ok.beta),
    vapply(ok.icv, sum, numeric(1))
  )

  list(
    beta = beta.boot,
    phi = phi.boot,
    km = km.boot,
    ok.beta = ok.beta,
    ok.phi = ok.phi,
    ok.km = ok.km,
    ok.icv = ok.icv,
    nvalid = nvalid,
    ndrop = B - nvalid
  )
}

# ---------------------------------------------------------------------------
# Matrix utilities
# ---------------------------------------------------------------------------
symmetrize <- function(A) (A + t(A)) / 2

# Stable inverse for an estimated covariance matrix.  The eigenvalue floor is
# a numerical safeguard and is applied before both the point-estimation and
# variance calculations.
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

# ---------------------------------------------------------------------------
# One Monte Carlo replication
# ---------------------------------------------------------------------------
Repeat <- function(r, rep.id, aux.info){
  data      <- Data$data[[r]]
  cens.rate <- Data$censor.rate[r]
  cure.rate <- Data$cure.rate[r]
  nn        <- nrow(data)
  fit0 <- get.est(data, aux.info)
  if (!isTRUE(fit0$ok)) return(NULL)

  beta.WEE <- fit0$est.beta

  bc <- boot.cov(data, aux.info, B)
  if (bc$nvalid[["WEE"]] < MIN_VALID_BOOT) return(NULL)

  beta.boot.WEE <- bc$beta[bc$ok.beta, , drop = FALSE]
  V11.WEE <- cov(beta.boot.WEE * sqrt(nn))

  if (any(!is.finite(V11.WEE))) return(NULL)
  se.WEE <- sqrt(pmax(diag(V11.WEE / nn), 0))

  out <- list(
    rep = rep.id,
    beta1.WEE = beta.WEE[1], beta2.WEE = beta.WEE[2],
    se1.WEE = se.WEE[1], se2.WEE = se.WEE[2],
    external.n = N_ex,
    rho.n.over.N = nn / N_ex
  )

  for (method in names(aux.info)) {
    phi.est <- fit0$est.phi[[method]]
    beta.ICV <- rep(NA_real_, length(beta.true))
    se.ICV <- rep(NA_real_, length(beta.true))
    method.ok <- FALSE
    min.eigen.VU <- NA_real_
    min.eigen.Sigma.phi <- NA_real_
    regularized.VU <- NA

    if (all(is.finite(phi.est)) &&
        bc$nvalid[[method]] >= MIN_VALID_BOOT) {
      ok <- bc$ok.icv[[method]]
      beta.boot <- bc$beta[ok, , drop = FALSE]
      phi.boot <- bc$phi[[method]][ok, , drop = FALSE]
      km.boot <- bc$km[[method]][ok, , drop = FALSE]

      V11 <- stats::cov(beta.boot * sqrt(nn))
      V22 <- symmetrize(stats::cov(phi.boot * sqrt(nn)))
      V12 <- stats::cov(
        phi.boot * sqrt(nn),
        beta.boot * sqrt(nn)
      )
      Sigma.phi <- symmetrize(stats::cov(km.boot * sqrt(nn)))

      m.num <- length(aux.info[[method]]$tstar)
      q <- length(phi.est)
      pai.vector <- rep(fit0$pai, times = m.num)
      if (length(pai.vector) != q) {
        stop("The dimension of D.pi does not match the auxiliary vector.")
      }
      D.pi <- diag(pai.vector, nrow = q)

      # Uncertainty-adjusted covariance used in both the ICV coefficient and
      # its estimated asymptotic variance.
      V22.U <- symmetrize(
        V22 + (nn / N_ex) * D.pi %*% Sigma.phi %*% D.pi
      )

      finite.matrices <-
        all(is.finite(V11)) &&
        all(is.finite(V12)) &&
        all(is.finite(V22)) &&
        all(is.finite(Sigma.phi)) &&
        all(is.finite(V22.U))

      if (finite.matrices) {
        inv.VU <- regularized_inverse(V22.U)
        min.eigen.VU <- inv.VU$min.eigen
        min.eigen.Sigma.phi <- min(
          eigen(Sigma.phi, symmetric = TRUE, only.values = TRUE)$values
        )
        regularized.VU <- inv.VU$regularized

        beta.candidate <- as.vector(
          beta.WEE - t(V12) %*% inv.VU$inverse %*% phi.est
        )
        var.candidate <- symmetrize(
          (V11 - t(V12) %*% inv.VU$inverse %*% V12) / nn
        )

        method.ok <-
          all(is.finite(beta.candidate)) &&
          all(is.finite(var.candidate)) &&
          all(diag(var.candidate) >= -1e-10)

        if (method.ok) {
          beta.ICV <- beta.candidate
          se.ICV <- sqrt(pmax(diag(var.candidate), 0))
        }
      }
    }

    out[[paste0("beta1.", method)]] <- beta.ICV[1]
    out[[paste0("beta2.", method)]] <- beta.ICV[2]
    out[[paste0("se1.", method)]] <- se.ICV[1]
    out[[paste0("se2.", method)]] <- se.ICV[2]
    out[[paste0("valid.", method)]] <- method.ok
    out[[paste0("min.eigen.VU.", method)]] <- min.eigen.VU
    out[[paste0("min.eigen.SigmaPhi.", method)]] <- min.eigen.Sigma.phi
    out[[paste0("regularized.VU.", method)]] <- regularized.VU
  }

  for (method in names(bc$nvalid)) {
    out[[paste0("nvalid.boot.", method)]] <- bc$nvalid[[method]]
    out[[paste0("ndrop.boot.", method)]] <- bc$ndrop[[method]]
  }

  c(
    out,
    list(
      cens.rate = cens.rate,
      cure.rate = cure.rate
    )
  )
}

# ---------------------------------------------------------------------------
# Simulation summary
# ---------------------------------------------------------------------------
GetSummary <- function(result_list, beta.true){
  options(digits = 3)

  valid <- Filter(
    function(x) !is.null(x) &&
      !inherits(x, "error") &&
      !is.null(x$beta1.WEE),
    result_list
  )

  if (length(valid) == 0L) stop("All iterations failed.")

  cens.rate <- vapply(valid, function(x) x$cens.rate, numeric(1))
  cure.rate <- vapply(valid, function(x) x$cure.rate, numeric(1))
  bt <- as.numeric(beta.true)
  methods <- c("WEE", icv.names)
  nvalid.boot <- vapply(methods, function(method) {
    key <- paste0("nvalid.boot.", method)
    vapply(valid, function(x) x[[key]], numeric(1))
  }, numeric(length(valid)))
  ndrop.boot <- vapply(methods, function(method) {
    key <- paste0("ndrop.boot.", method)
    vapply(valid, function(x) x[[key]], numeric(1))
  }, numeric(length(valid)))

  extract_pair <- function(prefix, method) {
    k1 <- paste0(prefix, "1.", method)
    k2 <- paste0(prefix, "2.", method)
    t(vapply(valid, function(x) c(x[[k1]], x[[k2]]), numeric(2)))
  }

  calc <- function(method){
    est <- extract_pair("beta", method)
    se <- extract_pair("se", method)
    ok <- complete.cases(est) & complete.cases(se)

    if (sum(ok) < 2L) {
      return(list(
        Bias = rep(NA_real_, 2),
        SD = rep(NA_real_, 2),
        SE = rep(NA_real_, 2),
        CP = rep(NA_real_, 2),
        MSE = rep(NA_real_, 2),
        ValidN = rep(sum(ok), 2)
      ))
    }

    est <- est[ok, , drop = FALSE]
    se <- se[ok, , drop = FALSE]
    bias <- colMeans(est) - bt
    SD   <- apply(est, 2, sd)
    SE   <- colMeans(se)
    truth.matrix <- matrix(bt, nrow(est), length(bt), byrow = TRUE)
    z <- qnorm(0.975)
    CP <- colMeans(
      (est - z * se <= truth.matrix) &
        (truth.matrix <= est + z * se)
    )

    list(
      Bias = bias,
      SD = SD,
      SE = SE,
      CP = CP,
      MSE = bias^2 + SD^2,
      ValidN = rep(sum(ok), 2)
    )
  }

  summary.by.method <- lapply(methods, calc)
  names(summary.by.method) <- methods

  result <- do.call(
    rbind,
    lapply(methods, function(method) {
      s <- summary.by.method[[method]]
      data.frame(
        Method = method,
        Parameter = c("beta1", "beta2"),
        Bias = s$Bias,
        SD = s$SD,
        SE = s$SE,
        CP = s$CP,
        MSE = s$MSE,
        ValidN = as.integer(s$ValidN),
        CureRate = mean(cure.rate),
        CensorRate = mean(cens.rate),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(result) <- NULL

  runtime <- difftime(Sys.time(), timestart, units = "mins")

  cat("======================= Simulation Summary =======================\n")
  cat("Run Time          :", round(runtime, 3), "mins\n")
  cat("Sample size       :", n, "\n")
  cat("Successful WEE    :", length(valid), "/", length(result_list), "\n")
  cat(
    "Valid ICV runs   :",
    paste(
      paste0(
        icv.names,
        "=",
        vapply(
          icv.names,
          function(method) summary.by.method[[method]]$ValidN[1],
          numeric(1)
        )
      ),
      collapse = ", "
    ),
    "\n"
  )
  cat("Censor rate       :", round(mean(cens.rate), 3), "\n")
  cat("Cure rate         :", round(mean(cure.rate), 3), "\n")
  cat(
    "Valid boot (mean) :",
    paste(
      paste0(methods, "=", round(colMeans(nvalid.boot), 1)),
      collapse = ", "
    ),
    "\n"
  )
  cat(
    "Dropped boot(mean):",
    paste(
      paste0(methods, "=", round(colMeans(ndrop.boot), 1)),
      collapse = ", "
    ),
    "\n"
  )
  cat("\n")
  result.print <- result
  numeric.columns <- vapply(result.print, is.numeric, logical(1))
  result.print[numeric.columns] <- lapply(
    result.print[numeric.columns],
    round,
    digits = 3
  )
  print(result.print, row.names = FALSE)
  cat("==================================================================\n")

  result
}

# ---------------------------------------------------------------------------
# External auxiliary information
# ---------------------------------------------------------------------------
set.seed(SEED_EXTERNAL_DATA)
sdata.external <- Get.random(N_ex)$data
aux.info <- setNames(vector("list", length(m.values)), icv.names)

for (j in seq_along(m.values)) {
  m.num <- m.values[j]
  method <- icv.names[j]
  # Keep the addICV011 construction unchanged: tstar and tpro are computed
  # once from the single external sample, with tpro obtained by protrue().
  tstar <- subtf(sdata.external, m.num)
  tpro <- protrue(sdata.external, tstar, m.num)

  aux.info[[method]] <- list(
    tstar = tstar,
    tpro = tpro
  )
}

# ---------------------------------------------------------------------------
# Internal Monte Carlo data
# ---------------------------------------------------------------------------
Data <- Get.Data(rep.ids)

# ---------------------------------------------------------------------------
# MPI parallel computation
# ---------------------------------------------------------------------------
packagelist <- c("survival", "MASS")

exports <- c(
  "Repeat", "get.est", "boot.cov", "symmetrize",
  "regularized_inverse", "getprob",
  "NWweight2", "S.Ti", "P.T", "KK", "subgf", "pai",
  "Data", "aux.info", "beta.true", "B", "N_ex",
  "k.num", "MIN_VALID_BOOT", "COV_EIGEN_FLOOR",
  "PARALLEL_SEED_BASE", "rep.ids"
)

cl <- startMPIcluster()
registerDoMPI(cl)

cat("MPI workers:", foreach::getDoParWorkers(), "\n")
cat("SETTING_ID:", SETTING_ID, "\n")
cat("Replications:", REP_START, "to", REP_END, "\n")
cat("External sample size:", N_ex, "(generated once)\n")
cat("Bandwidth: 3 * n^(-2/7) * sample SD(Z_j)\n")
cat("Output directory:", OUTPUT_DIR, "\n")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

tryCatch({

  MyResult <- foreach(
    r = seq_along(rep.ids),
    .packages = packagelist,
    .export = exports,
    .errorhandling = "pass",
    .inorder = TRUE
  ) %dopar% {

    rep.id <- rep.ids[r]
    set.seed(PARALLEL_SEED_BASE + rep.id)

    Repeat(
      r = r,
      rep.id = rep.id,
      aux.info = aux.info
    )
  }

  valid.output <- Filter(
    function(x) !is.null(x) && !inherits(x, "error") && !is.null(x$beta1.WEE),
    MyResult
  )

  if (length(valid.output) == 0L) stop("All iterations failed.")

  raw.file <- file.path(
    OUTPUT_DIR,
    sprintf("%s_rep%d-%d.csv", SETTING_ID, REP_START, REP_END)
  )
  summary.file <- file.path(
    OUTPUT_DIR,
    sprintf("%s_rep%d-%d_summary.csv", SETTING_ID, REP_START, REP_END)
  )

  write.csv(
    do.call(rbind, lapply(valid.output, as.data.frame)),
    file = raw.file,
    row.names = FALSE
  )
  
  MySummary <- GetSummary(MyResult, beta.true)
  write.csv(MySummary, file = summary.file, row.names = FALSE)
  print(MySummary)

}, finally = {

  closeCluster(cl)
  Rmpi::mpi.quit()
})
