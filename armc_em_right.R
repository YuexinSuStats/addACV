# EM estimation for a right-censored additive-risk mixture-cure model
#
# Susceptible hazard:
#   lambda_u(t | Z) = lambda_0(t) + Z^T beta, beta >= 0.
# Incidence model used for fitting:
#   P(eta = 1 | X) = plogis(X^T gamma).
#
# The baseline cumulative hazard is represented by a monotone Bernstein sieve.
# If B_{l,N} is the Bernstein basis and
#   theta_m = gamma_m - gamma_{m-1} >= 0,
# then
#   Lambda_0(t) = sum_m theta_m A_m(t),
# where A_m(t) is a beta CDF. Its derivative a_m(t) is a beta density.
# A second latent variable allocates every observed failure among the baseline
# components and the covariate-specific additive hazards. This gives closed-form
# M-step updates for theta and beta.
#
# Important: the Poisson-superposition EM below requires Z >= 0 and beta >= 0.
# It must not be applied unchanged to centered or negative latency covariates.
#
# Run one simulation setting from the shell, for example:
#   CENSOR_SETTING=L PAI_SETTING=S1 N=300 REPETITION=500 \
#     Rscript armc_em_right.R --run

log1pexp <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

logspace_add <- function(logx, logy) {
  m <- pmax(logx, logy)
  m + log(exp(logx - m) + exp(logy - m))
}

H_0 <- function(t, tau0 = 2, kappa = 1) {
  t <- as.numeric(t)
  ans <- rep(Inf, length(t))
  id <- is.finite(t) & t >= 0 & t < tau0
  ans[id] <- 2 * log1p(t[id]) - kappa * log1p(-t[id] / tau0)
  ans
}

make_pai <- function(pai.setting = c("S1", "S2", "S3")) {
  pai.setting <- match.arg(pai.setting)
  switch(
    pai.setting,
    S1 = function(u) plogis(u),
    S2 = function(u) plogis(0.4 * (u^3 - u^2 - 0.9 * u + 1)),
    S3 = function(u) exp(-exp(-u))
  )
}

simulation_parameters <- function(censor.setting = c("L", "H")) {
  censor.setting <- match.arg(censor.setting)
  if (censor.setting == "L") {
    list(lambdaC = 0.50, gamma.true = c(-1.5, 1, 1))
  } else {
    list(lambdaC = 1.00, gamma.true = c(-2.5, 1, 1))
  }
}

Get.random <- function(n,
                       censor.setting = c("L", "H"),
                       pai.setting = c("S1", "S2", "S3"),
                       seed = NULL,
                       beta.true = c(1, 1),
                       tau0 = 2,
                       kappa = 1) {
  censor.setting <- match.arg(censor.setting)
  pai.setting <- match.arg(pai.setting)
  if (!is.null(seed)) set.seed(seed)
  if (length(beta.true) != 2L || any(beta.true < 0)) {
    stop("This generator expects two nonnegative latency coefficients.")
  }

  pars <- simulation_parameters(censor.setting)
  pai <- make_pai(pai.setting)

  Z <- cbind(
    Z1 = runif(n, 0, 2),
    Z2 = runif(n, 0, 4)
  )
  X.inc <- cbind(`(Intercept)` = 1, Z)
  peta <- as.numeric(pai(drop(X.inc %*% pars$gamma.true)))
  eta <- rbinom(n, 1, peta)

  u <- runif(n)
  T0 <- numeric(n)
  for (i in seq_len(n)) {
    zi.beta <- sum(Z[i, ] * beta.true)
    root.fun <- function(t) {
      S0 <- exp(-H_0(t, tau0 = tau0, kappa = kappa) - zi.beta * t)
      1 - S0 - u[i]
    }
    T0[i] <- uniroot(
      root.fun,
      interval = c(0, tau0),
      tol = .Machine$double.eps^0.5
    )$root
  }

  TT <- ifelse(eta == 1, T0, Inf)
  C <- rexp(n, pars$lambdaC)
  Y <- pmin(TT, C)
  delta <- as.integer(eta == 1 & T0 <= C)

  dat <- data.frame(Y = Y, delta = delta, Z1 = Z[, 1], Z2 = Z[, 2])
  ord <- order(dat$Y)
  dat <- dat[ord, , drop = FALSE]
  rownames(dat) <- NULL

  list(
    data = dat,
    censor.rate = mean(delta == 0),
    cure.rate = mean(eta == 0),
    eta = eta[ord],
    peta = peta[ord],
    beta.true = beta.true,
    gamma.true = pars$gamma.true,
    lambdaC = pars$lambdaC,
    tau0 = tau0,
    pai.setting = pai.setting,
    censor.setting = censor.setting
  )
}

Get.Data <- function(rep.ids,
                     n,
                     censor.setting = c("L", "H"),
                     pai.setting = c("S1", "S2", "S3"),
                     SEED_INTERNAL_DATA = 20260921L) {
  censor.setting <- match.arg(censor.setting)
  pai.setting <- match.arg(pai.setting)
  out <- lapply(rep.ids, function(rep.id) {
    Get.random(
      n = n,
      censor.setting = censor.setting,
      pai.setting = pai.setting,
      seed = SEED_INTERNAL_DATA + as.integer(rep.id)
    )
  })
  list(
    data = lapply(out, `[[`, "data"),
    censor.rate = vapply(out, `[[`, numeric(1), "censor.rate"),
    cure.rate = vapply(out, `[[`, numeric(1), "cure.rate")
  )
}

bernstein_hazard_basis <- function(t, degree, tau) {
  if (degree < 1L || degree != as.integer(degree)) {
    stop("degree must be a positive integer.")
  }
  if (!is.finite(tau) || tau < max(t)) {
    stop("tau must be finite and no smaller than the largest supplied time.")
  }
  s <- pmin(pmax(as.numeric(t) / tau, 0), 1)
  n <- length(s)
  A <- matrix(0, nrow = n, ncol = degree)
  a <- matrix(0, nrow = n, ncol = degree)
  for (m in seq_len(degree)) {
    A[, m] <- pbeta(s, shape1 = m, shape2 = degree - m + 1)
    a[, m] <- dbeta(s, shape1 = m, shape2 = degree - m + 1) / tau
  }
  colnames(A) <- colnames(a) <- paste0("B", seq_len(degree))
  list(cumulative = A, hazard = a, scaled.time = s)
}

fractional_logit_fit <- function(X,
                                 response,
                                 start,
                                 tol = 1e-10,
                                 maxit = 50L,
                                 ridge = 1e-9) {
  X <- as.matrix(X)
  response <- as.numeric(response)
  gamma <- as.numeric(start)
  qfun <- function(g) {
    eta <- drop(X %*% g)
    sum(response * eta - log1pexp(eta))
  }

  for (iter in seq_len(maxit)) {
    eta <- drop(X %*% gamma)
    p <- plogis(eta)
    score <- drop(crossprod(X, response - p))
    information <- crossprod(X, X * (p * (1 - p)))
    information <- information + diag(ridge, ncol(X))
    step <- tryCatch(
      solve(information, score),
      error = function(e) qr.solve(information, score, tol = 1e-12)
    )

    old.q <- qfun(gamma)
    step.factor <- 1
    repeat {
      candidate <- gamma + step.factor * step
      if (qfun(candidate) >= old.q - 1e-12 || step.factor < 2^-20) break
      step.factor <- step.factor / 2
    }
    gamma <- candidate
    if (max(abs(step.factor * step)) < tol) break
  }
  gamma
}

observed_loglik_armc <- function(gamma,
                                 beta,
                                 theta,
                                 Y,
                                 delta,
                                 Z,
                                 X.inc,
                                 A,
                                 a,
                                 beyond.support = rep(FALSE, length(Y))) {
  linear.inc <- drop(X.inc %*% gamma)
  linear.lat <- drop(Z %*% beta)
  cumulative.hazard <- drop(A %*% theta) + Y * linear.lat
  baseline.hazard <- drop(a %*% theta)
  total.hazard <- baseline.hazard + linear.lat

  if (any(delta == 1 & beyond.support)) return(-Inf)
  if (any(total.hazard[delta == 1] <= 0) ||
      any(!is.finite(total.hazard[delta == 1]))) {
    return(-Inf)
  }

  log.pi <- -log1pexp(-linear.inc)
  log.one.minus.pi <- -log1pexp(linear.inc)
  ll.failure <- log.pi + log(total.hazard) - cumulative.hazard
  ll.censored <- logspace_add(log.one.minus.pi, log.pi - cumulative.hazard)
  ll.censored[beyond.support] <- log.one.minus.pi[beyond.support]
  sum(ifelse(delta == 1, ll.failure, ll.censored))
}

em_armc_right <- function(Y,
                          delta,
                          Z,
                          X.inc = NULL,
                          degree = NULL,
                          tau = NULL,
                          support.end = NULL,
                          start = NULL,
                          control = list()) {
  Y <- as.numeric(Y)
  delta <- as.integer(delta)
  Z <- as.matrix(Z)
  storage.mode(Z) <- "double"
  n <- length(Y)
  if (nrow(Z) != n || length(delta) != n) stop("Y, delta, and Z have incompatible dimensions.")
  if (any(!is.finite(Y)) || any(Y <= 0)) stop("All observed times must be positive and finite.")
  if (any(!delta %in% c(0L, 1L))) stop("delta must contain only 0 and 1.")
  if (sum(delta) == 0L) stop("At least one observed failure is required.")
  if (any(!is.finite(Z))) stop("Z contains a non-finite value.")
  if (any(Z < 0)) {
    stop("This closed-form EM requires nonnegative latency covariates Z.")
  }
  if (is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(ncol(Z)))

  if (is.null(X.inc)) X.inc <- cbind(`(Intercept)` = 1, Z)
  X.inc <- as.matrix(X.inc)
  storage.mode(X.inc) <- "double"
  if (nrow(X.inc) != n || any(!is.finite(X.inc))) {
    stop("X.inc must be a finite matrix with one row per subject.")
  }
  if (is.null(colnames(X.inc))) colnames(X.inc) <- paste0("X", seq_len(ncol(X.inc)))

  defaults <- list(
    maxit = 2000L,
    tol = 1e-6,
    min.iter = 5L,
    min.rate = 1e-10,
    logistic.tol = 1e-10,
    logistic.maxit = 50L,
    trace = FALSE,
    trace.every = 50L,
    monotone.tol = 1e-8
  )
  control <- modifyList(defaults, control)

  if (is.null(degree)) degree <- max(3L, floor(n^(1 / 4)))
  degree <- as.integer(degree)
  event <- delta == 1L
  if (!is.null(support.end)) {
    support.end <- as.numeric(support.end)
    if (length(support.end) != 1L || !is.finite(support.end) || support.end <= 0) {
      stop("support.end must be one positive finite number.")
    }
    if (any(Y[event] >= support.end)) {
      stop("Every observed failure must occur strictly before support.end.")
    }
    if (!is.null(tau) && abs(tau - support.end) > 1e-10 * (1 + support.end)) {
      stop("When support.end is supplied, tau must be NULL or equal to support.end.")
    }
    tau <- support.end
    basis.time <- pmin(Y, support.end)
    beyond.support <- Y >= support.end
  } else {
    if (is.null(tau)) tau <- max(Y) * (1 + 1e-8) + 1e-12
    if (tau <= max(Y)) stop("Without a finite support endpoint, tau must exceed max(Y).")
    basis.time <- Y
    beyond.support <- rep(FALSE, n)
  }
  basis <- bernstein_hazard_basis(basis.time, degree = degree, tau = tau)
  A <- basis$cumulative
  a <- basis$hazard

  p <- ncol(Z)
  q <- ncol(X.inc)
  if (is.null(start)) start <- list()
  if (is.null(start$gamma)) {
    event.rate <- min(max(mean(delta), 0.02), 0.98)
    gamma0 <- c(qlogis(event.rate), rep(0, q - 1L))
    gamma <- fractional_logit_fit(
      X.inc,
      delta,
      start = gamma0,
      tol = control$logistic.tol,
      maxit = control$logistic.maxit
    )
  } else {
    gamma <- as.numeric(start$gamma)
  }
  beta <- if (is.null(start$beta)) rep(0.5, p) else as.numeric(start$beta)
  theta <- if (is.null(start$theta)) rep(1 / degree, degree) else as.numeric(start$theta)

  if (length(gamma) != q || length(beta) != p || length(theta) != degree) {
    stop("The supplied start values have incorrect dimensions.")
  }
  if (any(beta <= 0) || any(theta <= 0)) {
    stop("All starting values for beta and theta must be strictly positive.")
  }

  ll.trace <- rep(NA_real_, control$maxit + 1L)
  ll.trace[1] <- observed_loglik_armc(
    gamma, beta, theta, Y, delta, Z, X.inc, A, a, beyond.support
  )
  converged <- FALSE
  monotonicity.violations <- 0L
  max.change <- Inf
  finite.change <- Inf
  baseline.change <- Inf

  for (iter in seq_len(control$maxit)) {
    linear.inc <- drop(X.inc %*% gamma)
    cumulative.hazard <- drop(A %*% theta) + Y * drop(Z %*% beta)

    # Posterior P(eta_i = 1 | observed data). It equals one for failures.
    susceptible.weight <- ifelse(
      event,
      1,
      ifelse(beyond.support, 0, plogis(linear.inc - cumulative.hazard))
    )

    # Conditional allocation of each failure to a nonnegative hazard component.
    base.parts <- sweep(a[event, , drop = FALSE], 2, theta, `*`)
    cov.parts <- sweep(Z[event, , drop = FALSE], 2, beta, `*`)
    event.hazard <- rowSums(base.parts) + rowSums(cov.parts)
    if (any(!is.finite(event.hazard)) || any(event.hazard <= 0)) {
      stop("A nonpositive fitted failure hazard occurred during the E-step.")
    }
    base.allocation <- base.parts / event.hazard
    cov.allocation <- cov.parts / event.hazard

    gamma.new <- fractional_logit_fit(
      X.inc,
      susceptible.weight,
      start = gamma,
      tol = control$logistic.tol,
      maxit = control$logistic.maxit
    )

    theta.denominator <- colSums(A * susceptible.weight)
    beta.denominator <- colSums(Z * (susceptible.weight * Y))
    if (any(theta.denominator <= 0) || any(beta.denominator <= 0)) {
      stop("A zero exposure denominator occurred during the M-step.")
    }
    theta.new <- colSums(base.allocation) / theta.denominator
    beta.new <- colSums(cov.allocation) / beta.denominator
    theta.new <- pmax(theta.new, control$min.rate)
    beta.new <- pmax(beta.new, control$min.rate)

    finite.change <- max(c(
      abs(gamma.new - gamma) / (1 + abs(gamma)),
      abs(beta.new - beta) / (1 + abs(beta))
    ))
    baseline.change <- max(abs(theta.new - theta) / (1 + abs(theta)))
    max.change <- max(finite.change, baseline.change)

    gamma <- gamma.new
    beta <- beta.new
    theta <- theta.new
    ll.trace[iter + 1L] <- observed_loglik_armc(
      gamma, beta, theta, Y, delta, Z, X.inc, A, a, beyond.support
    )
    relative.ll.change <- abs(ll.trace[iter + 1L] - ll.trace[iter]) /
      (1 + abs(ll.trace[iter]))

    if (ll.trace[iter + 1L] < ll.trace[iter] -
        control$monotone.tol * (1 + abs(ll.trace[iter]))) {
      monotonicity.violations <- monotonicity.violations + 1L
    }
    if (isTRUE(control$trace) &&
        (iter == 1L || iter %% control$trace.every == 0L)) {
      message(sprintf(
        "EM iteration %d: logLik = %.8f, max relative change = %.3e",
        iter, ll.trace[iter + 1L], max.change
      ))
    }
    if (iter >= control$min.iter &&
        finite.change < control$tol && relative.ll.change < control$tol) {
      converged <- TRUE
      break
    }
  }

  iterations <- iter
  ll.trace <- ll.trace[seq_len(iterations + 1L)]
  final.linear.inc <- drop(X.inc %*% gamma)
  final.cumulative.hazard <- drop(A %*% theta) + Y * drop(Z %*% beta)
  susceptible.weight <- ifelse(
    event,
    1,
    ifelse(beyond.support, 0, plogis(final.linear.inc - final.cumulative.hazard))
  )

  names(gamma) <- colnames(X.inc)
  names(beta) <- colnames(Z)
  names(theta) <- paste0("theta", seq_len(degree))

  ans <- list(
    call = match.call(),
    gamma = gamma,
    beta = beta,
    theta = theta,
    bernstein.coef = c(gamma0 = 0, cumsum(theta)),
    degree = degree,
    tau = tau,
    support.end = support.end,
    susceptible.probability = susceptible.weight,
    logLik = tail(ll.trace, 1),
    logLik.trace = ll.trace,
    converged = converged,
    iterations = iterations,
    max.change = max.change,
    finite.change = finite.change,
    baseline.change = baseline.change,
    monotonicity.violations = monotonicity.violations,
    n = n,
    n.event = sum(delta),
    data = list(
      Y = Y,
      delta = delta,
      Z = Z,
      X.inc = X.inc,
      A = A,
      a = a,
      beyond.support = beyond.support
    )
  )
  class(ans) <- "armc_em"
  ans
}

coef.armc_em <- function(object, ...) {
  setNames(
    c(unname(object$gamma), unname(object$beta)),
    c(
      paste0("gamma:", names(object$gamma)),
      paste0("beta:", names(object$beta))
    )
  )
}

logLik.armc_em <- function(object, ...) {
  out <- object$logLik
  attr(out, "df") <- length(object$gamma) + length(object$beta) + length(object$theta)
  attr(out, "nobs") <- object$n
  class(out) <- "logLik"
  out
}

print.armc_em <- function(x, ...) {
  cat("Right-censored additive-risk mixture-cure EM fit\n")
  cat(sprintf("n = %d; failures = %d; Bernstein degree = %d\n", x$n, x$n.event, x$degree))
  cat(sprintf(
    "converged = %s; iterations = %d; logLik = %.6f\n",
    x$converged, x$iterations, x$logLik
  ))
  cat("\nIncidence coefficients (gamma):\n")
  print(x$gamma)
  cat("\nLatency coefficients (beta):\n")
  print(x$beta)
  invisible(x)
}

predict.armc_em <- function(object,
                            newZ,
                            times,
                            newX.inc = NULL,
                            type = c(
                              "population_survival",
                              "susceptible_survival",
                              "incidence",
                              "cumulative_hazard",
                              "hazard"
                            ),
                            ...) {
  type <- match.arg(type)
  newZ <- as.matrix(newZ)
  times <- as.numeric(times)
  nr <- max(nrow(newZ), length(times))
  if (nrow(newZ) == 1L && nr > 1L) newZ <- newZ[rep(1L, nr), , drop = FALSE]
  if (length(times) == 1L && nr > 1L) times <- rep(times, nr)
  if (nrow(newZ) != nr || length(times) != nr) {
    stop("Rows of newZ and elements of times must be paired or recyclable from length one.")
  }
  if (ncol(newZ) != length(object$beta) || any(newZ < 0)) {
    stop("newZ has the wrong number of columns or contains negative values.")
  }
  if (any(times < 0)) stop("Prediction times must be nonnegative.")
  if (is.null(object$support.end) && any(times >= object$tau)) {
    stop("Without a finite support endpoint, prediction times must lie below fitted tau.")
  }
  if (is.null(newX.inc)) newX.inc <- cbind(`(Intercept)` = 1, newZ)
  newX.inc <- as.matrix(newX.inc)
  if (nrow(newX.inc) != nr || ncol(newX.inc) != length(object$gamma)) {
    stop("newX.inc has incompatible dimensions.")
  }

  beyond.support <- if (is.null(object$support.end)) {
    rep(FALSE, length(times))
  } else {
    times >= object$support.end
  }
  basis.times <- if (is.null(object$support.end)) times else pmin(times, object$support.end)
  basis <- bernstein_hazard_basis(basis.times, object$degree, object$tau)
  incidence <- plogis(drop(newX.inc %*% object$gamma))
  cumulative.hazard <- drop(basis$cumulative %*% object$theta) +
    times * drop(newZ %*% object$beta)
  hazard <- drop(basis$hazard %*% object$theta) + drop(newZ %*% object$beta)
  susceptible.survival <- exp(-cumulative.hazard)
  cumulative.hazard[beyond.support] <- Inf
  hazard[beyond.support] <- Inf
  susceptible.survival[beyond.support] <- 0

  switch(
    type,
    incidence = incidence,
    cumulative_hazard = cumulative.hazard,
    hazard = hazard,
    susceptible_survival = susceptible.survival,
    population_survival = 1 - incidence + incidence * susceptible.survival
  )
}

vcov_armc_em <- function(object, finite.only = TRUE, eigen.tol = 1e-8) {
  dat <- object$data
  q <- length(object$gamma)
  p <- length(object$beta)
  degree <- length(object$theta)
  par0 <- c(object$gamma, log(object$beta), log(object$theta))

  nll <- function(par) {
    gamma <- par[seq_len(q)]
    beta <- exp(par[q + seq_len(p)])
    theta <- exp(par[q + p + seq_len(degree)])
    -observed_loglik_armc(
      gamma, beta, theta,
      dat$Y, dat$delta, dat$Z, dat$X.inc, dat$A, dat$a,
      dat$beyond.support
    )
  }

  hessian <- optimHess(par0, nll)
  hessian <- (hessian + t(hessian)) / 2
  ee <- eigen(hessian, symmetric = TRUE)
  cutoff <- eigen.tol * max(abs(ee$values))
  if (any(!is.finite(ee$values)) || min(ee$values) <= cutoff) {
    warning("Observed information is singular or not positive definite; returning NA covariance.")
    d <- if (finite.only) q + p else q + p + degree
    out <- matrix(NA_real_, d, d)
    nm <- c(
      paste0("gamma:", names(object$gamma)),
      paste0("beta:", names(object$beta)),
      names(object$theta)
    )[seq_len(d)]
    dimnames(out) <- list(nm, nm)
    return(out)
  }

  covariance.logscale <- ee$vectors %*% diag(1 / ee$values) %*% t(ee$vectors)
  jacobian <- diag(c(rep(1, q), object$beta, object$theta))
  covariance <- jacobian %*% covariance.logscale %*% jacobian
  nm <- c(
    paste0("gamma:", names(object$gamma)),
    paste0("beta:", names(object$beta)),
    names(object$theta)
  )
  dimnames(covariance) <- list(nm, nm)
  if (finite.only) covariance <- covariance[seq_len(q + p), seq_len(q + p), drop = FALSE]
  covariance
}

fit_generated_data <- function(simulated,
                               degree = NULL,
                               compute.se = FALSE,
                               control = list()) {
  dat <- simulated$data
  fit <- em_armc_right(
    Y = dat$Y,
    delta = dat$delta,
    Z = as.matrix(dat[, c("Z1", "Z2")]),
    degree = degree,
    support.end = simulated$tau0,
    control = control
  )
  if (compute.se) {
    fit$vcov <- vcov_armc_em(fit)
    fit$se <- sqrt(diag(fit$vcov))
  }
  fit
}

run_armc_simulation <- function(n,
                                repetition,
                                censor.setting = c("L", "H"),
                                pai.setting = c("S1", "S2", "S3"),
                                degree = NULL,
                                seed.base = 20260921L,
                                compute.se = FALSE,
                                cores = 1L,
                                control = list()) {
  censor.setting <- match.arg(censor.setting)
  pai.setting <- match.arg(pai.setting)
  rep.ids <- seq_len(repetition)
  if (is.null(degree)) degree <- max(3L, floor(n^(1 / 4)))

  one.rep <- function(rep.id) {
    sim <- Get.random(
      n = n,
      censor.setting = censor.setting,
      pai.setting = pai.setting,
      seed = seed.base + rep.id
    )
    result <- tryCatch({
      fit <- fit_generated_data(sim, degree = degree, compute.se = compute.se, control = control)
      values <- c(
        setNames(fit$gamma, paste0("gamma", seq_along(fit$gamma) - 1L)),
        setNames(fit$beta, paste0("beta", seq_along(fit$beta)))
      )
      ses <- setNames(rep(NA_real_, length(values)), paste0("se_", names(values)))
      if (compute.se && !is.null(fit$se)) ses[] <- unname(fit$se[seq_along(values)])
      c(
        list(
          rep = rep.id,
          n = n,
          censor.setting = censor.setting,
          pai.setting = pai.setting,
          incidence.correctly.specified = pai.setting == "S1",
          censor.rate = sim$censor.rate,
          cure.rate = sim$cure.rate,
          converged = fit$converged,
          iterations = fit$iterations,
          logLik = fit$logLik,
          monotonicity.violations = fit$monotonicity.violations,
          error = NA_character_
        ),
        as.list(values),
        as.list(ses)
      )
    }, error = function(e) {
      value.names <- c(paste0("gamma", 0:2), paste0("beta", 1:2))
      c(
        list(
          rep = rep.id,
          n = n,
          censor.setting = censor.setting,
          pai.setting = pai.setting,
          incidence.correctly.specified = pai.setting == "S1",
          censor.rate = sim$censor.rate,
          cure.rate = sim$cure.rate,
          converged = FALSE,
          iterations = NA_integer_,
          logLik = NA_real_,
          monotonicity.violations = NA_integer_,
          error = conditionMessage(e)
        ),
        as.list(setNames(rep(NA_real_, length(value.names)), value.names)),
        as.list(setNames(rep(NA_real_, length(value.names)), paste0("se_", value.names)))
      )
    })
    as.data.frame(result, stringsAsFactors = FALSE, check.names = FALSE)
  }

  if (cores > 1L && .Platform$OS.type != "windows") {
    rows <- parallel::mclapply(rep.ids, one.rep, mc.cores = cores, mc.set.seed = FALSE)
  } else {
    rows <- lapply(rep.ids, one.rep)
  }
  do.call(rbind, rows)
}

summarise_armc_simulation <- function(results) {
  parameter.names <- c("gamma0", "gamma1", "gamma2", "beta1", "beta2")
  truth <- c(gamma0 = NA_real_, gamma1 = NA_real_, gamma2 = NA_real_, beta1 = 1, beta2 = 1)
  if (all(results$pai.setting == "S1")) {
    cs <- unique(results$censor.setting)
    if (length(cs) == 1L) truth[1:3] <- simulation_parameters(cs)$gamma.true
  }
  ok <- results$converged & is.na(results$error)
  out <- lapply(parameter.names, function(nm) {
    x <- results[[nm]][ok]
    data.frame(
      parameter = nm,
      truth = truth[nm],
      mean = if (length(x)) mean(x, na.rm = TRUE) else NA_real_,
      bias = if (length(x) && is.finite(truth[nm])) mean(x, na.rm = TRUE) - truth[nm] else NA_real_,
      SD = if (length(x) > 1L) sd(x, na.rm = TRUE) else NA_real_,
      valid_n = sum(is.finite(x)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

get_env_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  as.integer(value)
}

if ("--run" %in% commandArgs(trailingOnly = TRUE)) {
  censor.setting <- Sys.getenv("CENSOR_SETTING", unset = "L")
  pai.setting <- Sys.getenv("PAI_SETTING", unset = "S1")
  n <- get_env_integer("N", 300L)
  repetition <- get_env_integer("REPETITION", 1L)
  degree.env <- Sys.getenv("BERNSTEIN_DEGREE", unset = "")
  degree <- if (nzchar(degree.env)) as.integer(degree.env) else NULL
  seed.base <- get_env_integer("SEED_BASE", 20260921L)
  cores <- get_env_integer("CORES", 1L)
  compute.se <- identical(tolower(Sys.getenv("COMPUTE_SE", unset = "false")), "true")
  output.file <- Sys.getenv(
    "OUTPUT_FILE",
    unset = sprintf("armc_em_%s_%s_n%d.csv", pai.setting, censor.setting, n)
  )

  results <- run_armc_simulation(
    n = n,
    repetition = repetition,
    censor.setting = censor.setting,
    pai.setting = pai.setting,
    degree = degree,
    seed.base = seed.base,
    compute.se = compute.se,
    cores = cores
  )
  write.csv(results, output.file, row.names = FALSE)
  print(summarise_armc_simulation(results))
  message("Saved simulation results to: ", normalizePath(output.file, mustWork = FALSE))
}
