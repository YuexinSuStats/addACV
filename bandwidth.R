library(iterators)
library(foreach)
library(Rmpi)
library(doMPI)
library(survival)
library(MASS)


# Scenario 1 with n = 300 under moderate and high cure rates.
n <- 300L
repetition <- 500L
rep.ids <- seq_len(repetition)
SEED_INTERNAL_DATA <- 3L
beta.true <- c(1, 1)

scenario.settings <- list(
  S1_moderate_cure_n300 = list(
    lambdaC = 0.50,
    gamma.true = c(-1.5, 1, 1)
  ),
  S1_high_cure_n300 = list(
    lambdaC = 1.00,
    gamma.true = c(-2.5, 1, 1)
  )
)

bandwidth_grid <- seq(0.1, 2.0, by = 0.05)

pai <- function(u) plogis(u)

tau0 <- 2
kappa <- 1

H_0 <- function(t) {
  ans <- rep(Inf, length(t))
  id <- t < tau0
  ans[id] <- 2 * log1p(t[id]) -
    kappa * log1p(-t[id] / tau0)
  ans
}

Get.random <- function(n, lambdaC, gamma.true) {
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

Get.Data <- function(rep.ids, lambdaC, gamma.true){
  data <- vector("list", repetition)
  censor.rate <- cure.rate <- numeric(repetition)
  for (i in seq_along(rep.ids)) {
    set.seed(SEED_INTERNAL_DATA + rep.ids[i])
    tmp <- Get.random(n, lambdaC = lambdaC, gamma.true = gamma.true)
    data[[i]] <- tmp$data
    censor.rate[i] <- tmp$censor.rate
    cure.rate[i] <- tmp$cure.rate
  }
  list(data = data, censor.rate = censor.rate, cure.rate = cure.rate)
}

KK <- function(x, hn) 3/4 * (1 - (x/hn)^2) * (abs(x/hn) <= 1)

NWweight2 <- function(Z, Z.sample, hn, lam = rep(1, nrow(Z.sample))){
  m <- dim(Z)[1]
  p1 <- dim(Z)[2]
  nn <- dim(Z.sample)[1]
  Z <- as.data.frame(Z)
  Z.sample <- as.data.frame(Z.sample)
  MZ <- array(apply(Z, 2, function(x) rep(x, nn)), dim = c(m, nn, p1))
  tMZ <- aperm(
    array(apply(Z.sample, 2, function(x) rep(x, m)), dim = c(nn, m, p1)),
    c(2, 1, 3)
  )
  dMZ <- MZ - tMZ
  Kp <- array(NA, dim = c(m, nn, p1))
  # addICV.R uses two continuous covariates and the Epanechnikov kernel for
  # both components.
  Kp[,,1] <- KK(dMZ[,,1], hn[1])
  Kp[,,2] <- KK(dMZ[,,2], hn[2])
  prodKp <- apply(Kp, c(1, 2), prod)
  Lam <- matrix(lam, nrow = m, ncol = nn, byrow = TRUE)
  lamKp <- prodKp * Lam
  rs <- matrix(apply(lamKp, 1, sum), nrow = m, ncol = nn, byrow = FALSE)
  rs[rs == 0] <- 1e-10
  lamKp / rs
}

S.Ti <- function(data, w){
  m <- nrow(w)
  nn <- ncol(w)
  k <- cbind(rep(0, m), w[,1:(nn - 1)])
  cumsum.k <- t(apply(k, 1, cumsum))
  cum.w <- 1 - cumsum.k
  cum.w <- cum.w * (cum.w > w + 1e-5) + (w + 1e-4) * (cum.w <= w + 1e-5)
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

P.T <- function(data, w){
  m <- nrow(w)
  nn <- ncol(w)
  k <- cbind(rep(0, m), w[,1:(nn - 1)])
  cumsum.k <- t(apply(k, 1, cumsum))
  cum.w <- 1 - cumsum.k
  w <- replace(w, cum.w == 0, 0)
  cum.w <- replace(cum.w, cum.w == 0, 1)
  q <- 1 - w / cum.w
  delta <- matrix(data$delta, nrow = m, ncol = nn, byrow = TRUE)
  q <- replace(q, delta == 0, 1)
  q <- t(apply(q, 1, prod))
  as.numeric(replace(q, q <= 0, 0))
}

get.est <- function(data, bandwidth.c){
  data <- data[order(data[,1]), ]
  Y <- data$Y
  delta <- data$delta
  Z <- as.matrix(data[,-(1:2)])
  p <- ncol(Z)
  nn <- nrow(data)
  T1.max <- max(Y[delta == 1])
  # Same component-specific bandwidth and WEE construction as addICV.R;
  # bandwidth.c is varied only for the sensitivity curve.
  h <- bandwidth.c * nn^(-2/7) * apply(Z, 2L, stats::sd)
  if (any(!is.finite(h)) || any(h <= 0)) {
    stop("All bandwidths must be finite and strictly positive.")
  }

  wght <- NWweight2(Z = Z, Z.sample = Z, hn = h)
  S.beran <- S.Ti(data, wght)
  P.T_inf <- P.T(data, wght)
  S.safe <- S.beran
  S.safe[S.safe == 0] <- 1e-8
  case.wght <- delta +
    (1 - delta) * (S.beran - P.T_inf) / S.safe * (Y <= T1.max)

  dif <- diff(c(0, Y))
  rle.Y <- rle(Y)
  start.idx <- c(1, 1 + head(cumsum(rle.Y$lengths), -1))
  first.indices <- rep(start.idx, rle.Y$lengths)

  Zbarr <- matrix(0, nn, p)
  for (j in seq_len(nn)) {
    idx <- first.indices[j]:nn
    wj <- case.wght[idx]
    sw <- sum(wj)
    if (sw > 0) {
      Zbarr[j, ] <- colSums(Z[idx, , drop = FALSE] * wj) / sw
    }
  }

  fenzi <- colSums(case.wght * delta * (Z - Zbarr))
  fenmu <- matrix(0, p, p)
  for (j in seq_len(nn)) {
    if (dif[j] > 1e-8) {
      idx <- first.indices[j]:nn
      wj <- case.wght[idx]
      Zc <- sweep(Z[idx, , drop = FALSE], 2, Zbarr[j, ], "-")
      fenmu <- fenmu +
        (t(Zc) %*% sweep(Zc, 1, wj, "*")) * dif[j]
    }
  }

  beta.WEE <- as.vector(fenzi %*% solve(fenmu))
  if (any(!is.finite(beta.WEE))) {
    stop("Non-finite beta estimate.")
  }
  beta.WEE
}

Data <- lapply(
  scenario.settings,
  function(setting) {
    Get.Data(
      rep.ids = rep.ids,
      lambdaC = setting$lambdaC,
      gamma.true = setting$gamma.true
    )
  }
)
scenario.names <- names(Data)

packagelist <- c("survival", "MASS")

exports <- c(
  "get.est", "NWweight2", "S.Ti", "P.T", "KK",
  "Data", "scenario.names", "beta.true", "repetition"
)

# Draw once and export through two vector devices. The dimensions correspond
# to Wiley's large-figure width (180 mm).
draw_mse_plot <- function(plot.data) {
  par(mar = c(4.2, 4.2, 1, 1), mgp = c(2.4, 0.7, 0),
      cex.lab = 1.0, cex.axis = 0.9)
  plot(
    plot.data$bandwidth.c,
    plot.data$MSE_beta1,
    type = "l",
    lty = 1,
    lwd = 1.5,
    xlab = "Bandwidth coefficient",
    ylab = "MSE",
    ylim = c(
      0,
      max(plot.data$MSE_beta1, plot.data$MSE_beta2, na.rm = TRUE) * 1.1
    )
  )
  lines(plot.data$bandwidth.c, plot.data$MSE_beta2, lty = 2, lwd = 1.5)
  legend(
    "topright",
    legend = c(
      expression("Estimate " * hat(beta)[1]),
      expression("Estimate " * hat(beta)[2])
    ),
    lty = c(1, 2),
    lwd = 1.5,
    bty = "n",
    cex = 0.9
  )
}

figure_width_in <- 180 / 25.4
figure_height_in <- 4.2

# Use the same Rmpi/doMPI backend as addICV.R. The scheduler determines the
# number of MPI workers made available to startMPIcluster().
cl <- startMPIcluster()
registerDoMPI(cl)

cat("MPI workers:", foreach::getDoParWorkers(), "\n")

tryCatch({
  results.list <- foreach(
    bandwidth.c = bandwidth_grid,
    .packages = packagelist,
    .export = exports,
    .errorhandling = "pass",
    .inorder = TRUE
  ) %dopar% {
    scenario.results <- lapply(scenario.names, function(scenario.name) {
      est.mat <- matrix(NA_real_, nrow = repetition, ncol = length(beta.true))
      for (r in seq_len(repetition)) {
        fit <- tryCatch(get.est(Data[[scenario.name]]$data[[r]], bandwidth.c), error = function(e) rep(NA_real_, length(beta.true)))
        if (length(fit) == length(beta.true) && all(is.finite(fit))) est.mat[r, ] <- fit
      }
      data.frame(
        Scenario = scenario.name,
        bandwidth.c = bandwidth.c,
        MSE_beta1 = mean((est.mat[,1] - beta.true[1])^2, na.rm = TRUE),
        MSE_beta2 = mean((est.mat[,2] - beta.true[2])^2, na.rm = TRUE),
        n_fail = sum(!complete.cases(est.mat)),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, scenario.results)
  }
  results.df <- do.call(rbind, results.list)
  rownames(results.df) <- NULL
  print(results.df)
  for (scenario.name in scenario.names) {
    cat("Scenario:", scenario.name, "\n")
    cat("Censor rate:", round(mean(Data[[scenario.name]]$censor.rate), 3), "\n")
    cat("Cure rate:", round(mean(Data[[scenario.name]]$cure.rate), 3), "\n")
  }

  # Generate all figures before mpi.quit(), which terminates the master R
  # process after the MPI workers have been closed.
  for (scenario.name in scenario.names) {
    plot.data <- results.df[
      results.df$Scenario == scenario.name,
      ,
      drop = FALSE
    ]

    grDevices::cairo_pdf(
      paste0("MSEplot_", scenario.name, ".pdf"),
      width = figure_width_in,
      height = figure_height_in,
      pointsize = 10,
      family = "sans",
      bg = "white",
      fallback_resolution = 600
    )
    draw_mse_plot(plot.data)
    dev.off()

    grDevices::cairo_ps(
      paste0("MSEplot_", scenario.name, ".eps"),
      width = figure_width_in,
      height = figure_height_in,
      pointsize = 10,
      family = "sans",
      bg = "white",
      onefile = FALSE,
      fallback_resolution = 600
    )
    draw_mse_plot(plot.data)
    dev.off()
  }
}, finally = {
  closeCluster(cl)
  Rmpi::mpi.quit()
})

