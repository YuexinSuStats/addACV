library(iterators)
library(foreach)
library(doParallel)
library(parallel)
library(survival)
library(MASS)


set.seed(1)

beta.true <- c(1, 1)

lambdaC <- 0.50
gamma.true <- c(-1.5, 1, 1)

n <- 300
repetition <- 400
bandwidth_grid <- seq(0.1, 2.0, by = 0.05)

pai <- function(u){
  g <- exp(u) / (1 + exp(u))
  return(g)
}

H_0 <- function(t) 2 * log(t + 1)

Get.random <- function(n) {
  Z <- matrix(1, nrow = n, ncol = length(beta.true))
  Z[,1] <- runif(n, 0, 2)
  Z[,2] <- runif(n, 0, 4)

  T0 <- numeric(n)
  u <- runif(n, 0, 1)
  for (i in 1:n) {
    S0 <- function(t) exp(-H_0(t) - as.numeric(Z[i,] %*% beta.true) * t)
    T0[i] <- uniroot(function(t) 1 - S0(t) - u[i], c(0, 100))$root
  }

  Z.multi <- cbind(1, Z)
  peta <- pai(Z.multi %*% gamma.true)
  eta <- rbinom(n, 1, peta)
  TT <- eta * T0 + (1 - eta) * 100

  C <- rexp(n, lambdaC)
  Y <- pmin(TT, C)
  delta <- as.numeric(TT <= C)

  data <- data.frame(Y = Y, delta = delta, Z = Z)
  data <- data[order(data[,1]), ]
  list(data = data, censor.rate = mean(delta == 0), cure.rate = mean(eta == 0))
}

Get.Data <- function(){
  data <- vector("list", repetition)
  censor.rate <- cure.rate <- numeric(repetition)
  for (i in seq_len(repetition)) {
    tmp <- Get.random(n)
    data[[i]] <- tmp$data
    censor.rate[i] <- tmp$censor.rate
    cure.rate[i] <- tmp$cure.rate
  }
  list(data = data, censor.rate = censor.rate, cure.rate = cure.rate)
}

KK <- function(x, hn) 3/4 * (1 - (x/hn)^2) * (abs(x/hn) <= 1)

KK2 <- function(x, hn) hn^(ifelse(x != 0, 1, 0))

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
  for (i in seq_len(p1)) {
    if (sum(dMZ[,,i] == 1 | dMZ[,,i] == -1 | dMZ[,,i] == 0) == (nn * nn) &&
        any(dMZ[,,i] != 0)) {
      Kp[,,i] <- KK2(dMZ[,,i], hn[i])
    } else {
      Kp[,,i] <- KK(dMZ[,,i], hn[i])
    }
  }
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
  sd.Z <- apply(Z, 2, sd)
  h <- bandwidth.c * nn^(-2/7) * sd.Z

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

Data <- Get.Data()

packagelist <- c("survival", "MASS")

exports <- c(
  "get.est", "NWweight2", "S.Ti", "P.T", "KK", "KK2",
  "Data", "beta.true", "repetition"
)

num.cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)

tryCatch({
  results.list <- foreach(bandwidth.c = bandwidth_grid, .packages = packagelist, .export = exports, .errorhandling = "pass", .inorder = TRUE) %dopar% {
    est.mat <- matrix(NA_real_, nrow = repetition, ncol = length(beta.true))
    for (r in seq_len(repetition)) {
      fit <- tryCatch(get.est(Data$data[[r]], bandwidth.c), error = function(e) rep(NA_real_, length(beta.true)))
      if (length(fit) == length(beta.true) && all(is.finite(fit))) est.mat[r, ] <- fit
    }
    c(MSE_beta1 = mean((est.mat[,1] - beta.true[1])^2, na.rm = TRUE), MSE_beta2 = mean((est.mat[,2] - beta.true[2])^2, na.rm = TRUE), n_fail = sum(!complete.cases(est.mat)))
  }
  results <- do.call(rbind, results.list)
  rownames(results) <- NULL
  results.df <- data.frame(bandwidth.c = bandwidth_grid, results)
  print(round(results.df, 4))
  cat("Censor rate:", round(mean(Data$censor.rate), 3), "\n")
  cat("Cure rate:", round(mean(Data$cure.rate), 3), "\n")
}, finally = {
  parallel::stopCluster(cl)
})

png("MSEplot.png", width = 10, height = 6, units = "in", res = 500)
par(mar = c(4.2, 4.2, 1, 1))
plot(results.df$bandwidth.c, results.df$MSE_beta1, type = "l", lty = 1, lwd = 1.5, xlab = "Bandwidth coefficient", ylab = "MSE", ylim = c(0, max(results.df$MSE_beta1, results.df$MSE_beta2, na.rm = TRUE) * 1.1))
lines(results.df$bandwidth.c, results.df$MSE_beta2, lty = 2, lwd = 1.5)
legend("topright", legend = c(expression("Estimate " * hat(beta)[1]), expression("Estimate " * hat(beta)[2])), lty = c(1, 2), lwd = 1.5, bty = "n")
dev.off()

