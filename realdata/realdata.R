library(miscTools)
library(survival)
library(ggplot2)
library(gridExtra)
library(survminer)
library(MASS)
library(mvtnorm)
library(lars)
library(Matrix)



###prepare data---------------


load("clinical.RData")
# - a preparation
sdata.pre <- data.frame(
  yobs   = ifelse(!is.na(clinical[,'days_to_death']),
                  clinical[,'days_to_death'],
                  clinical[,'days_to_last_follow_up']) / 365,
  delta  = ifelse(!is.na(clinical[,'days_to_death']), 1, 0),
  
  Age    = clinical[,'age_at_diagnosis'] / 365,
  
  Race   = ifelse(clinical[,'race'] == 'black or african american', 'black',
                  ifelse(clinical[,'race'] == 'white', 'white', 'other')),
  
  Gender = ifelse(clinical[,'sex_at_birth'] == 'female', 'Female',
                  ifelse(clinical[,'sex_at_birth'] == 'male', 'Male', NA)),
  Laterality = ifelse(clinical[,'laterality'] == 'Left','Left','Right'),
  Stage  = sapply(
    clinical[,'ajcc_pathologic_stage'],
    function(x, pattern = 'Stage X|Stage IV|Stage [I]*') {
      ifelse(grepl(pattern, x),
             regmatches(x, regexpr(pattern, x)),
             NA)
    },
    USE.NAMES = FALSE
  ),
  
  ERstatus = clinical[,'ER_Status'],
  Day = clinical[,'days_to_last_follow_up']
  
)


# - extract covariates and remove undesiable subjects and NA
sdata.TCGA <- na.omit(
  sdata.pre[
    sdata.pre[,'yobs'] > 0
    & sdata.pre[,'Age'] <= 75
    & sdata.pre[,'Gender'] == "Female"
    & sdata.pre[,'Race'] %in% c('white')
    & sdata.pre[,'Laterality'] %in% c('Right')
    & sdata.pre[,'Stage'] %in% c('Stage I','Stage II','Stage III')
    & sdata.pre[,'Day'] > 0 ,
    c('yobs','delta','Age','ERstatus'),
  ]
)



####get small sample
rownames(sdata.TCGA) <- NULL
cens.Rate<-sum(sdata.TCGA$delta==0)/nrow(sdata.TCGA)
##cure test
Tmax<-max(sdata.TCGA$yobs)
T1max<-max(sdata.TCGA[sdata.TCGA$delta==1,1])
n<-nrow(sdata.TCGA)
subsetdata <- sdata.TCGA[
  sdata.TCGA$yobs > (2*T1max - Tmax) &
    sdata.TCGA$yobs <= T1max &
    sdata.TCGA$delta == 1,
]
curetest<-(1-nrow(subsetdata)/n)^n

maller_zhou_cure <- function(data) {
  Y     <- data$yobs
  delta <- data$delta
  
  
  km_fit <- survfit(Surv(Y, delta) ~ 1)
  
  cure_rate <- min(km_fit$surv)
  
  return(cure_rate)
}

cure.rate <- maller_zhou_cure(sdata.TCGA)
sum<-summary(sdata.TCGA)
# auxiliary information ---------------------------------------------------


# - rescale the Age variable and Stage variable
Age.Min <- min(sdata.TCGA$Age); Age.Max <- max(sdata.TCGA$Age)
sdata.TCGA$Age <- (sdata.TCGA$Age-Age.Min)/(Age.Max-Age.Min)


## fit the model with auxiliary information

# - prepare the auxiliary information
Age.Cut <- c(0,(c(45,65)-Age.Min)/(Age.Max-Age.Min),1)

##### calculate subgroup survival rates
k.num<-6; m.num<-2
q.num<-k.num*m.num
tstar<-c(5,10)
bandwidth.c <- 3
tpro<-matrix(c(0.934,0.869,0.769,0.788,0.724,0.653,
               0.735,0.661,0.701,0.716,0.635,0.580),k.num,m.num)

###auxiliary information tpro k*m,k for subgroup,m for time point
subgf<- function(X){
  rbind((X[,1] >= Age.Cut[1] & X[,1] <  Age.Cut[2] & X[,2] == 1),
        (X[,1] >= Age.Cut[2] & X[,1] <  Age.Cut[3] & X[,2] == 1),
        (X[,1] >= Age.Cut[3] & X[,1] <  Age.Cut[4] & X[,2] == 1),
        (X[,1] >= Age.Cut[1] & X[,1] <  Age.Cut[2] & X[,2] == 0),
        (X[,1] >= Age.Cut[2] & X[,1] <  Age.Cut[3] & X[,2] == 0),
        (X[,1] >= Age.Cut[3] & X[,1] <  Age.Cut[4] & X[,2] == 0))}

getprob<-function(data,tstar){
  pro <- matrix(0, nrow = k.num, ncol = m.num)
  group<-subgf(as.matrix(data[,-(1:2)]))
  for (i in 1:k.num) {
    group_indices <- which(group[i, ])
    group_data <- data[group_indices, ]
    time<-group_data[,1]
    status<-group_data[,2]
    fit <- survival::survfit(survival::Surv(time, status) ~ 1)##K_M method
    pro[i,] <- summary(fit, times = tstar, extend = TRUE)$surv
  }
  return(pro)
}

# Estimation ---------------------------------------------------------------

N_ex <- 400L
MIN_VALID_BOOT <- 30L
COV_EIGEN_FLOOR <- 1e-8

KK <- function(x, hn) 3/4 * (1 - (x/hn)^2) * (abs(x/hn) <= 1)
KK2 <- function(x, hn) hn^(ifelse(x != 0, 1, 0))

NWweight2 <- function(Z, Z.sample, hn, lam = rep(1, nrow(Z.sample))){
  m <- nrow(Z); p1 <- ncol(Z); nn <- nrow(Z.sample)
  Z <- as.data.frame(Z); Z.sample <- as.data.frame(Z.sample)
  MZ <- array(apply(Z, 2, function(x) rep(x, nn)), dim = c(m, nn, p1))
  tMZ <- aperm(array(apply(Z.sample, 2, function(x) rep(x, m)), dim = c(nn, m, p1)), c(2, 1, 3))
  dMZ <- MZ - tMZ
  Kp <- array(NA, dim = c(m, nn, p1))
  for (i in seq_len(p1)) {
    is.binary <- all(dMZ[,,i] %in% c(-1, 0, 1)) && any(dMZ[,,i] != 0)
    if (is.binary) Kp[,,i] <- KK2(dMZ[,,i], hn[i])
    else Kp[,,i] <- KK(dMZ[,,i], hn[i])
  }
  prodKp <- apply(Kp, c(1, 2), prod)
  lamKp <- prodKp * matrix(lam, nrow = m, ncol = nn, byrow = TRUE)
  rs <- matrix(rowSums(lamKp), nrow = m, ncol = nn, byrow = FALSE)
  rs[rs == 0] <- 1e-10
  lamKp / rs
}

S.Ti <- function(data, w){
  m <- nrow(w); nn <- ncol(w)
  k <- cbind(rep(0, m), w[, seq_len(nn - 1), drop = FALSE])
  cum.w <- 1 - t(apply(k, 1, cumsum))
  cum.w <- cum.w * (cum.w > w + 1e-5) + (w + 1e-4) * (cum.w <= w + 1e-5)
  q <- 1 - w / cum.w
  delta <- matrix(data$delta, nrow = m, ncol = nn, byrow = TRUE)
  q <- replace(q, delta == 0, 1)
  q <- t(apply(q, 1, cumprod))
  r <- rle(data$yobs)
  end.index <- cumsum(r$lengths)
  tied.end <- rep(end.index, r$lengths)
  q[cbind(seq_len(m), tied.end[seq_len(m)])]
}

P.T <- function(data, w){
  m <- nrow(w); nn <- ncol(w)
  k <- cbind(rep(0, m), w[, seq_len(nn - 1), drop = FALSE])
  cum.w <- 1 - t(apply(k, 1, cumsum))
  w <- replace(w, cum.w == 0, 0)
  cum.w <- replace(cum.w, cum.w == 0, 1)
  q <- 1 - w / cum.w
  delta <- matrix(data$delta, nrow = m, ncol = nn, byrow = TRUE)
  q <- replace(q, delta == 0, 1)
  q <- t(apply(q, 1, prod))
  as.numeric(replace(q, q <= 0, 0))
}

get.est <- function(data, tstar, tpro){
  data <- data[order(data$yobs), ]
  Y <- data$yobs
  delta <- data$delta
  Z <- model.matrix(~ . - 1, data = data[, c("Age", "ERstatus")])
  p <- ncol(Z); nn <- nrow(data)
  T1.max <- max(Y[delta == 1])
  h <- c(bandwidth.c * nn^(-2/7) * sd(Z[,1]), 0.01)
  
  wght <- NWweight2(Z, Z, h)
  S.beran <- S.Ti(data, wght)
  P.T_inf <- P.T(data, wght)
  S.safe <- S.beran
  S.safe[S.safe == 0] <- 1e-8
  case.wght <- delta + (1 - delta) * (S.beran - P.T_inf) / S.safe * (Y <= T1.max)
  
  dif <- diff(c(0, Y))
  rle.Y <- rle(Y)
  start.index <- c(1, 1 + head(cumsum(rle.Y$lengths), -1))
  tied.start <- rep(start.index, rle.Y$lengths)
  
  Zbarr <- matrix(0, nn, p)
  for (j in seq_len(nn)) {
    risk.index <- tied.start[j]:nn
    risk.weight <- case.wght[risk.index]
    sum.weight <- sum(risk.weight)
    if (sum.weight > 0) {
      Zbarr[j, ] <- colSums(Z[risk.index, , drop = FALSE] * risk.weight) / sum.weight
    }
  }
  
  fenzi <- colSums(case.wght * delta * (Z - Zbarr))
  fenmu <- matrix(0, p, p)
  for (j in seq_len(nn)) {
    if (dif[j] > 1e-8) {
      risk.index <- tied.start[j]:nn
      risk.weight <- case.wght[risk.index]
      Z.centered <- sweep(Z[risk.index, , drop = FALSE], 2, Zbarr[j, ], "-")
      fenmu <- fenmu + t(Z.centered) %*% sweep(Z.centered, 1, risk.weight, "*") * dif[j]
    }
  }
  
  beta.WEE <- as.vector(fenzi %*% solve(fenmu))
  if (any(!is.finite(beta.WEE))) return(list(ok = FALSE, reason = "non-finite beta"))
  
  lambest <- function(ts){
    eps <- 1e-12
    risk.weight <- rev(cumsum(rev(case.wght)))
    zbeta <- as.vector(Z %*% beta.WEE)
    risk.weight.z <- rev(cumsum(rev(case.wght * zbeta)))
    
    event.number <- delta * case.wght * (Y <= ts)
    event.risk <- risk.weight[tied.start]
    if (any(event.number != 0 & event.risk <= eps)) {
      stop("Zero risk-set denominator at an event time.")
    }
    term1 <- numeric(nn)
    valid.event <- event.risk > eps
    term1[valid.event] <- event.number[valid.event] / event.risk[valid.event]
    result1 <- sum(term1)
    
    point <- findInterval(ts, Y)
    if (point == 0L) {
      if (risk.weight[1] <= eps) stop("Zero initial risk-set denominator.")
      result2 <- ts * risk.weight.z[1] / risk.weight[1]
    } else {
      index <- seq_len(point)
      if (any(dif[index] > 0 & risk.weight[index] <= eps)) {
        stop("Zero risk-set denominator in baseline-hazard integral.")
      }
      valid.interval <- dif[index] > 0 & risk.weight[index] > eps
      result2 <- sum(dif[index][valid.interval] * risk.weight.z[index][valid.interval] /
                       risk.weight[index][valid.interval])
      if (point < nn && ts > Y[point]) {
        if (risk.weight[point + 1L] <= eps) {
          stop("Zero risk-set denominator after the last observed time.")
        }
        result2 <- result2 + (ts - Y[point]) * risk.weight.z[point + 1L] /
          risk.weight[point + 1L]
      }
    }
    
    result <- result1 - result2
    if (!is.finite(result)) stop("Non-finite baseline cumulative hazard.")
    result
  }
  
  alpha <- tryCatch(vapply(tstar, lambest, numeric(1)), error = function(e) NULL)
  if (is.null(alpha) || any(!is.finite(alpha))) {
    return(list(ok = FALSE, reason = "invalid alpha"))
  }
  
  phi <- matrix(0, k.num, m.num)
  subdata <- cbind(1 - P.T_inf, Z)
  group <- subgf(Z)
  pai.hat <- apply(group, 1, mean)
  for (m in seq_len(m.num)) {
    for (k in seq_len(k.num)) {
      group.index <- which(group[k, ])
      if (length(group.index) > 0) {
        group.data <- subdata[group.index, , drop = FALSE]
        result <- 0
        for (i in seq_len(nrow(group.data))) {
          result <- result + (1 - group.data[i,1]) + group.data[i,1] *
            exp(-alpha[m] - as.numeric(beta.WEE %*% group.data[i,-1]) * tstar[m]) - tpro[k,m]
        }
        phi[k,m] <- result / nn
      }
    }
  }
  
  est.phi <- as.vector(phi)
  if (any(!is.finite(est.phi))) return(list(ok = FALSE, reason = "non-finite phi"))
  list(ok = TRUE, est.beta = beta.WEE, est.phi = est.phi, pai = pai.hat)
}

boot.cov <- function(data, tstar, tpro, B){
  nn <- nrow(data); p <- 2; q <- q.num
  beta.boot <- matrix(NA_real_, B, p)
  phi.boot <- matrix(NA_real_, B, q)
  km.boot <- matrix(NA_real_, B, q)
  
  for (b in seq_len(B)) {
    bootstrap.sample <- data[sample.int(nn, nn, replace = TRUE), ]
    fit <- tryCatch(get.est(bootstrap.sample, tstar, tpro), error = function(e) list(ok = FALSE))
    if (!isTRUE(fit$ok)) next
    beta.boot[b, ] <- fit$est.beta
    phi.boot[b, ] <- fit$est.phi
    km.boot[b, ] <- tryCatch(
      as.vector(getprob(bootstrap.sample, tstar)),
      error = function(e) rep(NA_real_, q)
    )
  }
  
  ok.beta <- complete.cases(beta.boot)
  ok.icv <- ok.beta & complete.cases(phi.boot) & complete.cases(km.boot)
  list(beta = beta.boot,
       phi = phi.boot,
       km = km.boot,
       ok.beta = ok.beta,
       ok.icv = ok.icv,
       nvalid = sum(ok.icv), ndrop = B - sum(ok.icv))
}

symmetrize <- function(A) (A + t(A)) / 2

regularized_inverse <- function(A, rel.floor = COV_EIGEN_FLOOR) {
  A <- symmetrize(A)
  ee <- eigen(A, symmetric = TRUE)
  scale <- max(abs(ee$values), 1)
  floor.value <- rel.floor * scale
  values.used <- pmax(ee$values, floor.value)
  A.inv <- ee$vectors %*%
    diag(1 / values.used, nrow = length(values.used)) %*%
    t(ee$vectors)
  symmetrize(A.inv)
}

acv.fit <- function(beta, phi, pai.q, V.inv, V.inv.root, V12, nn){
  if (length(phi) != length(pai.q)) {
    stop("phi and pai.q have incompatible dimensions.")
  }
  if (any(!is.finite(phi)) || any(!is.finite(pai.q)) || any(pai.q <= 0)) {
    stop("Non-finite auxiliary quantities or an empty subgroup.")
  }
  
  y.tilde <- as.vector(V.inv.root %*% phi)
  X.tilde <- V.inv.root %*% diag(pai.q, length(pai.q), length(pai.q))
  penalty.weight <- pai.q / pmax(abs(phi), .Machine$double.eps)
  penalty.weight[!is.finite(penalty.weight)] <- 1e8
  penalty.weight <- pmin(penalty.weight, 1e8)
  X.star <- sweep(X.tilde, 2L, penalty.weight, "/")
  
  fit <- lars::lars(X.star, y.tilde, type = "lasso", normalize = FALSE,
                    intercept = FALSE, trace = FALSE)
  tau.path <- t(as.matrix(fit$beta))
  tau.path <- sweep(tau.path, 1L, penalty.weight, "/")
  fitted.path <- X.tilde %*% tau.path
  RSS <- colSums(sweep(fitted.path, 1L, y.tilde, "-")^2)
  cardinality <- colSums(abs(tau.path) > 1e-10)
  IC <- RSS + log(nn) / nn * cardinality
  selected.step <- which.min(IC)
  tau <- tau.path[, selected.step]
  active.set <- which(abs(tau) > 1e-10)
  
  beta.ACV <- beta - as.vector(t(V12) %*% V.inv %*% (phi - pai.q * tau))
  list(beta = beta.ACV, tau = tau, active.set = active.set)
}

getest <- function(sdata.TCGA, tpro, B, tstar, k.num, q.num){
  fit0 <- tryCatch(get.est(sdata.TCGA, tstar, tpro), error = function(e) list(ok = FALSE))
  if (!isTRUE(fit0$ok)) stop("The WEE fit failed.")
  
  beta.WEE <- fit0$est.beta
  phi.est <- fit0$est.phi
  pai.hat <- fit0$pai
  nn <- nrow(sdata.TCGA)
  p <- length(beta.WEE)
  
  bootstrap <- boot.cov(sdata.TCGA, tstar, tpro, B)
  if (sum(bootstrap$ok.beta) < MIN_VALID_BOOT ||
      bootstrap$nvalid < MIN_VALID_BOOT) {
    stop("Fewer than 30 valid bootstrap samples.")
  }
  
  beta.boot.WEE <- bootstrap$beta[bootstrap$ok.beta, , drop = FALSE]
  V11.WEE <- cov(beta.boot.WEE * sqrt(nn))
  se.WEE <- sqrt(pmax(diag(V11.WEE / nn), 0))
  
  ok <- bootstrap$ok.icv
  beta.boot <- bootstrap$beta[ok, , drop = FALSE]
  phi.boot <- bootstrap$phi[ok, , drop = FALSE]
  km.boot <- bootstrap$km[ok, , drop = FALSE]
  V11 <- cov(beta.boot * sqrt(nn))
  V22 <- symmetrize(cov(phi.boot * sqrt(nn)))
  V12 <- cov(phi.boot * sqrt(nn), beta.boot * sqrt(nn))
  Sigma.phi <- symmetrize(cov(km.boot * sqrt(nn)))
  
  pai.q <- rep(pai.hat, times = m.num)
  if (length(phi.est) != q.num || length(pai.q) != q.num) {
    stop("The auxiliary vector must have dimension k.num * m.num.")
  }
  D.pi <- diag(pai.q, nrow = q.num)
  V22.U <- symmetrize(V22)
  
  finite.matrices <-
    all(is.finite(V11.WEE)) &&
    all(is.finite(V11)) &&
    all(is.finite(V12)) &&
    all(is.finite(V22)) &&
    all(is.finite(Sigma.phi)) &&
    all(is.finite(V22.U))
  if (!finite.matrices) stop("Non-finite bootstrap covariance matrix.")
  
  VU.inv <- regularized_inverse(V22.U)
  beta.ICV <- as.vector(beta.WEE - t(V12) %*% VU.inv %*% phi.est)
  var.ICV <- symmetrize((V11 - t(V12) %*% VU.inv %*% V12) / nn)
  if (any(!is.finite(beta.ICV)) || any(!is.finite(var.ICV)) ||
      any(diag(var.ICV) < -1e-10)) {
    stop("Invalid ICV estimate or variance.")
  }
  se.ICV <- sqrt(pmax(diag(var.ICV), 0))
  
  ee.inv <- eigen(VU.inv, symmetric = TRUE)
  V.inv.root <- ee.inv$vectors %*%
    diag(sqrt(pmax(ee.inv$values, 0)), nrow = q.num) %*%
    t(ee.inv$vectors)
  fit.ACV <- acv.fit(beta.WEE, phi.est, pai.q, VU.inv, V.inv.root, V12, nn)
  beta.ACV <- fit.ACV$beta
  inactive.set <- setdiff(seq_len(q.num), fit.ACV$active.set)
  if (length(inactive.set) > 0L) {
    Gamma2.hat <- V12[inactive.set, , drop = FALSE]
    V22.U.hat <- V22.U[inactive.set, inactive.set, drop = FALSE]
    inv.V22.U.hat <- regularized_inverse(V22.U.hat)
    var.ACV <- symmetrize(
      (V11 - t(Gamma2.hat) %*% inv.V22.U.hat %*% Gamma2.hat) / nn
    )
  } else {
    var.ACV <- symmetrize(V11 / nn)
  }
  if (any(!is.finite(beta.ACV)) || any(!is.finite(var.ACV)) ||
      any(diag(var.ACV) < -1e-10)) {
    stop("Invalid ACV estimate or variance.")
  }
  se.ACV <- sqrt(pmax(diag(var.ACV), 0))
  
  z.critical <- qnorm(0.975)
  CI.WEE <- cbind(beta.WEE - z.critical * se.WEE, beta.WEE + z.critical * se.WEE)
  CI.ICV <- cbind(beta.ICV - z.critical * se.ICV, beta.ICV + z.critical * se.ICV)
  CI.ACV <- cbind(beta.ACV - z.critical * se.ACV, beta.ACV + z.critical * se.ACV)
  RE.WEE <- rep(1, p)
  RE.ICV <- se.WEE^2 / se.ICV^2
  RE.ACV <- se.WEE^2 / se.ACV^2
  test.value <- as.vector(nn * t(phi.est) %*% VU.inv %*% phi.est)
  test.pvalue <- 1 - pchisq(test.value, q.num)
  compatibility.test <- c(value = test.value, df = q.num, pvalue = test.pvalue)
  format.CI <- function(CI) sprintf("[%.3f, %.3f]", CI[,1], CI[,2])
  result <- data.frame(EST_WEE = round(beta.WEE, 3), SE_WEE = round(se.WEE, 3), CI_WEE = format.CI(CI.WEE), RE_WEE = round(RE.WEE, 3), 
                       EST_ICV = round(beta.ICV, 3), SE_ICV = round(se.ICV, 3), CI_ICV = format.CI(CI.ICV), RE_ICV = round(RE.ICV, 3), 
                       EST_ACV = round(beta.ACV, 3), SE_ACV = round(se.ACV, 3), CI_ACV = format.CI(CI.ACV), RE_ACV = round(RE.ACV, 3), 
                       row.names = c("Age", "ER Status(+)"))
  list(result = result, tau.ACV = fit.ACV$tau,
       cure.rate = cure.rate,cure.test = curetest,
       compatibility.test = compatibility.test)
}

set.seed(1)
realdata.result <- getest(sdata.TCGA = sdata.TCGA, tpro = tpro, B = 400,
                          tstar = tstar, k.num = k.num, q.num = q.num)
sink("1.txt")
print(realdata.result)
sink()

# KMplot ------------------------------------------------------------------


survfit_to_df <- function(survfit_obj, group_name = "group") {
  df <- data.frame(
    time = survfit_obj$time,
    surv = survfit_obj$surv,
    lower = survfit_obj$lower,
    upper = survfit_obj$upper,
    group = rep(names(survfit_obj$strata), survfit_obj$strata)
  )
  return(df)
}

fit_overall <- survfit(Surv(yobs, delta) ~ 1, data = sdata.TCGA)

surv_df <- data.frame(
  time = fit_overall$time,
  surv = fit_overall$surv,
  lower = fit_overall$lower,
  upper = fit_overall$upper
)

p1 <- ggplot(surv_df, aes(x = time, y = surv)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "gray80",
    alpha = 0.3
  ) +
  geom_step(linewidth = 0.5, color = "black") +
  coord_cartesian(xlim = c(0, 25), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 24, by = 2)) +
  labs(
    x = "Survival Time (in Years)",
    y = "Survival Probability"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13)
  )

ggsave(
  filename = "kmplot.png",
  plot = p1,
  width = 8,
  height = 4.5,
  units = "in",
  dpi = 500
)
