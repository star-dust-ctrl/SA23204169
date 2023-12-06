## ----setup, include=FALSE-----------------------------------------------------
## Global options
knitr::opts_chunk$set(cache = TRUE)

## -----------------------------------------------------------------------------
library(boot)
library(bootstrap)
library(DAAG)

## -----------------------------------------------------------------------------
b.mean <- function(x,i) mean(x[i])
obj <- boot(aircondit$hours, statistic = b.mean, R=1e4)
theta.hat <- mean(aircondit$hours)
se.hat <- sd(obj$t[,1])
cat("Standard normal CI:[", theta.hat-qnorm(0.975)*se.hat, ",", theta.hat+qnorm(1-0.05/2)*se.hat, "]\n", sep = "")
cat("Basic CI:[", 2*theta.hat-quantile(obj$t[,1], 0.975), ",", 2*theta.hat-quantile(obj$t[,1], 0.025), "]\n", sep = "")
cat("Percentile CI:[", quantile(obj$t[,1], 0.025), ",", quantile(obj$t[,1], 0.975), "]\n", sep = "")
jack <- jackknife(aircondit$hours, mean)
L <- mean(jack$jack.values) - jack$jack.values
a.hat <- sum(L^3)/(6 * sum(L^2)^1.5)
z0.hat <- qnorm(mean(obj$t <theta.hat))
alpha1 = pnorm(z0.hat + (z0.hat + qnorm(0.025))/(1-a.hat*(z0.hat + qnorm(0.025))))
alpha2 = pnorm(z0.hat + (z0.hat + qnorm(0.975))/(1-a.hat*(z0.hat + qnorm(0.975))))
cat("BCa CI:[", quantile(obj$t[,1], alpha1), ",", quantile(obj$t[,1], alpha2), "]\n", sep = "")

## -----------------------------------------------------------------------------
eig <- eigen(cov(scor))$values
theta.hat <- eig[1]/sum(eig)
n = dim(scor)[1]
theta.jack <- numeric(n)
for (i in 1:n) {
  dat = scor[-i,]
  eig.temp <- eigen(cov(dat))$values
  theta.jack[i] <- eig.temp[1]/sum(eig.temp)
}
bias = (n-1)*(mean(theta.jack)-theta.hat)
se = (n-1)^2 / n * var(theta.jack)
cat("Bias estimation:", bias, "\n")
cat("Standard error estimation:", se, "\n")

## -----------------------------------------------------------------------------
n = dim(ironslag)[1]
e1 <- e2 <- e3 <- e4 <- numeric(n*(n-1))
k <- 1
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    exclude <- c(i,j)
    dat <- ironslag[-exclude,]
    x <- dat$chemical
    y <- dat$magnetic
    
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * ironslag$chemical[exclude]
    e1[(2*k-1):(2*k)] <- ironslag$magnetic[exclude] - yhat1
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * ironslag$chemical[exclude] +
    J2$coef[3] * ironslag$chemical[exclude]^2
    e2[(2*k-1):(2*k)] <- ironslag$magnetic[exclude] - yhat2
    
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * ironslag$chemical[exclude]
    yhat3 <- exp(logyhat3)
    e3[(2*k-1):(2*k)] <- ironslag$magnetic[exclude] - yhat3

    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(ironslag$chemical[exclude])
    yhat4 <- exp(logyhat4)
    e4[(2*k-1):(2*k)] <- ironslag$magnetic[exclude] - yhat4
    
    k = k + 1
  }
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

