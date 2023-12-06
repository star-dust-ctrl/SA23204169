## ----setup, include=FALSE-----------------------------------------------------
## Global options
knitr::opts_chunk$set(cache = TRUE)

## -----------------------------------------------------------------------------
M = 1000
m = 1000
m1 = m * 0.95
m2 = m * 0.05
alpha = 0.1
bon.FWER <- 0
bon.V  <- 0
bon.S <- 0
bon.T <- 0
BH.FWER <- 0
BH.V  <- 0
BH.S <- 0
BH.T <- 0
set.seed(0)
for (i in 1:M) {
  p1 <- runif(m1)
  p2 <- rbeta(m2, 0.1, 1)
  p <- c(p1, p2)
  p.bonferroni <- p.adjust(p, method = "bonferroni")
  p.BH <- p.adjust(p, method = "BH")
  bon.FWER <- bon.FWER + (any(p.bonferroni[1:m1] <= alpha))/M
  bon.V <- bon.V + sum((p.bonferroni[1:m1] <= alpha))
  bon.S <- bon.S + sum((p.bonferroni[(m1+1):m] <= alpha))
  bon.T <- bon.T + sum((p.bonferroni[(m1+1):m] > alpha))
  BH.FWER <- BH.FWER + (any(p.BH[1:m1] <= alpha))/M
  BH.V <- BH.V + sum((p.BH[1:m1] <= alpha))
  BH.S <- BH.S + sum((p.BH[(m1+1):m] <= alpha))
  BH.T <- BH.T + sum((p.BH[(m1+1):m] > alpha))
}
bon.FDR <- bon.V/(bon.V + bon.S) 
bon.TPR <- bon.S/(bon.T + bon.S)
BH.FDR <- BH.V/(BH.V + BH.S) 
BH.TPR <- BH.S/(BH.T + BH.S)
result <- data.frame(Method = c("Bonferroni","BH"), FWER = c(bon.FWER,BH.FWER), FDR = c(bon.FDR,BH.FDR), TPR = c(bon.TPR, BH.TPR))
knitr::kable(result)

## -----------------------------------------------------------------------------
lambda = 2
B = 1000
m = 1000

## ---- fig.height=5, fig.width=5-----------------------------------------------
# n = 5
n = 5
set.seed(0)
bias.r <- numeric(m)
sd.r <- numeric(m)
for (i in 1:m) {
  x <- rexp(n, lambda)
  temp <- numeric(B)
  for (j in 1:B) {
    x.sample <- sample(x, size = n, replace = TRUE)
    temp[j] <- 1/mean(x.sample)
  }
  bias.r[i] <- mean(temp) - lambda
  sd.r[i] <- sd(temp)
}
par(mfrow = c(2,1))
hist(bias.r,breaks = 30);abline(v=lambda/(n-1),col='red',lwd=2)
hist(sd.r, breaks = 30);abline(v=lambda * n/((n-1)*sqrt(n-2)),col='red',lwd=2)

## ---- fig.height=5, fig.width=5-----------------------------------------------
# n = 10
n = 10
set.seed(0)
bias.r <- numeric(m)
sd.r <- numeric(m)
for (i in 1:m) {
  x <- rexp(n, lambda)
  temp <- numeric(B)
  for (j in 1:B) {
    x.sample <- sample(x, size = n, replace = TRUE)
    temp[j] <- 1/mean(x.sample)
  }
  bias.r[i] <- mean(temp) - lambda
  sd.r[i] <- sd(temp)
}
par(mfrow = c(2,1))
hist(bias.r,breaks = 30);abline(v=lambda/(n-1),col='red',lwd=2)
hist(sd.r, breaks = 30);abline(v=lambda * n/((n-1)*sqrt(n-2)),col='red',lwd=2)

## ---- fig.height=5, fig.width=5-----------------------------------------------
# n = 20
n = 20
set.seed(0)
bias.r <- numeric(m)
sd.r <- numeric(m)
for (i in 1:m) {
  x <- rexp(n, lambda)
  temp <- numeric(B)
  for (j in 1:B) {
    x.sample <- sample(x, size = n, replace = TRUE)
    temp[j] <- 1/mean(x.sample)
  }
  bias.r[i] <- mean(temp) - lambda
  sd.r[i] <- sd(temp)
}
par(mfrow = c(2,1))
hist(bias.r,breaks = 30);abline(v=lambda/(n-1),col='red',lwd=2)
hist(sd.r, breaks = 30);abline(v=lambda * n/((n-1)*sqrt(n-2)),col='red',lwd=2)

## -----------------------------------------------------------------------------
library(bootstrap)
B <- 200
n <- nrow(law)
R <- numeric(B)
set.seed(0)
for (b in 1:B) {
  i <- sample(1:n, size = n, replace = TRUE)
  LSAT <- law$LSAT[i]
  GPA <- law$GPA[i]
  R[b] <- cor(LSAT, GPA)
}
mean.R <- mean(R)
se.R <- sd(R)
t <- qt(p = 0.975, df = n-1)
cat("The bootstrap t confidence interval for correlation statistic is [", 
    mean.R - t * se.R, ",", mean.R + t * se.R, "]")

