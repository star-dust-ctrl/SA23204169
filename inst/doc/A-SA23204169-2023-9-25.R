## ----setup, include=FALSE-----------------------------------------------------
## Global options
knitr::opts_chunk$set(cache = TRUE)

## -----------------------------------------------------------------------------
# set random seed
set.seed(0)

# settings
N = 1e6
d = 1
rho = c(0.3,0.6,1)
for (l in rho) {
  X <- matrix(runif(N*100, 0, d/2),nrow = 100)
  Y <- matrix(runif(N*100, 0, pi/2),nrow = 100)
  pi.hat = 2*l/d/rowMeans(l/2*sin(Y)>X)
  cat("The variance of estimator when rho =",l,"is",var(pi.hat),"\n")
  rm(list = c("X","Y"))
  gc()
}

## -----------------------------------------------------------------------------
# settings
N = 1e5
cyc = 100
set.seed(0)

U1 = matrix(runif(N*cyc), nrow = cyc)
U2.pre = matrix(runif(N/2*cyc),nrow = cyc)
U2 = cbind(U2.pre, 1-U2.pre)
theta1=rowMeans(exp(U1))
theta2=rowMeans(exp(U2))

cat("The output of Simple MC:", mean(theta1), "Variance:", var(theta1),"\n",
    "The output of Antithetic variables:", mean(theta2), "Variance:", var(theta2),"\n",
    "The variance is decreased to", var(theta2)/var(theta1), "namely by", 1-var(theta2)/var(theta1))

