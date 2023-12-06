## ----setup, include=FALSE-----------------------------------------------------
## Global options
knitr::opts_chunk$set(cache = TRUE)

## -----------------------------------------------------------------------------
library(SA23204169)
library(mvtnorm)
set.seed(1)
covmatrix <- matrix(0.3, ncol = 5, nrow = 5)
diag(covmatrix) <- 1
meanvalue <- rep(0,5)
X <- rmvnorm(100, meanvalue, covmatrix)
beta <- runif(5, min = 0.5, max = 1)
mu <- sample(c(-2, 2), size = 100, replace = TRUE)
eps <- rnorm(100, 0, 0.5)
y <- as.vector(X %*% beta) + mu + eps

## -----------------------------------------------------------------------------
# initialize the solver
solver <- SubgroupSolver(X,y)

## -----------------------------------------------------------------------------
solver$fit("MCP", 0.6, 3, 1)
solver$mu

## -----------------------------------------------------------------------------
library(microbenchmark)
set.seed(0)
dat <- as.matrix(iris[,1:2])
ts <- microbenchmark(
  kmeansR = kmeansR(dat,3),
  kmeansCpp = kmeansCpp(dat,3),
  kmeans = kmeans(dat, 3)
)
print(summary(ts)[,c(1,3,5,6)])

