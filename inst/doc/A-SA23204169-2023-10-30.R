## ----setup, include=FALSE-----------------------------------------------------
## Global options
knitr::opts_chunk$set(cache = TRUE)

## -----------------------------------------------------------------------------
# ecdf
F.n <- function(para, vec){
    result <- numeric(length(vec))
    for (i in 1:length(vec)) {
      result[i] <- mean(vec[i]>=para)
    }
    return(result)
  }
# completion of Two-Sample Cramer-von Mises test
cmv.test <- function(x,y){
  n = length(x)
  m = length(y)
  res = m * n / (m+n)^2 * (sum((F.n(x,x)-F.n(y,x))^2)+sum((F.n(x,y)-F.n(y,y))^2))
  return(res)
}

## -----------------------------------------------------------------------------
# permutation
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)
D0 <- cmv.test(x,y)
R <- 999
z <- c(x, y)
K <- 1:26
D <- numeric(R)
set.seed(0)
for (i in 1:R) {
  k <- sample(K, size = 14, replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k]
  D[i] <- cmv.test(x1, y1)
}
p <- mean(c(D0, D) >= D0)
p

## -----------------------------------------------------------------------------
count5 <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx, outy))))
}

counttest <- function(x,y, repsize = 1000, seed = 0){
  n1 <- length(x)
  n2 <- length(y)
  z <- c(x, y)
  N <- n1 + n2
  K <- 1:N
  D0 <- count5(x, y)
  D <- numeric(repsize - 1)
  set.seed(seed)
  for (i in 1:(repsize - 1)) {
    k <- sample(K,size = n1, replace = FALSE)
    x1 <- z[k]
    y1 <- z[-k]
    D[i] <- count5(x1, y1)
  }
  p <- mean(c(D0, D) >= D0)
  return(p)
}

## -----------------------------------------------------------------------------
# simulate when sigma1 = sigma2
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- 1
sigma2 <- 1
m <- 1000
set.seed(0)
pvalue <- replicate(m, expr={
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
counttest(x, y)
})
table(pvalue)

## -----------------------------------------------------------------------------
# simulate when sigma1 != sigma2
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- 1
sigma2 <- 2
m <- 1000
set.seed(0)
pvalue <- replicate(m, expr={
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
counttest(x, y)
})
table(pvalue)

