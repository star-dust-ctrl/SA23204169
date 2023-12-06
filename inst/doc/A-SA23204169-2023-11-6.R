## ----setup, include=FALSE-----------------------------------------------------
## Global options
knitr::opts_chunk$set(cache = TRUE)

## -----------------------------------------------------------------------------
# the solve function
eva.a <- function(N, b1, b2, b3, f0){
  set.seed(0)
  x1 <- rpois(N, 1)
  x2 <- rexp(N, 1)
  x3 <- sample(c(1,0),size = N, replace = TRUE)
  f <- function(a, x1, x2, x3, f0){
    temp = exp(a+b1*x1+b2*x2+b3*x3)
    return(mean(1/(1+temp))-f0)
  }
  return(uniroot(f, interval = c(-10,10), x1, x2, x3, f0)$root)
}

## -----------------------------------------------------------------------------
# call the function
N = 1e6; b1 = 0; b2=1; b3=-1; f=c(0.1,0.01,0.001,0.0001);
a = numeric(4)
for (i in 1:4) {
  a[i] <- eva.a(N, b1, b2, b3, f[i])
}

## ---- fig.height=5, fig.width=5-----------------------------------------------
### draw the picture
neglogf <- -log(f)
plot(neglogf, a, type="o")

## -----------------------------------------------------------------------------
dlaplace <- function(x){
  return(1/2*exp(-abs(x)))
}

# Laplace chain
laplace.Metropolis <- function(sigma, x0, N){
  set.seed(0)
  x <- numeric(N)
  x[1] <- x0
  U <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (U[i] <= dlaplace(y)/dlaplace(x[i-1])){
      x[i] <- y
    }else{
      x[i] <- x[i-1]
      k <- k + 1
    }
  }
  return(list(x=x, k=k))
}

## -----------------------------------------------------------------------------
# Generate chain
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
mc1 <- laplace.Metropolis(sigma[1], x0, N)
mc2 <- laplace.Metropolis(sigma[2], x0, N)
mc3 <- laplace.Metropolis(sigma[3], x0, N)
mc4 <- laplace.Metropolis(sigma[4], x0, N)

## ---- fig.height=5, fig.width=5-----------------------------------------------
# Plot chain
par(mfrow = c(2,2))
mc <- cbind(mc1$x, mc2$x, mc3$x, mc4$x)
for (j in 1:4) {
  plot(mc[,j], type = "l",
       xlab = bquote(sigma == .(round(sigma[j],3))),
       ylab = "X", ylim = range(mc[,j]))
}
par(mfrow = c(1,1))

## -----------------------------------------------------------------------------
# print acceptance rate
print(1-c(mc1$k, mc2$k, mc3$k, mc4$k)/N)

## ---- fig.height=5, fig.width=5-----------------------------------------------
#initialize
N <- 5000 
burn <- 1000 
X <- matrix(0, N, 2) #the chain, a bivariate sample
rho <- 0.9 #correlation
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

X[1, ] <- c(mu1, mu2)
for (i in 2:N) {
  x2 <- X[i-1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
  X[i, 2] <- rnorm(1, m2, s2)
}
b <- burn + 1
x <- X[b:N, ]
plot(x, main="", cex=.5, xlab=bquote(X[1]),
ylab=bquote(X[2]), ylim=range(x[,2]))

## ---- fig.height=5, fig.width=5-----------------------------------------------
lr = lm(x[,2]~x[,1])
barplot(lr$residuals, space = 0)

## ---- fig.height=5, fig.width=5-----------------------------------------------
# check with manual code
Gelman.Rubin <- function(psi) {
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)
  B <- n * var(psi.means) 
  psi.w <- apply(psi, 1, "var") 
  W <- mean(psi.w) 
  v.hat <- W*(n-1)/n + (B/n) 
  r.hat <- v.hat / W 
  return(r.hat)
}

f <- function(x, sigma) {
  if (any(x < 0)) return (0)
  stopifnot(sigma > 0)
  return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}

gen.chain <- function(x0, N, sigma){
  set.seed(0)
  x = numeric(N)
  x[1] <- x0
  k <- 0
  u = runif(N)
  for (i in 2:N) {
    xt = x[i-1]
    y = rchisq(1, df=xt)
    if(u[i] <= (f(y,sigma)*dchisq(xt, df=y))/(f(xt,sigma)*dchisq(y,df=xt))){
      x[i] = y
      k <- k + 1
    }else{
      x[i] = xt
    }
  }
  return(list(x=x, k=k))
}


N <- 1e4
sigma <- 4
b <- 1e3
x0 <- c(1,2,4,8)
chains <- matrix(0, nrow = 4, ncol = N)
for (i in 1:4) {
  chains[i,] <- gen.chain(x0[i], N, sigma)$x
}
psi = t(apply(chains, 1, cumsum))
for (i in 1:nrow(psi)){
  psi[i,] <- psi[i,] / (1:ncol(psi))
}
rhat <- rep(0, N)
for (j in (b+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## ---- fig.height=5, fig.width=5-----------------------------------------------
# check with coda
library(coda)
mc1 <- mcmc(data=chains[1,])
mc2 <- mcmc(data=chains[2,])
mc3 <- mcmc(data=chains[3,])
mc4 <- mcmc(data=chains[4,])
mcs <- mcmc.list(list(mc1,mc2,mc3,mc4))
gelman.plot(mcs)

