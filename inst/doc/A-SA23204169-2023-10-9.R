## ----setup, include=FALSE-----------------------------------------------------
## Global options
knitr::opts_chunk$set(cache = TRUE)

## -----------------------------------------------------------------------------
set.seed(0)
g <- function(x) 
  return(x^2*exp(-0.5*x^2)/sqrt(2*pi)*(x>1))
f <- function(x) 
  return(1/sqrt(4*pi)*exp(-0.25*x^2))
#C = 1 / integrate(f, 1, Inf)$value
#e <- function(x)
  #return(g(x) / (C * f(x)))
x = rnorm(1e5, mean = 0, sd=sqrt(2))
est = mean(g(x)/f(x))
v = var((g(x)/f(x))[x>1])
cat('Mean: ', est,"; Variance: ", v)

## -----------------------------------------------------------------------------
gdivf <- function(x){
  g <- exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
  f <- 5 * exp(-x)/ (1-exp(-1))
  return(g/f)
}
gen_data <- function(n,left,right){
  dat = NULL
  while(length(dat)<n){
    u = runif(n)
    x <- - log(1 - u * (1 - exp(-1)))
    ind <- x>left & x<right
    x <- x[ind]
    dat <- c(dat, x[1:min(length(x), n-length(dat))])
  }
  return(dat)
}
set.seed(0)
theta = NULL
for (i in 1:5) {
  x = gen_data(2000,(i-1)/5,i/5)
  temp = gdivf(x)
  theta = rbind(theta, temp)
}
theta0 = colSums(theta)
cat("Mean:", mean(theta0),";Se:", sd(theta0), "\n")
cat("In example 5.10, Mean:", 0.52506988, ";Se:", 0.09658794)

## -----------------------------------------------------------------------------
set.seed(0)
n = 20
result = rep(FALSE, 10000)
alpha = qt(p = 0.975,df=19)
for (i in 1:10000) {
  x = rchisq(n, df=2)
  mu = mean(x)
  se = sd(x)
  left = mu - alpha * se/sqrt(20)
  right = mu + alpha * se/sqrt(20)
  if(left < 2 && right > 2){
    result[i] = TRUE
  }
}
cat("The probability that t-interval can cover real value is",mean(result), ".")

## -----------------------------------------------------------------------------
# simulation example 6.4
set.seed(0)
n = 20
result = rep(FALSE, 10000)
alpha = qchisq(0.05, df=n-1)
for (i in 1:10000) {
  x <- rnorm(n, mean=0, sd=2)
  UCL <- (n-1) * var(x) / alpha
  if(UCL > 4){
    result[i] = TRUE
  }
}
cat("The probability that 95% confidence interval can cover real value is",mean(result), ".")

## -----------------------------------------------------------------------------
# chi-square distribution
set.seed(0)
n = 20
result = rep(FALSE, 10000)
alpha = qt(p = 0.975,df=19)
for (i in 1:10000) {
  x = rchisq(n, df=1)
  mu = mean(x)
  se = sd(x)
  left = mu - alpha * se/sqrt(20)
  right = mu + alpha * se/sqrt(20)
  if(left < 1 && right > 1){
    result[i] = TRUE
  }
}
cat("The probability that t-interval can cover real value in chi-square distribution is",mean(result), ".")

## -----------------------------------------------------------------------------
# uniform distribution
set.seed(0)
n = 20
result = rep(FALSE, 10000)
alpha = qt(p = 0.975,df=19)
for (i in 1:10000) {
  x = runif(n, 0, 2)
  mu = mean(x)
  se = sd(x)
  left = mu - alpha * se/sqrt(20)
  right = mu + alpha * se/sqrt(20)
  if(left < 1 && right > 1){
    result[i] = TRUE
  }
}
cat("The probability that t-interval can cover real value in uniform distribution is",mean(result), ".")

## -----------------------------------------------------------------------------
# Exponential distribution
set.seed(0)
n = 20
result = rep(FALSE, 10000)
alpha = qt(p = 0.975,df=19)
for (i in 1:10000) {
  x = rexp(n, rate=1)
  mu = mean(x)
  se = sd(x)
  left = mu - alpha * se/sqrt(20)
  right = mu + alpha * se/sqrt(20)
  if(left < 1 && right > 1){
    result[i] = TRUE
  }
}
cat("The probability that t-interval can cover real value in exponential distribution is",mean(result), ".")

