## ----setup, include=FALSE-----------------------------------------------------
## Global options
knitr::opts_chunk$set(cache = TRUE)

## ----message = FALSE----------------------------------------------------------
library(tidyverse)
library(boot)
library(bootstrap)
library(DAAG)
library(microbenchmark)
library(Rcpp)
library(coda)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("rmdformats")

## -----------------------------------------------------------------------------
grade = data.frame(name = c("小明","小红","小方"), chinese = c(90,89,93), math = c(88,95,89), english = c(94,95,91))

## -----------------------------------------------------------------------------
knitr::kable(grade)

## -----------------------------------------------------------------------------
grade$total = grade$chinese + grade$math + grade$english
knitr::kable(grade)

## -----------------------------------------------------------------------------
mpg %>% 
  ggplot(aes(x=class, y=hwy, fill=class)) + 
  geom_boxplot()

## -----------------------------------------------------------------------------
# implement of sample() when replace = TRUE
sample.me <- function(x, size, prob = NULL){
  if(!is.null(prob)){
    cump = cumsum(prob)
  }else{
      prob = rep(1/len(x),len(x)) # assign default value
      cump = cumsum(prob)
    }
  U = runif(size)
  return(x[findInterval(U,cump)+1]) # discrete inverse transformation
}

## -----------------------------------------------------------------------------
x = c(1,2,3,4)
prob = c(0.2,0.2,0.3,0.3)
samples = sample.me(x = x, size = 10000, prob = prob)
table(samples)

## -----------------------------------------------------------------------------
# The function to generate samples from standard Laplace distribution
rlaplace <- function(n){
  U = runif(n)
  samples = ifelse(U < 1/2, log(2*U),-log(2-2*U)) # continuous inverse transformation
  return(samples)
}

## -----------------------------------------------------------------------------
# generate samples and qqplot
n = 1000
samples = rlaplace(1000)
sampleq = seq(from = 1/(2*n),to = 1- 1/(2*n),by = 1/n)
samplet = sort(ifelse(samples < 0, 1/2*exp(samples), 1-1/2*exp(-samples)))
dat = data.frame(`q.therotical` = samplet, `q.sample` = sampleq)
ggplot(dat, aes(x = `q.therotical`, y = `q.sample`)) +
  geom_point() +
  geom_abline(slope = 1,intercept = 0) +
  ggtitle('qqplot of generated samples') +
  theme(plot.title = element_text(hjust=0.5,size=rel(1.5)))

## -----------------------------------------------------------------------------
# return the pdf of Beta(a,b) distribution at x
pdfBeta <- function(a,b,x){
  density = 1/beta(a,b)*x^(a-1)*(1-x)^(b-1)
  return(density)
}

# Use acceptance-rejection method to generate samples from Beta(a,b)
myBeta <- function(a,b,n){
  if((a <= 1)|(b <= 1)){
    stop('Input a,b must > 1.')
  }else{
    samples = NULL
    len = 0
    maxnum = pdfBeta(a,b,(a-1)/(a+b-2))
    # loop until there are enough samples
    while(len < n){
      U = runif(n)
      X = runif(n)
      sample.temp = X[U<=pdfBeta(a,b,X)/maxnum]
      samples = c(samples, sample.temp)
      len = len + length(sample.temp)
    }
    return(samples[1:n])
  }
}

## -----------------------------------------------------------------------------
#sample and draw
samples = myBeta(3,2,1000)
xticks = seq(from = 0, to = 1, length.out = 1000)
yticks = pdfBeta(3,2,xticks)
dat = data.frame(samples, xticks, yticks)
ggplot(dat)+
  geom_histogram(mapping = aes(x=samples,y=after_stat(density)),binwidth = 0.02,fill='grey')+
  geom_line(mapping = aes(x=xticks,y=yticks),color='red')

## -----------------------------------------------------------------------------
# write a function to sample from Rescaled Epanechnikov kernel
repane <- function(n){
  u1 = runif(n,-1,1)
  u2 = runif(n,-1,1)
  u3 = runif(n,-1,1)
  logic1 = (abs(u3)>=abs(u2))&(abs(u3)>=abs(u1))
  logic2 = !logic1
  samples = logic1*u2 + logic2*u3
  return(samples)
}

## -----------------------------------------------------------------------------
# generate samples
samples = repane(10000)
xticks = seq(from = -1, to = 1, length.out = 1000)
yticks = 3/4*(1-xticks^2)
dat = data.frame(samples, xticks, yticks)
ggplot(dat)+
  geom_histogram(mapping = aes(x=samples,y=after_stat(density)),binwidth = 0.02,fill='grey')+
  geom_line(mapping = aes(x=xticks,y=yticks),color='red')

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
mc1 <- mcmc(data=chains[1,])
mc2 <- mcmc(data=chains[2,])
mc3 <- mcmc(data=chains[3,])
mc4 <- mcmc(data=chains[4,])
mcs <- mcmc.list(list(mc1,mc2,mc3,mc4))
gelman.plot(mcs)

## -----------------------------------------------------------------------------
# MLE
u <- c(11,8,27,13,16,0,23,10,24,2)
v <- c(12,9,28,14,17,1,24,11,25,3)
Likelihood <- function(lambda){
  top <- -u * exp(-lambda * u) + v * exp(-lambda * v)
  down <- exp(-lambda * u) - exp(-lambda * v)
  return(sum(top/down))
}
lambda.star1 <- uniroot(Likelihood, interval = c(0,10))$root

lambda.star1

## -----------------------------------------------------------------------------
n <- length(u)
EM <- function(lambda0){
  up <- u * exp(-lambda0 * u)-v * exp(-lambda0 * v)
  down <- exp(-lambda0 * u) - exp(-lambda0 * v)
  un <- n/lambda0 + sum(up/down)
  return(n/un)
}
lambda <- 1
while (TRUE) {
  lambda.star2 <- EM(lambda)
  if(abs(lambda.star2-lambda)<0.0001){
    break
  }
  lambda <- lambda.star2
}

lambda.star2

## -----------------------------------------------------------------------------
solve.game <- function(A) {
m <- nrow(A)
n <- ncol(A)
it <- n^3
a <- c(rep(0, m), 1) 
A1 <- -cbind(t(A), rep(-1, n)) 
b1 <- rep(0, n)
A3 <- t(as.matrix(c(rep(1, m), 0))) 
b3 <- 1
sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=TRUE, n.iter=it)
a <- c(rep(0, n), 1) 
A1 <- cbind(A, rep(-1, m)) 
b1 <- rep(0, m)
A3 <- t(as.matrix(c(rep(1, n), 0))) 
b3 <- 1
sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=FALSE, n.iter=it)
soln <- list("A" = A,
"x" = sx$soln[1:m],
"y" = sy$soln[1:n],
"v" = sx$soln[m+1])
soln
}

## -----------------------------------------------------------------------------
# for game A
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)
s <- solve.game(A)
# solution
cat(round(cbind(s$x, s$y), 7))
# value
cat(s$v)

## -----------------------------------------------------------------------------
s <- solve.game(A+2)
# solution
round(cbind(s$x, s$y), 7)
# value
cat(s$v)

## -----------------------------------------------------------------------------
x <- c(1,2,3)
dim(x)

## -----------------------------------------------------------------------------
m <- matrix(1:4, ncol = 2)
is.array(m)

## -----------------------------------------------------------------------------
x <- data.frame(a=c(1,2),b=c("a","b"))
as.matrix(x)

## -----------------------------------------------------------------------------
A = data.frame(matrix(0, ncol=0, nrow=2))
A

## -----------------------------------------------------------------------------
B = data.frame(matrix(0, ncol=2, nrow=0))
A

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
X <- data.frame(a = c(1,2,3,4), b = c(2,1,4,5))
lapply(X, scale01)

## -----------------------------------------------------------------------------
Y <- data.frame(a = c(1,2,3,4), b = c("a","b","c","d"))
lapply(Y[sapply(Y, is.numeric)], scale01)

## -----------------------------------------------------------------------------
# a)
vapply(X, sd, numeric(1))

## -----------------------------------------------------------------------------
# b)
vapply(Y[vapply(Y, is.numeric, logical(1))], sd, numeric(1))

## -----------------------------------------------------------------------------
dir_cpp <- './Rcpp/'

source(paste0(dir_cpp,'GibbsR.R'))
sourceCpp(paste0(dir_cpp,'GibbsC.cpp'))
ts <- microbenchmark(gibbR=gibbsR(10, 1, 3, 100,10), gibbC=gibbsC(10, 1, 3, 100,10))
summary(ts)[,c(1,3,5,6)]

