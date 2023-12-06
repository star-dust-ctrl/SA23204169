## ----setup, include=FALSE-----------------------------------------------------
## Global options
knitr::opts_chunk$set(cache = TRUE)

## -----------------------------------------------------------------------------
# load plot library
library(ggplot2)

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

