## ----setup, include=FALSE-----------------------------------------------------
## Global options
knitr::opts_chunk$set(cache = TRUE)

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
library(boot)
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

