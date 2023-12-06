#' @import MASS
#' @import Rcpp
#' @import methods
#' @import mvtnorm
library(MASS)
library(Rcpp)

pclip <- function(x, a, b){
  s = pmax(a, pmin(x, b))
  return(s)
}

ST <- function(t,lamb){
  s = sign(t) * pclip(abs(t)-lamb, 0, Inf)
  return(s)
}

MCP_eta <- function(lamb, gamma, theta, delta){
  index = (abs(delta) <= gamma * lamb)
  eta <- ifelse(index, ST(delta, lamb/theta)/(1-1/(gamma*theta)), delta)
  return(eta)
}

SCAD_eta <- function(lamb, gamma, theta, delta){
  eta <- ST(delta, gamma*lamb/((gamma-1)*theta))/(1-1/((gamma-1)*theta))
  idx1 <- abs(delta) <= (lamb + lamb/theta)
  idx2 <- abs(delta) > gamma*lamb
  eta[idx1] <- ST(delta, lamb/theta)[idx1]
  eta[idx2] <- delta[idx2]
  return(eta)
}

L1_eta <- function(lamb, theta, omega, delta){
  eta = ST(delta, lamb*omega/theta)
  return(eta)
}


#' @title A subgroup solver based on ADMM
#' @description Here I implement the method proposed by Shujie Ma & Jian Huang in 2016.
#' @references Shujie Ma & Jian Huang(2016).A concave pairwise fusion approach to subgroup analysis
#' @examples
#' set.seed(1)
#' covmatrix <- matrix(0.3, ncol = 5, nrow = 5)
#' diag(covmatrix) <- 1
#' meanvalue <- rep(0,5)
#' X <- mvtnorm::rmvnorm(100, meanvalue, covmatrix)
#' beta <- runif(5, min = 0.5, max = 1)
#' mu <- sample(c(-2, 2), size = 100, replace = TRUE)
#' eps <- rnorm(100, 0, 0.5)
#' y <- as.vector(X %*% beta) + mu + eps
#' solver <- SubgroupSolver(X,y)
#' solver$fit("MCP",0.6, 3, 1)
#' solver$mu
#' @exportClass SubgroupSolver
#' @export SubgroupSolver
SubgroupSolver <- setRefClass('SubgroupSolver', fields = c(X="ANY", y="ANY", Qx="ANY", DELTA="ANY", n="ANY", p="ANY", 
                                                           epoch="ANY", mu="ANY", beta="ANY", eta="ANY", v="ANY"),
                              methods = list(
                                initialize = function(X,y){
                                "
                                Initialize a SubgroupSolver.
                                
                                This function creates and initializes a new instance of the SubgroupSolver class. 
                                
                                Here X is a numeric matrix and y is a numeric vector, we apply \\code{solver = SubgroupSolver(X,y)} to initialize a object.
                                "
                                X <<- X
                                y <<- y
                                n <<- dim(X)[1]
                                p <<- dim(X)[2]
                                I <- diag(n)
                                delta <- matrix(rep(0,n^2*(n-1)/2),nrow=n)
                                k <- 1
                                for (i in 1:(n-1)) {
                                  for (j in (i+1):n) {
                                    delta[,k] <- I[,i] - I[, j]
                                    k <- k + 1
                                  }
                                }
                                Qx <<- X %*% ginv(t(X) %*% X) %*% t(X)
                                DELTA <<- t(delta)
                                epoch <<- 0
                                rg <- lm(y~X)
                                beta <<- unname(rg$coefficients)[2:(p+1)]
                                mu <<- as.vector(y - X %*% beta)
                                muij <- matrix(rep(mu,n), ncol = n)
                                etamatrix = muij - t(muij)
                                eta <<- t(etamatrix)[lower.tri(t(etamatrix))]
                                v <<- rep(0, n*(n-1)/2)
                              },
                              
                              update = function(method, lamb, gamma, theta, omega=NULL){
                                "Update the solver under some given parameters, usually useless for user.
                                "
                                epoch <<- epoch + 1
                                In = diag(n)
                                mu <<- as.vector(ginv(theta * t(DELTA) %*% DELTA + In - Qx) %*% ((In - Qx)%*%y + theta * t(DELTA) %*% (eta - v/theta)))
                                beta <<- as.vector(ginv(t(X) %*% X) %*% t(X) %*% (y - mu))
                                muij = matrix(rep(mu, n), nrow=n)
                                etamatrix = muij - t(muij)
                                temp <- t(etamatrix)[lower.tri(t(etamatrix))]
                                deltaij = t(etamatrix)[lower.tri(t(etamatrix))] + v / theta
                                if(method == 'MCP'){
                                  eta <<- MCP_eta(lamb, gamma, theta, deltaij)
                                }else{
                                  if(method == 'SCAD'){
                                    eta <<- SCAD_eta(lamb, gamma, theta, deltaij)
                                  }else{
                                    eta <<- L1_eta(lamb, theta, omega, deltaij)
                                  }
                                }
                                v <<- v + theta * (temp - eta)
                              },
                              
                              
                              fit = function(method, lamb, gamma, theta, eps=1e-3, omega = NULL){
                                "Giving the solver the method and some parameters and start to solve the subgroup problem.
                                
                                 The parameter \\code{method} must be one of ('MCP','SCAD' and 'L1'), which applies to different penalties in the paper.
                                 
                                 MCP:
                                 
                                 \\deqn{p_{\\gamma}(t,\\lambda)=\\lambda\\int_0^t(1-x/(\\gamma\\lambda))_+ dx}
                                 
                                 SCAD:
                                 
                                 \\deqn{p_{\\gamma}(t,\\lambda)=\\lambda\\int_0^t\\min\\{1,(\\gamma-x/\\lambda)_+ /(\\gamma-1)\\} dx}
                                 
                                 L1:
                                 
                                 \\deqn{p_{\\omega}(t,\\lambda)=\\lambda \\omega_{ij} |t|}
                                 
                                 and \\code{theta} is another parameter in the process of ADMM algorithm.
                                "
                                epoch <<- 0
                                rg <- lm(y~X)
                                beta <<- unname(rg$coefficients)[2:(p+1)]
                                mu <<- as.vector(y - X %*% beta)
                                muij <- matrix(rep(mu,n), ncol = n)
                                etamatrix = muij - t(muij)
                                eta <<- t(etamatrix)[lower.tri(t(etamatrix))]
                                v <<- rep(0, n*(n-1)/2)
                                if (!(method %in% c('MCP','SCAD','L1'))){
                                  stop("Parameter 'method' must be one of 'MCP', 'SCAD' and 'L1'.")
                                }else{
                                  if(method == 'L1' & is.null(omega)){
                                    stop("L1 must be applied with omega.")
                                  }
                                }
                                while(TRUE){
                                  update(method, lamb, gamma, theta, omega)
                                  r <- as.vector(DELTA %*% mu - eta)
                                  if(norm(r,"2") < eps){
                                    break
                                  }
                                }
                              }))
