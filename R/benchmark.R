#' @import microbenchmark
#' @title Benchmark R and Rcpp functions
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of R function \code{kmeansR} and Cpp function \code{kmeansCpp}.
#' @examples
#' \dontrun{
#' dat <- as.matrix(iris[,1:2])
#' ts <- microbenchmark::microbenchmark(
#'   kmeansR = kmeansR(dat,3),
#'   kmeansCpp = kmeansCpp(dat,3),
#'   kmeans = kmeans(dat, 3)
#' )
#' print(summary(ts)[,c(1,3,5,6)])
#'#       expr       lq   median       uq
#'#1   kmeansR 49994.55 50236.30 52378.50
#'#2 kmeansCpp   128.35   137.85   143.85
#'#3    kmeans   113.85   141.65   249.65
#'}
#' @details Using Rcpp while implementing K-means improves the efficiency, however the built-in function of R is still faster. We also
#' note that I use lots of loops in the R implementation, so there exactly may be better solution, however, the speed-up effect of Rcpp is obvious though.
NULL