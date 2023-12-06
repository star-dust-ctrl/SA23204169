#' @title K-Means clustering using R
#' @description Perform k-means clustering on a data matrix using R.
#' @param data numeric matrix of data, each row representing a sample
#' @param k the number of clusters
#' @param maxIter the maximum number of iterations allowed
#' @return an integer cluster vector, with each integer representing a group.
#' @examples
#' dat <- as.matrix(iris[,1:2])
#' kmeansR(dat, k=3)
#' @export
kmeansR <- function(data, k, maxIter = 100) {
  n <- nrow(data)
  p <- ncol(data)
  
  centroids <- matrix(0, nrow = k, ncol = p)
  for (i in 1:k) {
    idx <- sample(1:n, size = 1)
    centroids[i, ] <- data[idx, ]
  }
  
  clusterAssign <- numeric(n)
  distances <- numeric(n)
  
  for (iter in 1:maxIter) {
    for (i in 1:n) {
      minDist <- Inf
      closestCluster <- -1
      
      for (j in 1:k) {
        dist <- sum((data[i,] - centroids[j,])^2)
        
        if (dist < minDist) {
          minDist <- dist
          closestCluster <- j
        }
      }
      
      clusterAssign[i] <- closestCluster
      distances[i] <- minDist
    }
    
    for (i in 1:k) {
      count <- sum(clusterAssign == i)
      if (count > 0) {
        centroids[i, ] <- colMeans(data[clusterAssign == i, ])
      } else {
        idx <- sample(1:n, size = 1)
        centroids[i, ] <- data[idx, ]
      }
    }
  }
  
  return(clusterAssign)
}