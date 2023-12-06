//' @import Rcpp
//' @useDynLib SA23204169

#include <Rcpp.h>
#include <random>
using namespace Rcpp;
//' @title K-Means clustering using Rcpp
//' @description Perform k-means clustering on a data matrix using Rcpp.
//' @param data numeric matrix of data, each row representing a sample
//' @param k the number of clusters
//' @param maxIter the maximum number of iterations allowed
//' @return an integer cluster vector, with each integer representing a group.
//' @examples
//' dat <- as.matrix(iris[,1:2])
//' kmeansCpp(dat, k=3)
//' @export
// [[Rcpp::export]]
IntegerVector kmeansCpp(NumericMatrix data, int k, int maxIter = 100) {
  int n = data.nrow();
  int p = data.ncol();
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, n - 1);
  

  NumericMatrix centroids(k, p);
  for (int i = 0; i < k; ++i) {
    int idx = dis(gen) % n;
    for (int j = 0; j < p; ++j) {
      centroids(i, j) = data(idx, j);
    }
  }
  
  IntegerVector clusterAssign(n);
  NumericVector distances(n);
  
  for (int iter = 0; iter < maxIter; ++iter) {
    for (int i = 0; i < n; ++i) {
      double minDist = std::numeric_limits<double>::max();
      int closestCluster = -1;
      
      for (int j = 0; j < k; ++j) {
        double dist = 0.0;
        for (int l = 0; l < p; ++l) {
          dist += pow(data(i, l) - centroids(j, l), 2);
        }
        
        if (dist < minDist) {
          minDist = dist;
          closestCluster = j;
        }
      }
      
      clusterAssign[i] = closestCluster;
      distances[i] = minDist;
    }

    for (int i = 0; i < k; ++i) {
      for (int j = 0; j < p; ++j) {
        double sum = 0.0;
        int count = 0;
        for (int l = 0; l < n; ++l) {
          if (clusterAssign[l] == i) {
            sum += data(l, j);
            count++;
          }
        }
        if (count > 0) {
          centroids(i, j) = sum / count;
        }
      }
    }
  }
  
  IntegerVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = clusterAssign[i] + 1;
  }
  return result;
}