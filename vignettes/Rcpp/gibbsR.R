gibbsR <- function(n, a, b, N, thin) {
  mat <- matrix(nrow = N, ncol = 2)
  x <- y <- 0
  for (i in 1:N) {
    for (j in 1:thin) {
      x <- rbinom(1, n, y)
      y <- rbeta(1, x + a, n - x + b)
    }
    mat[i, ] <- c(x, y)
  }
  mat
}
