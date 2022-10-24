VAR <- function(x, p) {
  n <- nrow(x)
  Y <- x[(p + 1):n, ]
  X <- c()
  for (i in 1:p) {
    X <- cbind(X, x[(p + 1 - i):(n - i), ])
  }
  A <- solve(t(X) %*% X) %*% t(X) %*% Y
  res <- Y - X %*% A
  
  return(list(Y = Y, X = X, A = A, res = res))
}
