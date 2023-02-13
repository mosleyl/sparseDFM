#' Fit a VAR(p) model 
#' 
#' @param x numeric data matrix 
#' @param p number of lags 
#' @noRd


VAR <- function(x, p) {
  n <- nrow(x)
  Y <- x[(p + 1):n, ]
  X <- c()
  for (i in 1:p) {
    X <- cbind(X, x[(p + 1 - i):(n - i), ])
  }
  A <- solve(t(X) %*% X) %*% t(X) %*% Y
  A[A >= 1] = 0.99
  res <- Y - X %*% A
  
  return(list(Y = Y, X = X, A = A, res = res))
}
