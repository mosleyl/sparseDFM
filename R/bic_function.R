#' BIC calculation for LASSO tuning 
#'
#' @param X: n x p data matrix 
#' @param factors: n x r factor estimates 
#' @param loadings: p x r loading estimates 
#' 
#' @importFrom stats na.omit
#' 
#' @noRd



bic_function <- function(X, factors, loadings){
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  W = !is.na(X)
  
  loglik = 0 
  
  for(t in 1:n){
    
    new_Lambda <- loadings[W[t,],]
    error <- as.numeric(na.omit(X[t,])) - new_Lambda %*% factors[t,]
    forb <- rowSums(t(error) %*% error)
    loglik = loglik + forb 
    
  }
  
  bic = log(loglik/(n*p)) + sum(as.vector(loadings) != 0)*log(n*p)/(n*p)
  
  return(bic)
  
}




