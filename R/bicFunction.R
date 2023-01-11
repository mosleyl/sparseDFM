
bicFunction <- function(X, factors, loadings){
  
  ## Input:
  ##
  ##  X: n x p data matrix 
  ##  factors: n x r factor estimates 
  ##  loadings: p x r loading estimates 
  
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




