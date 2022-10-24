
bic_function <- function(X, model_fit){
  
  ## Input:
  ##
  ##  X: n x N data matrix 
  ##  model_fit: Sparse_DFM fit
  
  n = dim(X)[1]
  N = dim(X)[2]
  
  Lambda = model_fit$Lambda
  factors = model_fit$factors.QML
  
  W = !is.na(X)
  
  loglik = 0 
  
  for(t in 1:n){
    
    new_Lambda <- Lambda[W[t,],]
    error <- as.numeric(na.omit(X[t,])) - new_Lambda %*% factors[t,]
    forb <- rowSums(t(error) %*% error)
    loglik = loglik + forb 
    
  }
  
  bic = log(loglik/(n*N)) + sum(as.vector(Lambda) != 0)*log(N*n)/(N*n)
  
  return(bic)
  
}
