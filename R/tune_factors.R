# BAI AND NG (2002) INFORMATION CRITERIA FOR SELECTING THE NUMBER OF FACTORS 

tune_factors <- function(X, type = 2, standardize = TRUE, r.max = min(15,ncol(X)-1), plot = FALSE){
  
  if(standardize){
    X = scale(X)
  }
  
  if(anyNA(X)){
    message("Data contains missing values: imputing data with fill_NA()")
    X = fill_NA(X)
  }
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  evd = eigen(cov(X))
  loadings = evd$vectors 
  factors = X %*% loadings 
  
  IC = c()
  
  for(r in 1:r.max){
    
    eps = X - tcrossprod(factors[,1:r], loadings[,1:r])
    logV = log(sum(colSums(eps^2)/(n*p)))
    
    if(type == 1){
      IC[r] = logV + r*((n+p)/(n*p))*log((n*p)/(n+p))
    }else if(type == 2){
      IC[r] = logV + r*((n+p)/(n*p))*log(min(n,p))
    }else{
      IC[r] = logV + r*log(min(n,p))/(min(n,p))
    }
    
  }
  
  if(plot){
    graphics::plot(IC, main = paste("IC",type), xlab = "Number of factors", ylab = "Index")
    graphics::points(which.min(IC), IC[which.min(IC)], pch = 19, col ="red")
  }
  
  
  return(paste("The chosen number of factors using criteria type ", type, " is ", which.min(IC)))
  
  
}
