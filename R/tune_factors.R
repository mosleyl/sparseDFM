#' Tune the number of factors to use 
#' 
#' Uses Bai and Ng (2002) informaion criteria approach. Missing data is interpolated using fill_NA() function. 
#' 
#' @param X A \code{n x p} numeric data matrix or data frame of (stationary) time series.
#' @param type Character. Option for which information criteria to use. Default is 2. 
#' @param standardize Logical. Standardize the data before estimating the model. Default is \code{TRUE}.
#' @param r.max Integer. Maximum number of factors to search for. Default is min(15,ncol(X)-1). 
#' @param plot Logical. Make a plot showing the IC value for each of the number of factors considered. Default is \code{FALSE}.
#'
#' @importFrom stats cov 
#' @importFrom graphics plot points 
#' 
#' @export


tune_factors <- function(X, type = 2, standardize = TRUE, r.max = min(15,ncol(X)-1), plot = FALSE){
  
  X = as.matrix(X)
  
  if(standardize){
    X = scale(X)
  }
  
  if(anyNA(X)){
    message("Data contains missing values: imputing data with fill_NA()")
    X = fill_NA(X)
    X = X$X
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
