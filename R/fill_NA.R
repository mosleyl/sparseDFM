#' Interpolation of missing data for initPCA
#'
#' Internal missing data is filled in using a cubic spline. 
#' Start and end of sample missing data is filled in using the median of the series and then
#' smoothed with a MA(3) process. 
#' 
#' @param X n x p numeric matrix of (stationary) time series 
#' 
#' @importFrom stats spline filter median 
#' 
#' @noRd


fill_NA <-function(X){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  k <- 3
  idx.na <- is.na(X)
  
  for (i in 1:p){  
    x = X[,i]
    na_x = is.na(x)
    t1 = min(which(!na_x))
    t2 = max(which(!na_x))
    
    x1 <- stats::spline(x[t1:t2],xout = 1:(t2-t1+1))
    xx <- x1$y
    x[t1:t2] <- x1$y
    na_x <- is.na(x)
    x[na_x] <- median(x,na.rm = T)
    
    x_MA3 <- stats::filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
    
    x_MA3 = x_MA3[(2*k+1):length(x_MA3)]
    x[idx.na[,i]] = x_MA3[idx.na[,i]]
    X[,i] = x
  }
  
  
  return(list(X = X, idx.na=idx.na))
  
}
