#' Transform data to make it stationary 
#' 
#' @description 
#' Methods to transform the data to make it stationary. Input a \eqn{n \times p}{n x p} numeric data matrix and what transform is required for each data series. Returns a \eqn{n \times p}{n x p} matrix of the transformed data.  
#' 
#' @param X n x p numeric data matrix 
#' @param stationary_transform p-dimensional vector filled with numbers from \eqn{\{1,2,3,4,5,6,7\}}{\{1,2,3,4,5,6,7\}} representing:
#'    \tabular{llll}{
#' \code{1} \tab\tab no change \cr\cr
#' \code{2} \tab\tab first difference \eqn{X_{i,t} - X_{i,t-1}}{X_{i,t} - X_{i,t-1}} \cr\cr
#' \code{3} \tab\tab second difference \eqn{(X_{i,t} - X_{i,t-1}) - (X_{i,t-1} - X_{i,t-2})}{(X_{i,t} - X_{i,t-1}) - (X_{i,t-1} - X_{i,t-2})} \cr\cr
#' \code{4} \tab\tab log first difference \eqn{log(X_{i,t}) - log(X_{i,t-1})}{log(X_{i,t}) - log(X_{i,t-1})} \cr\cr
#' \code{5} \tab\tab log second difference \eqn{(log(X_{i,t}) - log(X_{i,t-1})) - (log(X_{i,t-1}) - log(X_{i,t-2}))}{(log(X_{i,t}) - log(X_{i,t-1})) - (log(X_{i,t-1}) - log(X_{i,t-2}))} \cr\cr
#' \code{6} \tab\tab growth rate \eqn{(X_{i,t} - X_{i,t-1})/X_{i,t-1}}{(X_{i,t} - X_{i,t-1})/X_{i,t-1}} \cr\cr
#' \code{7} \tab\tab log growth rate \eqn{(log(X_{i,t}) - log(X_{i,t-1}))/log(X_{i,t-1})}{(log(X_{i,t}) - log(X_{i,t-1}))/log(X_{i,t-1})} \cr\cr
#' }
#' 
#' @importFrom stats is.ts ts start frequency 
#' 
#' @returns 
#' Transformed stationary version of \eqn{\bm{X}}{X}. 
#' 
#' @export 

transformData <- function(X, stationary_transform){
  
  
  X = as.matrix(X)
  
  newX = matrix(NA, NROW(X), NCOL(X))
  
  for(i in 1:NCOL(newX)){
   
    transf = stationary_transform[i]
    
    # No change
    if(transf==1){
      newX[,i] = X[,i]
    }
    
    # First diff
    if(transf==2){
      newX[,i] = c(NA, diff(X[,i]))
    }
    
    # Second diff
    if(transf==3){
      newX[,i] = c(NA,NA, diff(diff(X[,i])))
    }
    
    # Log first diff
    if(transf==4){
      newX[,i] = c(NA, diff(log(X[,i])))
    }
    
    # Log second diff 
    if(transf==5){
      newX[,i] = c(NA,NA, diff(diff(log(X[,i]))))
    }
    
    # Growth rate 
    if(transf==6){
      newX[,i] = c(NA, diff(X[,i])/X[1:(NROW(X)-1),i])
    }
    
    # Log growth rate 
    if(transf==7){
      newX[,i] = c(NA, diff(log(X[,i]))/log(X[1:(NROW(X)-1),i]))
    }
    
  }
  colnames(newX) = colnames(X)
  
  if(is.ts(X)){
    
    newX = ts(newX, start = start(X), frequency = frequency(X))
    
  }
  
  
  return(newX)
}
