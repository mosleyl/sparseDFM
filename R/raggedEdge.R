#' Generate a ragged edge structure for a data matrix 
#' 
#' @param X numeric data matrix 
#' @param lags vector of integers representing publication lag of each variable 
#' 
#' @return ragged edge version of X
#' 
#' @examples
#' data = matrix(rnorm(100),ncol=10)
#' pub_lags = c(rep(2,5),rep(1,3),rep(0,2))
#' new_data = raggedEdge(data, pub_lags)
#' 
#' @export


raggedEdge <- function(X, lags){
  
  Y = X
  nr = nrow(X)
  nc = ncol(X)
  
  for(i in 1:nc){
    lg = lags[i]
    if(lg>0){
      Y[(nr-lg+1):nr,i] = NA
    }
  }
  return(Y)
}