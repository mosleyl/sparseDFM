#' @name summary.SparseDFM
#' @aliases print.SparseDFM
#' @aliases summary.SparseDFM
#' 
#' @title
#' SparseDFM Summary Outputs 
#' 
#' @description 
#' Summary and print outputs for class 'SparseDFM'.
#' 
#' @param x an object of class 'SparseDFM'
#' 
#' @returns 
#' Information on the model fitted.
#' 
#' @export

print.SparseDFM <- function(x){
  
  X = x$data$X
  A = x$params$A
  n = dim(X)[1]
  p = dim(X)[2]
  r = dim(A)[1]
  
  if(x$data$method=='PCA'){
    typeFM = 'Static'
  }else if(x$data$method=='EM-sparse'){
    typeFM = 'Sparse Dynamic'
  }else{
    typeFM = 'Dynamic'
  }
  
  cat('\n', typeFM, 'Factor Model using', x$data$method, 'with: \n\n n =', n,
      'observations, \n p =', p, 'variables, \n r =', r, 'factors, \n err =',
      x$data$err, 'idiosyncratic errors. \n\n Call summary() for estimation details.')
  
  
}

#' @rdname summary.SparseDFM
#' @param x an object of class 'SparseDFM'
#' @returns 
#' Summary information on estimation details. 
#' 
#' @export

summary.SparseDFM <- function(x){
    
  X = x$data$X
  A = x$params$A
  n = dim(X)[1]
  p = dim(X)[2]
  r = dim(A)[1]
  
  if(x$data$method=='PCA'){
    typeFM = 'Static'
  }else if(x$data$method=='EM-sparse'){
    typeFM = 'Sparse Dynamic'
  }else{
    typeFM = 'Dynamic'
  }
  
  cat(typeFM, 'Factor Model using', x$data$method, 'with: \n\n n =', n,
      'observations, \n p =', p, 'variables, \n r =', r, 'factors, \n err =',
      x$data$err)
  
  cat('\n\nCall: \n\n', paste(deparse(x$data$call)))
  
}






