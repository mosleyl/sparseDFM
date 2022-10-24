softThresh <- function(X,thresh){
  # Takes matrix X and soft-thresholds entries
  
  return( apply(X,1:2,function(x) sign(x)*max(abs(x)-thresh,0)))
  
}