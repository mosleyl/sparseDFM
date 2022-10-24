rmse <- function(true, est){
  
  N = dim(true)[1]
  r = dim(true)[2]
  e = true - est 
  return(sqrt(psych::tr(t(e)%*%e)/(N*r)))
  
}

mae <- function(true, est){
  
  e = true - est
  return(mean(abs(e)))
  
}