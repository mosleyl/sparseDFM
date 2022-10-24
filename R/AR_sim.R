### Simulate an AR(1) process using this function. Make sure phi is strictly between -1 and 1 
### for stationarity. 


AR.sim <- function(n, phi, sigma) {
  
  x <- rep(0,n)
  eps <- rnorm(n,0,sigma)
  x[1] <- rnorm(1,0,sigma/(1-phi^2))
  
  for(t in 2:n) {
    x[t] <- phi*x[t-1] + eps[t]
  }
  return(x)
}



