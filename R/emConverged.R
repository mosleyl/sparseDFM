emConverged <- function(loglik, previous_loglik, threshold=1e-4) {
  ### EM converge function
  
  converged <- FALSE
  
  delta_loglik <- abs(loglik - previous_loglik)
  avg_loglik <- (abs(loglik) + abs(previous_loglik) + .Machine$double.eps)/2
  
  if ((delta_loglik/avg_loglik) < threshold) {
    converged <- TRUE
  }
  return(converged)
  
}