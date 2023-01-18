#' Convergence of EM algorithm check function 
#' 
#' See Doz et al. (2011) for details 
#' 
#' @param loglik current log-likelihood of EM iteration
#' @param previous_loglik previous log-likelihood of EM iteration 
#' @param threshold EM algorithm comvergence threshold. Default is 1e-4. 
#' @noRd



emConverged <- function(loglik, previous_loglik, threshold=1e-4) {

  converged <- FALSE
  
  delta_loglik <- abs(loglik - previous_loglik)
  avg_loglik <- (abs(loglik) + abs(previous_loglik) + .Machine$double.eps)/2
  
  if ((delta_loglik/avg_loglik) < threshold) {
    converged <- TRUE
  }
  return(converged)
  
}