# aggregation constraint matrix C for temporal disaggregation / benchmarking 
# M is the number of annual observations
# s is equal to 4 for quarters, 12 for months, etc...

agg.mat <- function(M,s) { 
  
  idmat <- diag(1,nrow=M)
  agvec <- rep(1,s)
  agmat <- t(kronecker(idmat,agvec))
  return(agmat)
  
}

agg.mat.diff <- function(M,s) {
  
  n = M*s + 2
  z = c(1,2,3,2,1,rep(0,n-5+3))
  Z = rep(z, M)
  M = matrix(Z[1:(n*M)], nrow = M, ncol = n, byrow = TRUE)
  
}

first_difference <- function(n) {
  diags <- list(rep(1, times = n), rep(-1, times = n-1))
  Delta_t <- bandSparse(n, k = 0:1, diagonals = diags, symmetric = FALSE)
  Delta <- t(Delta_t) 
  return(Delta)
}

Q_mat <- function(n,s){
  z = c(0,1,2,3,2,1, rep(0,3*(n-2)), rep(0,3))
  Z = rep(z, n-1)
  M = matrix(Z[1:(3*n*(n-1))], nrow = n-1, ncol = 3*n, byrow = T)
  Q = rbind(c(3,2,1, rep(0,3*n-3)), M)
  return(as.matrix(Q))
}
  

