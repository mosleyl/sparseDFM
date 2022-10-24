
initPCA <- function(X,r,p) {
  # Function to get initial parameter estimates 

  
  # fill in NA's using spline plus moving average (see Gianonne (2008))
  # rows with more than 80% of NAs are removed (I have removed this)
  # similar function used in nowcasting R package to initialise 
  data_fill = fill_NA(X)
  xBal<- data_fill$X
  idx.na <- data_fill$idx.na
  
  n <- dim(xBal)[1] # time n of series
  N <- dim(xBal)[2] # number of series
  
  # restore original missingness
  xNaN = xBal
  for(i in 1:N){
    xNaN[idx.na[,i],i] <- NA
  }

  A <- rbind(matrix(0, nrow=r, ncol=r*p),
             diag(1, nrow=r*(p-1), ncol=r*p))
  
  Sig_u <- matrix(0, nrow=p*r, ncol=p*r)
  Sig_u[1:r, 1:r] <- diag(1, r)
  
  
  evd = eigen(cov(xBal))
  loadings.pca = as.matrix(evd$vectors[,1:r])
  factors.pca = xBal %*% loadings.pca
  factors = factors.pca 
  
  
  inv_error = xBal - factors %*% t(loadings.pca)
  for(jj in 1:dim(idx.na)[2]){
    inv_error[idx.na[,jj],jj] <- NA
  }
  
  # covariance of measurement error without missing data
  Sig_e = diag(apply(inv_error, 2, stats::var, na.rm = T))
  
  fit <- VAR(factors, p)
  A[1:r, 1:(r*p)] <- t(fit$A)
  Sig_u <- cov(fit$res)
  
  f_0 = factors[1,]
  V_0 = matrix(corpcor::pseudoinverse(diag(dim(kronecker(A,A))[1])- kronecker(A,A)) %*% matrix(Sig_u, ncol = 1), r, r)
  loadings.pca <- cbind(loadings.pca, matrix(0, nrow=N, ncol=r*(p-1)))
  
  output = list('f_0' = f_0, 'V_0' = V_0, 'factors.pca' = factors.pca,
                'loadings.pca' = loadings.pca, 'A' = A, 'Sig_u' = Sig_u,
                'Sig_e' = Sig_e)
  
  return(output)
  
}
