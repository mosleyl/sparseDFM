


Sparse_DFM <- function(X, r, p, alpha=NA, standardize = TRUE, max_iter=500, threshold=1e-4) {
  
  # Implements Sparse DFM via EM
  # Sparsity is on the factor loadings implemented in M-step
  
  # Inputs:
  #
  #   X: n x N matrix of (stationary) data 
  #   r: number of factors
  #   p: VAR(p) model for factors
  #   alpha: level of shrinkage on factors >=0 - default = NA 
  #   standardize: option to standardise data - default is TRUE  
  #   max_iter: maximum number of EM iterations
  #   threshold: tolerance on EM iterates
  
  library(Matrix)
  
  
  n = dim(X)[1]
  N = dim(X)[2]
  
  if(standardize){
    X = scale(X)
    X.mean = attr(X, "scaled:center")
    X.sd = attr(X, "scaled:scale")
  }
  
  ## Initialise with PCA and VAR(1)
  
  initialise <- initPCA(X,r,p)
  
  x_init = initialise$f_0
  V_init = initialise$V_0
  A = initialise$A
  Lambda = initialise$loadings.pca
  Sig_e = initialise$Sig_e
  Sig_u = initialise$Sig_u
  factors.PCA = initialise$factors.pca
  loadings.PCA = initialise$loadings.pca
  
  ## Kalman Filter and Smoother (Doz (2011) 2-step) 
  
  
  KFS <- kalmanCpp(X, x_init, V_init, A, Lambda, Sig_e, Sig_u)
  
  factors.KF = as.matrix(KFS$factors.KF[2:(n+1),])
  covariance.KF = KFS$covariance.KF[,,2:(n+1)] 
  factors.KS = as.matrix(KFS$factors.KS[2:(n+1),])
  covariance.KS = KFS$covariance.KS[,,2:(n+1)]
  
  
  ## EM Algorithm (Doz (2012) QML)
  
  previous_loglik = -.Machine$double.xmax
  loglik = 0
  num_iter = 0
  LL = c()
  converged = 0
  
  # data availability - n x p - 0 for missing, 1 for available 
  W = 1*!is.na(X)
  
  
  while ((num_iter < max_iter) & !converged) {
    
    
    # E-step - notation from Shumway and Stoffer (1982)
    KFS <- kalmanCpp(X, x_init, V_init, A, Lambda, Sig_e, Sig_u)
    xt_n = KFS$xt_n   # r x (n+1)
    Vt_n = KFS$Vt_n   # r x r x (n+1)
    Vt_tlag_n = KFS$Vt_tlag_n   # r x r x n
    loglik = KFS$logl
    
    num_iter = num_iter + 1
    
    
    if(r==1){
      P_t = t(xt_n[1:r,2:(n+1)]) %*% xt_n[1:r,2:(n+1)] + sum(Vt_n[1:r,1:r,2:(n+1)])
      P_tlag = t(xt_n[1:r,1:n]) %*% xt_n[1:r,1:n] + sum(Vt_n[1:r,1:r,1:n])
      Pt_tlag = t(xt_n[1:r,2:(n+1)]) %*% xt_n[1:r,1:n] + sum(Vt_tlag_n[1:r,1:r,])
      Ptlag_t = t(xt_n[1:r,1:n]) %*% xt_n[1:r,2:(n+1)] + sum(Vt_tlag_n[1:r,1:r,])
      P_0 = sum(Vt_n[1:r,1:r,1])
    }else{
      P_t = xt_n[1:r,2:(n+1)] %*% t(xt_n[1:r,2:(n+1)]) + apply(Vt_n[1:r,1:r,2:(n+1)],c(1,2),sum)
      P_tlag = xt_n[1:r,1:n] %*% t(xt_n[1:r,1:n]) + apply(Vt_n[1:r,1:r,1:n],c(1,2),sum)
      Pt_tlag = xt_n[1:r,2:(n+1)] %*% t(xt_n[1:r,1:n]) + apply(Vt_tlag_n[1:r,1:r,],c(1,2),sum)
      Ptlag_t = xt_n[1:r,1:n] %*% t(xt_n[1:r,2:(n+1)]) + apply(Vt_tlag_n[1:r,1:r,],c(1,2),sum)
      P_0 = apply(Vt_n[1:r,1:r,1],c(1,2),sum)
    }
    
    
    ##  M-step - not involving X
    
    x_init = as.matrix(xt_n[,1])
    V_init = as.matrix(P_0)
    A = Pt_tlag %*% solve(P_tlag) 
    Sig_u = (P_t - A %*% Ptlag_t)/n 
    
    ## M-step - involving X - deal with missing data 
    
    y = t(X)
    y[is.na(y)]=0
    
    # Loadings - Lambda update 
    Sig_e_inv <- 1/diag(Sig_e)
    At <- array(NA, dim=c(r, r, n))
    Bt <- matrix(NA, nrow = N, ncol = n)
    C_n = 0
    for(t in 1:n){
      At[,,t] = xt_n[1:r,t+1]%*%t(xt_n[1:r,t+1])+Vt_n[1:r,1:r,t+1]
      Bt[,t] = W[t,]*Sig_e_inv*W[t,]
      C_n = C_n + as.matrix(Bt[,t]*y[,t])%*%xt_n[1:r,t+1]
    }
    
    
    if(!(is.na(alpha))){
      
      D_cube = solveCube(At, Bt, nu = 1)
      sol = solveLambda(D_cube, C_n, nu=1, alpha=alpha)
      Lambda = sol$Z
      
      
    }else{
      
      D_cube = solveCube(At, Bt, nu = 0)
      Lambda = fastLambda(D_cube, C_n)
    }
    
    
    # Sig_e update - make diagonal (exact DFM assumed)
    
    Sig_e_new = 0
    I = rep(1,N)
    for(t in 1:n){
      
      Sig_e_new = Sig_e_new + (tcrossprod(W[t,]*y[,t]) - tcrossprod(W[t,]*y[,t],W[t,]*Lambda%*%xt_n[1:r,t+1])
                               - tcrossprod(W[t,]*Lambda%*%xt_n[1:r,t+1],W[t,]*y[,t]) + tcrossprod(W[t,]*Lambda%*%(At[,,t]),W[t,]*Lambda)
                               + diag((I - W[t,])*diag(Sig_e)*(I - W[t,])))
    }
    Sig_e = Sig_e_new/n
    Sig_e = diag(diag(Sig_e))


    #   # check convergence
    converged <- emConverged(loglik, previous_loglik, threshold)
    previous_loglik <- loglik
    LL <- c(LL, loglik)
    
  }
  
  if (converged == TRUE)
  {
    cat("Converged after", num_iter, "iterations.\n")
  } else
  {
    cat("Maximum number of iterations reached.\n")
  }
  
  # # final run of KFS
  KFS <- kalmanCpp(X, x_init, V_init, A, Lambda, Sig_e, Sig_u)
  factors.QML = t(matrix(KFS$xt_n[,2:(n+1)], nrow = r))
  covariance.QML = KFS$Vt_n[,,2:(n+1)]
  
  # # # full X matrix estimate (scaled and unscaled)
  fore_X = factors.QML %*% t(Lambda)
  if(standardize){
    fore_X = kronecker(t(X.sd),rep(1,n))*fore_X + kronecker(t(X.mean),rep(1,n))
  }
  
  output = list()
  output$factors.PCA = factors.PCA
  output$loadings.PCA = loadings.PCA
  output$factors.KF = factors.KF
  output$covariance.KF = covariance.KF
  output$factors.KS = factors.KS
  output$covariance.KS = covariance.KS
  output$factors.QML = factors.QML
  output$covariance.QML = covariance.QML
  output$fore_X = fore_X
  output$Lambda = Lambda
  output$A = A
  output$Sig_e = Sig_e
  output$Sig_u = Sig_u
  output$numIter = num_iter
  
  return(output)
  
  
}








