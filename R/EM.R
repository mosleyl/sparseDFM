#' Sparsified EM algorithm 
#'
#' E-step: Kalman filter and smoother using parameters estimated in M-step to obtain
#'        state mean, covariance and lagged-covariance. Log-likelihood calculated here.  
#' M-step: Maximisation of expected log-likelihood equations using E-step output. 
#'          With LASSO regularisation applied to the loadings matrix. 
#' Convergence: Uses the log-likelihood convergence rule from from Doz (2012)
#'
#' @param X n x p, (stationary) time series.
#' @param a0_0 k x 1, initial state mean vector 
#' @param P0_0 k x k, initial state covariance matrix
#' @param A.tilde k x k, initial state transition matrix
#' @param Lambda.tilde p x k, initial measurement matrix 
#' @param Sigma.eta p x p, initial measurement equation residuals covariance matrix (diagonal)
#' @param Sigma.u.tilde k x k, initial state equation residuals covariance matrix
#' @param alpha.lasso lasso tuning parameter value (>= 0). Default set to 0.1. 
#' @param err idiosyncratic error structure, 'AR1' or 'IID'. Default is 'IID'.
#' @param kalman KFS equations, 'multivariate' or 'univariate'. Default is 'univariate'. 
#' @param sparse if TRUE, lasso regularisation applied to the loadings. Default is TRUE. 
#' @param max_iter maximum number of iterations for the EM algorithm. Default is 500.
#' @param threshold threshold for log-likelihood convergence check. Default is 1e-4. 
#' @noRd


EM <- function(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde, alpha.lasso = 0.1, err = 'IID', kalman = 'univariate', sparse = TRUE, max_iter=100, threshold=1e-4) {

  
  if(err == 'AR1'){
    
    ## Initialise 
    
    n = dim(X)[1]
    p = dim(X)[2]
    k = dim(A.tilde)[1]
    r = k-p
    
    previous_loglik = -.Machine$double.xmax   # at least 2 iterations
    num_iter = 0        # counter begins at 0
    converged = 0       # convergence initialised to FALSE
    loglik.store = c()  # store log likelihoods 
    
    W = 1*!is.na(X)   # data availability - n x p - 0 for missing, 1 for observed
    
    
    ## EM iterations 
    
    while ((num_iter < max_iter) & !converged) {
      
      ## E-step 
      if(kalman == 'univariate'){
        KFS <- kalmanUnivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
      }else{
        KFS <- kalmanMultivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
      }
      
      at_n = KFS$at_n             # state mean: k x n matrix (t=1,...,n)
      Pt_n = KFS$Pt_n             # state covariance: k x k x n array (t=1,...,n) 
      Pt_tlag_n = KFS$Pt_tlag_n   # state covariance with lag: k x k x n array (t=1,...,n)
      
      loglik = KFS$logl           # log-likelihood value 

      num_iter = num_iter + 1     # new iteration
      
      
      ## M-step 
      
      ## Initial state mean and covariance update
      
      a0_0 = as.matrix(at_n[,1])
      P_F = as.matrix(Pt_n[1:r,1:r,1])
      P_E = diag(diag(as.matrix(Pt_n[(r+1):k,(r+1):k,1])))
      P0_0 = blkdiag(P_F, P_E)        
      
      ## State transition equation parameter updates: A.tilde and Sigma.u.tilde 
      
      ## required sums for the update equations
      
      if(r == 1){
        
        E_Ft_Ftlag = crossprod(at_n[1:r,2:n], at_n[1:r,1:(n-1)]) + sum(Pt_tlag_n[1:r,1:r,]) 
        E_Ftlag_Ftlag = crossprod(at_n[1:r,1:(n-1)], at_n[1:r,1:(n-1)]) + sum(Pt_n[1:r,1:r,1:(n-1)])
        E_Ft_Ft = crossprod(at_n[1:r,2:n], at_n[1:r,2:n]) + sum(Pt_n[1:r,1:r,2:n])
        E_Ftlag_Ft = crossprod(at_n[1:r,1:(n-1)], at_n[1:r,2:n]) + sum(Pt_tlag_n[1:r,1:r,])
        
      }else{
        
        E_Ft_Ftlag = tcrossprod(at_n[1:r,2:n], at_n[1:r,1:(n-1)]) + rowSums(Pt_tlag_n[1:r,1:r,], dims = 2L) 
        E_Ftlag_Ftlag = tcrossprod(at_n[1:r,1:(n-1)], at_n[1:r,1:(n-1)]) + rowSums(Pt_n[1:r,1:r,1:(n-1)], dims = 2L)
        E_Ft_Ft = tcrossprod(at_n[1:r,2:n], at_n[1:r,2:n]) + rowSums(Pt_n[1:r,1:r,2:n], dims = 2L)
        E_Ftlag_Ft = tcrossprod(at_n[1:r,1:(n-1)], at_n[1:r,2:n]) + rowSums(Pt_tlag_n[1:r,1:r,], dims = 2L)
        
      }
      
      E_et_etlag = diag(diag(tcrossprod(at_n[(r+1):k,2:n], at_n[(r+1):k,1:(n-1)])) + diag(rowSums(Pt_tlag_n[(r+1):k,(r+1):k,],dims = 2L)))
      E_etlag_etlag = diag(diag(tcrossprod(at_n[(r+1):k,1:(n-1)],at_n[(r+1):k,1:(n-1)])) + diag(rowSums(Pt_n[(r+1):k,(r+1):k,1:(n-1)],dims = 2L)))
      E_et_et = diag(diag(tcrossprod(at_n[(r+1):k,2:n],at_n[(r+1):k,2:n])) + diag(rowSums(Pt_n[(r+1):k,(r+1):k,2:n], dims = 2L)))
      E_etlag_et = diag(diag(tcrossprod(at_n[(r+1):k,1:(n-1)],at_n[(r+1):k,2:n])) + diag(rowSums(Pt_tlag_n[(r+1):k,(r+1):k,], dims = 2L)))
      
      ## A update 
      
      A = E_Ft_Ftlag %*% solve(E_Ftlag_Ftlag)
      
      ## Sig_u update 
      
      Sig_u = (E_Ft_Ft - A %*% E_Ftlag_Ft)/n 
      
      ## Phi update 
      
      Phi = E_et_etlag %*% diag(1 / diag(E_etlag_etlag))
      
      ## Sig_epsilon update 
      
      Sig_epsilon = (E_et_et - Phi %*% t(E_etlag_et)) / n 
      
      ## A.tilde construction 
      
      A.tilde = blkdiag(A, Phi)
      
      ## Sigma.u.tilde construction 
      
      Sigma.u.tilde = blkdiag(Sig_u, Sig_epsilon)
      
      
      ## Measurement equation parameter updates: Lambda.tilde and Sigma.eta 
      
      ## Sigma.eta update (just the same)
      
      kappa = 1e-4
      Sigma.eta = kappa*diag(p)
      
      ## Lambda.tilde update (ADMM algorithm for lasso regularisation - see paper)
      
      y = t(X)
      y[is.na(y)] = 0
      
      ## Calculate components A, B and C in the fast algorithm for Lambda computation (see paper)
      ##
      ## At: r x r x n array 
      ## Bt: p x p matrix 
      ## Ct: p x r matrix 
      
      At = array(NA, dim=c(r, r, n))
      Bt = t(W)
      Ct = 0
      
      for(t in 1:n){
        
        At[,,t] = at_n[1:r,t] %*% t(at_n[1:r,t]) + Pt_n[1:r,1:r,t]
        Ct = Ct + (as.matrix(Bt[,t]*y[,t]) %*% at_n[1:r,t]) - (diag(Bt[,t]) %*% (at_n[(r+1):k,t] %*% t(at_n[1:r,t]) + Pt_n[(r+1):k,1:r,t]))
        
      }
      
      ## Estimating loadings: with or without sparsity 
      
      if(sparse){
        
        ## Solve fast algorithm for Lambda 
        
        D_cube = solveCube(At, Bt, nu = kappa)
        sol = solveLambda(D_cube, Ct, nu = kappa, alpha = alpha.lasso)
        Lambda = sol$Z
        
        ## Lambda.tilde construction 
        
        Lambda.tilde = cbind(Lambda, diag(p))
        
      }else{
        
        D_cube = solveCube(At, Bt, nu = 0)
        Lambda = fastLambda(D_cube, Ct)
        
        Lambda.tilde = cbind(Lambda, diag(p))
        
      }
      
      
      ## Check convergence
      
      converged <- emConverged(loglik, previous_loglik, threshold)
      previous_loglik <- loglik
      loglik.store[num_iter] = loglik 
      
      
    }
    
    
    output <- list('a0_0' = a0_0, 'P0_0' = P0_0, 'A.tilde' = A.tilde, 'Sigma.u.tilde' = Sigma.u.tilde,
                   'Lambda.tilde' = Lambda.tilde, 'Sigma.eta' = Sigma.eta, 'loglik.store' = loglik.store,
                   'converged' = converged, 'num_iter' = num_iter)
    
    return(output)
    
  }else {         # IID case 
    
    ## Initialise 
    
    n = dim(X)[1]
    p = dim(X)[2]
    k = dim(A.tilde)[1]
    r = k
    
    previous_loglik = -.Machine$double.xmax   # at least 2 iterations
    num_iter = 0        # counter begins at 0
    converged = 0       # convergence initialised to FALSE
    loglik.store = c()  # store log likelihoods 
    
    W = 1*!is.na(X)   # data availability - n x p - 0 for missing, 1 for observed
    
    
    ## EM iterations 
    
    while ((num_iter < max_iter) & !converged) {
      
      ## E-step 
      if(kalman == 'univariate'){
        KFS <- kalmanUnivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
      }else{
        KFS <- kalmanMultivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
      }
      
      
      at_n = KFS$at_n             # state mean: k x n matrix (t=1,...,n)
      Pt_n = KFS$Pt_n             # state covariance: k x k x n array (t=1,...,n) 
      Pt_tlag_n = KFS$Pt_tlag_n   # state covariance with lag: k x k x n array (t=1,...,n)
      
      loglik = KFS$logl           # log-likelihood value 
      
      
      num_iter = num_iter + 1     # new iteration

      ## M-step 
      
      ## Initial state mean and covariance update
      
      a0_0 = as.matrix(at_n[,1])
      P0_0 = as.matrix(Pt_n[,,1])
      
      ## State transition equation parameter updates: A.tilde and Sigma.u.tilde 
      
      ## required sums for the update equations
      
      if(r == 1){
        
        E_Ft_Ftlag = crossprod(at_n[1:r,2:n], at_n[1:r,1:(n-1)]) + sum(Pt_tlag_n[1:r,1:r,]) 
        E_Ftlag_Ftlag = crossprod(at_n[1:r,1:(n-1)], at_n[1:r,1:(n-1)]) + sum(Pt_n[1:r,1:r,1:(n-1)])
        E_Ft_Ft = crossprod(at_n[1:r,2:n], at_n[1:r,2:n]) + sum(Pt_n[1:r,1:r,2:n])
        E_Ftlag_Ft = crossprod(at_n[1:r,1:(n-1)], at_n[1:r,2:n]) + sum(Pt_tlag_n[1:r,1:r,])
        
      }else{
        
        E_Ft_Ftlag = tcrossprod(at_n[1:r,2:n], at_n[1:r,1:(n-1)]) + rowSums(Pt_tlag_n[1:r,1:r,], dims = 2L) 
        E_Ftlag_Ftlag = tcrossprod(at_n[1:r,1:(n-1)], at_n[1:r,1:(n-1)]) + rowSums(Pt_n[1:r,1:r,1:(n-1)], dims = 2L)
        E_Ft_Ft = tcrossprod(at_n[1:r,2:n], at_n[1:r,2:n]) + rowSums(Pt_n[1:r,1:r,2:n], dims = 2L)
        E_Ftlag_Ft = tcrossprod(at_n[1:r,1:(n-1)], at_n[1:r,2:n]) + rowSums(Pt_tlag_n[1:r,1:r,], dims = 2L)
        
      }
      
      
      
      ## A update 
      
      A = E_Ft_Ftlag %*% solve(E_Ftlag_Ftlag)
      A.tilde = A
      
      ## Sig_u update 
      
      Sig_u = (E_Ft_Ft - A %*% E_Ftlag_Ft)/n 
      Sig.u.tilde = Sig_u
      
      
      ## Measurement equation parameter updates: Lambda.tilde and Sigma.eta 
      
      
      ## Lambda.tilde update (ADMM algorithm for lasso regularisation - see paper)
      
      y = t(X)
      y[is.na(y)] = 0
      
      ## Calculate components A, B and C in the fast algorithm for Lambda computation (see paper)
      ##
      ## At: r x r x n array 
      ## Bt: p x p matrix 
      ## Ct: p x r matrix 
      
      At = array(NA, dim=c(r, r, n))
      Bt = matrix(NA, nrow = p, ncol = n)
      Ct = 0
      Sigma.eta.inv = 1/diag(Sigma.eta)
      
      for(t in 1:n){
        
        At[,,t] = at_n[1:r,t] %*% t(at_n[1:r,t]) + Pt_n[1:r,1:r,t]
        Bt[,t] = W[t,]*Sigma.eta.inv*W[t,]
        Ct = Ct + as.matrix(Bt[,t]*y[,t]) %*% at_n[1:r,t]
        
      }
      
      ## Estimating loadings: with or without sparsity 
      
      if(sparse){
        
        ## Solve fast algorithm for Lambda 
        
        D_cube = solveCube(At, Bt, nu = 1)
        sol = solveLambda(D_cube, Ct, nu = 1, alpha = alpha.lasso)
        Lambda = sol$Z
        
        ## Lambda.tilde construction 
        
        Lambda.tilde = Lambda
        
      }else{
        
        D_cube = solveCube(At, Bt, nu = 0)
        Lambda = fastLambda(D_cube, Ct)
        
        Lambda.tilde = Lambda
        
      }
      
      
      
      ## Sigma.eta update 
      
      Sig_e_new = 0
      I = rep(1,p)
      for(t in 1:n){
        
        Sig_e_new = Sig_e_new + (tcrossprod(W[t,]*y[,t]) - tcrossprod(W[t,]*y[,t],W[t,]*Lambda%*%at_n[1:r,t])
                                 - tcrossprod(W[t,]*Lambda%*%at_n[1:r,t],W[t,]*y[,t]) + tcrossprod(W[t,]*Lambda%*%(At[,,t]),W[t,]*Lambda)
                                 + diag((I - W[t,])*diag(Sigma.eta)*(I - W[t,])))
      }
      Sig_e = Sig_e_new/n
      Sig_e = diag(Sig_e)
      Sig_e[Sig_e < 1e-7] = 1e-7
      Sig_e = diag(Sig_e)
      
      Sigma.eta = Sig_e
      converged <- emConverged(loglik, previous_loglik, threshold)
      previous_loglik <- loglik
      loglik.store[num_iter] = loglik
      
    }
    
    output <- list('a0_0' = a0_0, 'P0_0' = P0_0, 'A.tilde' = A.tilde, 'Sigma.u.tilde' = Sigma.u.tilde,
                   'Lambda.tilde' = Lambda.tilde, 'Sigma.eta' = Sigma.eta, 'loglik.store' = loglik.store,
                   'converged' = converged, 'num_iter' = num_iter)
    
    return(output)
    
    
  }
  
}

