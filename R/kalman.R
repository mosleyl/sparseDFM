kalman <- function(X, a0_0, P0_0, A, Lambda, Sig_e, Sig_u) {
  
  ##  Kalman Filter and Smoother equations from Shumway and Stoffer (1982)
  ##
  ##  Works for a general state space model of the form:
  ##
  ##  X_t = Lambda*F_t + e_t,  e_t ~ N(0,Sig_e)
  ##  F_t = A*F_{t-1} + u_t,   u_t ~ N(0,Sig_u)

  ##  Inputs:
  ##
  ##  X: n x p, matrix of (stationary) time series 
  ##  a0_0: 1 x k, initial state mean vector 
  ##  P0_0: k x k, initial state covariance matrix
  ##  A: k x k, state matrix
  ##  Lambda: p x k, measurement matrix 
  ##  Sig_e: p x p, measurement equation residuals covariance matrix (diagonal)
  ##  Sig_u: k x k, state equation residuals covariance matrix  
  ##
  ##  NOTE: For the DFM with AR(1) errors we have that k = r + p where r is 
  ##        the number of factors and p is the number of variables.

  ## Initialise 
  
  n = dim(X)[1]
  p = dim(X)[2]
  k = dim(A)[1]
  
  at_tlag = matrix(0, k, n)         # state mean prediction
  Pt_tlag = array(0, c(k, k, n))    # state covariance prediction 
  at_t = matrix(0, k, n+1)          # state mean update 
  Pt_t = array(0, c(k, k, n+1))     # state covariance update 
  at_n = matrix(0, k, n+1)          # Kalman smoothed state mean
  Pt_n = array(0, c(k, k, n+1))     # Kalman smoothed state covariance 
  Pt_tlag_n = array(0, c(k, k, n))  # Kalman smoothed state covariance with lag
  
  at_t[,1] = a0_0     # initial state mean at t=0
  Pt_t[,,1] = P0_0    # initial state covariance at t=0

  W = !is.na(X)   # used to remove rows/cols in X, Lambda and Sig_e that are missing at t 
  y = t(X)        # work with p x n matrix y
  
  logl = 0        # log-likelihood required for convergence check in EM
  
  
  ## Kalman Filter loop for t = 1,...,n
  
  for(t in 1:n) {
    
    Lambda_t = Lambda[W[t,],]             # remove row of Lambda if X is missing at t 
    Sig_e_t = Sig_e[W[t,], W[t,]]         # remove row and column of Sig_e if X is missing at t
    y_t = as.numeric(na.omit(y[,t]))      # omit missing values at t 
    
    Sig_e_inv = diag(1/diag(Sig_e_t), ncol = length(diag(Sig_e_t)))   # Sig_e inverse 
    
    ## prediction equations 
    
    at_tlag[,t] = A %*% at_t[,t]
    Pt_tlag[,,t] = A %*% Pt_t[,,t] %*% t(A) + Sig_u
    
    ## calculating kalman gain from observed data 
    
    inov_res = y_t - Lambda_t %*% at_tlag[,t]
    inov_cov = Lambda_t %*% Pt_tlag[,,t] %*% t(Lambda_t) + Sig_e_t
    GG = t(Lambda_t) %*% Sig_e_inv %*% Lambda_t
    inov_cov_inv = Sig_e_inv - Sig_e_inv %*% Lambda_t %*% corpcor::pseudoinverse(diag(k) + Pt_tlag[,,t] %*% GG) %*% Pt_tlag[,,t] %*% t(Lambda_t) %*% Sig_e_inv   # uses Harvey (1990) Woodbury identity trick
    KG = Pt_tlag[,,t] %*% t(Lambda_t) %*% inov_cov_inv
    
    ## update equations
    
    at_t[,t+1] = at_tlag[,t] + KG %*% inov_res
    Pt_t[,,t+1] = Pt_tlag[,,t] - KG %*% Lambda_t %*% Pt_tlag[,,t]
    
    ## calculate log-likelihood 
    
    d = length(as.numeric(inov_res))
    inov_cov_det = prod(diag(Sig_e_t)) %*% det(diag(k) + Pt_tlag[,,t] %*% GG)  # again see Harvey (1990)
    denom = (2*pi)^(d/2)*sqrt(abs(inov_cov_det))
    mahal = rowSums(t(inov_res) %*% inov_cov_inv %*% inov_res)
    logl = logl + -0.5*mahal - log(denom)
    
    
  }
  
  
  ## output for kalman filter 
  
  factors.KF = t(at_t[,2:n+1])
  covariance.KF = Pt_t[,,2:n+1]
  
  ## initialise kalman smoother with t=n of filtered mean and covariance 

  at_n[,n+1] = at_t[,n+1]
  Pt_n[,,n+1] = Pt_t[,,n+1]
  
  ## initialise kalman smooth covariance with lag and smoother gain 
  
  Pt_tlag_n[,,n] = (diag(k) - KG %*% Lambda_t) %*% A %*% Pt_t[,,n]    
  
  J_2 = Pt_t[,,n] %*% t(A) %*% corpcor::pseudoinverse(Pt_tlag[,,n])  
  
  ## Kalman smoother loop for t = n,...,1
  
  for(t in n:1) {
    
    J_1 = J_2
    
    at_n[,t] = at_t[,t] + J_1 %*% (at_n[,t+1] - at_tlag[,t])
    Pt_n[,,t] = Pt_t[,,t] + J_1 %*% (Pt_n[,,t+1] - Pt_tlag[,,t]) %*% t(J_2)
    
    if(t > 1){
      
      J_2 <- Pt_t[,,t-1] %*% t(A) %*% corpcor::pseudoinverse(Pt_tlag[,,t-1]) 
      Pt_tlag_n[,,t-1] = Pt_t[,,t] %*% t(J_2) + J_1 %*% (Pt_tlag_n[,,t] - A %*% Pt_t[,,t]) %*% t(J_2)
      
    }
    
  }
  
  ## output for kalman smoother 
  
  factors.KS = t(at_n)
  covariance.KS = Pt_n
  
  
  ## output 
  
  output = list('at_n' = at_n, 'Pt_n' = Pt_n, 'Pt_tlag_n' = Pt_tlag_n, 'logl' = logl,
                'factors.KF' = factors.KF, 'covariance.KF' = covariance.KF, 
                'factors.KS' = factors.KS, 'covariance.KS' = covariance.KS)
  return(output)
  
}



