#' Estimate a Sparse Dynamic Factor Model
#'
#' Main function to allow estimation of a standard DFM or a sparse DFM (with sparse loadings) with:
#' an option for IID or AR(1) idiosyncratic errors
#' an option for Kalman Filter/Smoother estimation using standard multivariate equations or fast univariate filtering equations
#'
#' @param X A \code{n x p} numeric data matrix or data frame of (stationary) time series.
#' @param r Integer. Number of factors. 
#' @param q Integer. The first q series (columns of X) should not be made sparse. Default q = 0.
#' @param alphas Numeric vector or value of LASSO regularisation parameters. Default is alphas = logspace(-2,3,100).
#' @param alg Character. Option for estimation algorithm. Default is \code{"EM-sparse"}. Options are:
#'    \tabular{llll}{
#' \code{"PCA"} \tab\tab principle components analysis (PCA) for static factors seen in Stock and Watson (2002). \cr\cr
#' \code{"2Stage"} \tab\tab the two-stage framework of PCA plus Kalman filter/smoother seen in Giannone et al. (2008) and Doz et al. (2011). \cr\cr
#' \code{"EM"} \tab\tab the quasi-maximum likelihood approach using the EM algorithm to handle arbitrary patterns of missing data seen in Banbura and Modugno (2014). \cr\cr
#' \code{"EM-sparse"} \tab\tab the novel sparse EM approach allowing LASSO regularisation on factor loadings seen in (cite our paper). \cr\cr
#' }
#' @param err Character. Option for idiosyncratic errors. Default is \code{"IID"}. Options are:
#'    \tabular{llll}{
#' \code{"IID"} \tab\tab errors are IID white noise. \cr\cr
#' \code{"AR1"} \tab\tab errors follow an AR(1) process. \cr\cr
#' }
#' @param kalman Character. Option for Kalman filter and smoother equations. Default is \code{"univariate"}. Options are:
#'    \tabular{llll}{
#' \code{"multivariate"} \tab\tab classic Kalman filter and smoother equations seen in Shumway and Stoffer (1982). \cr\cr
#' \code{"univaraite"} \tab\tab univariate treatment (sequential processing) of the multivariate equations for fast Kalman filter and smoother seen in Koopman and Durbin (2000). \cr\cr
#' }
#' @param standardize Logical. Standardize the data before estimating the model. Default is \code{TRUE}.
#' @param max_iter Integer. Maximum number of EM iterations. Default is 100. 
#' @param threshold Numeric. Tolerance on EM iterates. Default is 1e-4. 
#'
#' @details 
#' Add details of the model here. 
#'  
#' @return item to be returned 
#'  
#' @useDynLib SparseDFM, .registration = TRUE
#'  
#' @export

SparseDFM <- function(X, r, q = 0, alphas = logspace(-2,3,100), alg = 'EM-sparse', err = 'IID', kalman = 'univariate', standardize = TRUE, max_iter=100, threshold=1e-4) {
  
  
  ## Correct input checks
  
  if(alg != 'PCA' && alg != '2Stage' && alg != 'EM' && alg != 'EM-sparse'){
    stop("Incorrect alg input")
  }
  if(err != 'IID' && err != 'AR1'){
    stop("Incorrect err input")
  }
  if(kalman != 'multivariate' && kalman != 'univariate'){
    stop("Incorrect kalman input")
  }
  if(!is.numeric(r) || r <= 0L){
    stop("r needs to be an integer > 0")
  }
  if(!is.numeric(q) || q < 0){
    stop("q needs to be an integer >= 0")
  }
  if(!is.numeric(alphas) || anyNA(alphas)){
    stop("alphas must be a numeric vector with no missing values")
  }
  if(!is.numeric(max_iter) || max_iter <= 0){
    stop("max_iter must be an integer > 0")
  }
  if(!is.numeric(threshold) || threshold <= 0){
    stop("threshold must be > 0")
  }
  
  
  # make sure lasso parameter grid is low to high 
  alphas = sort(alphas)
  
  # dimensions 
  n = dim(X)[1]
  p = dim(X)[2]
  k = r + p 
  
  # standardize if TRUE 
  X.raw = X
  X.scale = scale(X)
  X.mean = attr(X.scale, "scaled:center")
  X.sd = attr(X.scale, "scaled:scale")
  
  if(standardize){
    X = X.scale
  }
  
  
  ## Apply algorithm: PCA, 2Stage, EM, EM-sparse

  if(alg == 'PCA'){       # PCA algorithm applied 
  
    ## PCA and VAR(1) 
  
      initialise <- initPCA(X,r,err)
    
        a0_0 = initialise$a0_0
        P0_0 = initialise$P0_0
        A.tilde = initialise$A.tilde
        Lambda.tilde = initialise$Lambda.tilde
        Sigma.u.tilde = initialise$Sigma.u.tilde
        Sigma.eta = initialise$Sigma.eta
        factors.PCA = initialise$factors.pca
        loadings.PCA = initialise$loadings.pca
        
        fore_X = factors.PCA %*% t(loadings.PCA)
        if(standardize){
          fore_X = kronecker(t(X.sd),rep(1,n))*fore_X + kronecker(t(X.mean),rep(1,n))
        }
        
        errors.PCA = initialise$X.bal - factors.PCA %*% t(loadings.PCA)
        
    ## Output for PCA - depends on if err = 'AR1' or 'IID'
        
      if(err == 'AR1'){
        
        output = list(data = list(X = X.raw,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  predict = fore_X),
                      params = list(A = A.tilde[1:r,1:r],
                                    Phi = A.tilde[(r+1):k,(r+1):k],
                                    Lambda = Lambda.tilde[,1:r],
                                    Sigma_u = Sigma.u.tilde[1:r,1:r],
                                    Sigma_epsilon = Sigma.u.tilde[(r+1):k,(r+1):k]),
                      state = list(factors = factors.PCA,
                                   errors = errors.PCA,
                                   factors.cov = P0_0[1:r,1:r],
                                   errors.cov = P0_0[(r+1):k,(r+1):k]))
        
        return(output)
        
      }else {
        
        output = list(data = list(X = X.raw,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  predict = fore_X),
                      params = list(A = A.tilde,
                                    Lambda = Lambda.tilde,
                                    Sigma_u = Sigma.u.tilde,
                                    Sigma_epsilon = Sigma.eta),
                      state = list(factors = factors.PCA,
                                   factors.cov = P0_0[1:r,1:r]))
        
        return(output)
        
      }
        
  }else if(alg == '2Stage'){          # 2 stage algorithm applied (Doz, 2011)
    
    ## Initialise with PCA and VAR(1) 
    
      initialise <- initPCA(X,r,err)
    
        a0_0 = initialise$a0_0
        P0_0 = initialise$P0_0
        A.tilde = initialise$A.tilde
        Lambda.tilde = initialise$Lambda.tilde
        Sigma.u.tilde = initialise$Sigma.u.tilde
        Sigma.eta = initialise$Sigma.eta
        factors.PCA = initialise$factors.pca
        loadings.PCA = initialise$loadings.pca
  
    ## Kalman Filter and Smoother (Doz (2011) 2-step) 
  
      if(kalman == 'univariate'){
        KFS <- kalmanUnivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
      }else{
        KFS <- kalmanCpp(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
      }
        
        state.KF = t(as.matrix(KFS$at_t))
        covariance.KF = KFS$Pt_t 
        state.KS = t(as.matrix(KFS$at_n))
        covariance.KS = KFS$Pt_n
        
    ## Fill in missing data in X - KF and KS 
        
        fore_X_KF = state.KF %*% t(Lambda.tilde)
        if(standardize){
          fore_X_KF = kronecker(t(X.sd),rep(1,n))*fore_X_KF + kronecker(t(X.mean),rep(1,n))
        }
        
        fore_X_KS = state.KS %*% t(Lambda.tilde)
        if(standardize){
          fore_X_KS = kronecker(t(X.sd),rep(1,n))*fore_X_KS + kronecker(t(X.mean),rep(1,n))
        }
        
        
    ## Output for 2Stage - depends on if err = 'AR1' or 'IID'
    
      if(err == 'AR1'){
        
        output = list(data = list(X = X.raw,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  predict.KF = fore_X_KF,
                                  predict.KS = fore_X_KS),
                      params = list(A = A.tilde[1:r,1:r],
                                    Phi = A.tilde[(r+1):k,(r+1):k],
                                    Lambda = Lambda.tilde[,1:r],
                                    Sigma_u = Sigma.u.tilde[1:r,1:r],
                                    Sigma_epsilon = Sigma.u.tilde[(r+1):k,(r+1):k]),
                      state = list(factors.KF = state.KF[,1:r],
                                   errors.KF = state.KF[,(r+1):k],
                                   factors.KS = state.KS[,1:r],
                                   errors.KS = state.KS[,(r+1):k],
                                   factors.KF.cov = covariance.KF[1:r,1:r,],
                                   errors.KF.cov = covariance.KF[(r+1):k,(r+1):k,],
                                   factors.KS.cov = covariance.KS[1:r,1:r,],
                                   errors.KS.cov = covariance.KS[(r+1):k,(r+1):k,]))
        
        return(output)
        
      }else {
        
        output = list(data = list(X = X.raw,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  predict.KF = fore_X_KF,
                                  predict.KS = fore_X_KS),
                      params = list(A = A.tilde,
                                    Lambda = Lambda.tilde,
                                    Sigma_u = Sigma.u.tilde,
                                    Sigma_epsilon = Sigma.eta),
                      state = list(factors.KF = state.KF,
                                   factors.KS = state.KS,
                                   factors.KF.cov = covariance.KF,
                                   factors.KS.cov = covariance.KS))
        
        return(output)
        
      }
        
        
      
  }else if(alg == 'EM'){             # EM algorithm applied (Banbura and Modugno, 2014)
    
    ## Initialise with PCA and VAR(1) 
    
      initialise <- initPCA(X,r,err)
    
        a0_0 = initialise$a0_0
        P0_0 = initialise$P0_0
        A.tilde = initialise$A.tilde
        Lambda.tilde = initialise$Lambda.tilde
        Sigma.u.tilde = initialise$Sigma.u.tilde
        Sigma.eta = initialise$Sigma.eta
        factors.PCA = initialise$factors.pca
        loadings.PCA = initialise$loadings.pca
        
        
    ## EM Algorithm
          
      # EM iterations function
        EM.fit <- EM(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde, 
                   err = err, kalman = kalman, sparse = FALSE, max_iter = max_iter, threshold = threshold)
                   
      
      # Optimal parameters
        
        a0_0 = EM.fit$a0_0
        P0_0 = EM.fit$P0_0
        A.tilde = EM.fit$A.tilde
        Lambda.tilde = EM.fit$Lambda.tilde
        Sigma.u.tilde = EM.fit$Sigma.u.tilde
        Sigma.eta = EM.fit$Sigma.eta
        loglik.store = EM.fit$loglik.store
        converged = EM.fit$converged 
        num_iter = EM.fit$num_iter 
        
      # run KFS on final parameter estimates 
      
        if(kalman == 'univariate'){
          best.KFS <- kalmanUnivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
        }else{
          best.KFS <- kalmanCpp(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
        }
      
      state.EM = t(best.KFS$at_n)
      covariance.EM = best.KFS$Pt_n
      
      
      # fill in missing data in X
      
      fore_X = state.EM %*% t(Lambda.tilde)
      if(standardize){
        fore_X = kronecker(t(X.sd),rep(1,n))*fore_X + kronecker(t(X.mean),rep(1,n))
      }
      
      ## Output for EM - depends on if err = 'AR1' or 'IID'
      
      if(err == 'AR1'){
        
        output = list(data = list(X = X.raw,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  predict = fore_X),
                      params = list(A = A.tilde[1:r,1:r],
                                    Phi = A.tilde[(r+1):k,(r+1):k],
                                    Lambda = Lambda.tilde[,1:r],
                                    Sigma_u = Sigma.u.tilde[1:r,1:r],
                                    Sigma_epsilon = Sigma.u.tilde[(r+1):k,(r+1):k]),
                      state = list(factors = state.EM[,1:r],
                                   errors = state.EM[,(r+1):k],
                                   factors.cov = covariance.EM[1:r,1:r,],
                                   errors.cov = covariance.EM[(r+1):k,(r+1):k,]),
                      converged = list(converged = converged,
                                       loglik = loglik.store,
                                       num_iter = num_iter,
                                       tol = threshold,
                                       max_iter = max_iter))
        
        return(output)
        
      }else {
        
        output = list(data = list(X = X.raw,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  predict = fore_X),
                      params = list(A = A.tilde,
                                    Lambda = Lambda.tilde,
                                    Sigma_u = Sigma.u.tilde,
                                    Sigma_epsilon = Sigma.eta),
                      state = list(factors = state.EM,
                                   factors.cov = covariance.EM),
                      converged = list(converged = converged,
                                       loglik = loglik.store,
                                       num_iter = num_iter,
                                       tol = threshold,
                                       max_iter = max_iter))
        
        return(output)
        
      }

    
    
  }else {         # sparse EM algorithm applied (Mosley et al, 2022)
    
    
    ## Initialise with PCA and VAR(1) 
    
      initialise <- initPCA(X,r,err)
      
        a0_0 = initialise$a0_0
        P0_0 = initialise$P0_0
        A.tilde = initialise$A.tilde
        Lambda.tilde = initialise$Lambda.tilde
        Sigma.u.tilde = initialise$Sigma.u.tilde
        Sigma.eta = initialise$Sigma.eta
        factors.PCA = initialise$factors.pca
        loadings.PCA = initialise$loadings.pca
    
    
    ## Sparsified EM Algorithm
      
      ## Loop over alphas and calculate BIC until column of Lambda becomes 0 
        
        bic <- c()
        num_iter = c()
        best.bic <- .Machine$double.xmax
        
        for(alphas.index in 1:length(alphas)) {
          

          # lasso regularisation parameter 
            
            alpha.value = alphas[alphas.index]
          
          # make into a matrix 
           
            alpha.value = matrix(alpha.value, nrow = p, ncol = r)
            
          # alpha parameter set to 0 for no regularisation on first q series
            
           if(q > 0){
             alpha.value[1:q,] = 0
           }
          
          # EM iterations function
            
            EM.fit <- EM(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde, 
                       alpha.lasso = alpha.value, err = err, kalman = kalman, sparse = TRUE, 
                       max_iter = max_iter, threshold = threshold)
          
          # store number of iterations for each alpha
            
            num_iter[alphas.index] = EM.fit$num_iter
          
          # check if a column of Lambda has been set entirely to 0
            
            if(any(colSums(EM.fit$Lambda.tilde[(q+1):p,1:r]) == 0)){
              break
            }
          
          # Update parameters - used for warm start of the EM algorithm  
            
            a0_0 = EM.fit$a0_0
            P0_0 = EM.fit$P0_0
            A.tilde = EM.fit$A.tilde
            Lambda.tilde = EM.fit$Lambda.tilde 
            Sigma.u.tilde = EM.fit$Sigma.u.tilde
            Sigma.eta = EM.fit$Sigma.eta
          
          
          # run KFS on final parameter estimates  
            
            if(kalman == 'univariate'){
              KFS <- kalmanUnivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
            }else{
              KFS <- kalmanCpp(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
            }
          

          # calculate BIC 
            
            bic[alphas.index] = bic_function(X, t(KFS$at_n[1:r,]), Lambda.tilde[,1:r])
          
          
          # store estimates if BIC improved 
          
            if(bic[alphas.index] < best.bic){
            
            best.EM = EM.fit
            best.KFS = KFS
            best.bic = bic[alphas.index]
            
          }
          
          
        } 
        
        
      ## Store the optimal outputs  
      
        alphas.used = alphas[1:length(bic)]
        best.alpha = alphas[which.min(bic)]
        loglik.store = best.EM$loglik.store 
        converged = best.EM$converged 

        A.tilde = best.EM$A.tilde
        Lambda.tilde = best.EM$Lambda.tilde
        Sigma.u.tilde = best.EM$Sigma.u.tilde
        Sigma.eta = best.EM$Sigma.eta
        
        state.EM = t(best.KFS$at_n)
        covariance.EM = best.KFS$Pt_n
        
        
      ## Fill in missing data in X
        
        fore_X = state.EM %*% t(Lambda.tilde)
        if(standardize){
          fore_X = kronecker(t(X.sd),rep(1,n))*fore_X + kronecker(t(X.mean),rep(1,n))
        }
        
      ## Output for EM-sparse - depends on if err = 'AR1' or 'IID'
        
        if(err == 'AR1'){
          
          output = list(data = list(X = X.raw,
                                    standardize = standardize,
                                    X.mean = X.mean, 
                                    X.sd = X.sd,
                                    X.bal = initialise$X.bal,
                                    eigen = initialise$eigen,
                                    predict = fore_X),
                        params = list(A = A.tilde[1:r,1:r],
                                      Phi = A.tilde[(r+1):k,(r+1):k],
                                      Lambda = Lambda.tilde[,1:r],
                                      Sigma_u = Sigma.u.tilde[1:r,1:r],
                                      Sigma_epsilon = Sigma.u.tilde[(r+1):k,(r+1):k]),
                        state = list(factors = state.EM[,1:r],
                                     errors = state.EM[,(r+1):k],
                                     factors.cov = covariance.EM[1:r,1:r,],
                                     errors.cov = covariance.EM[(r+1):k,(r+1):k,]),
                        converged = list(converged = converged,
                                         alpha_grid = alphas.used,
                                         alpha_opt = best.alpha,
                                         bic = bic,
                                         loglik = loglik.store,
                                         num_iter = num_iter,
                                         tol = threshold,
                                         max_iter = max_iter))
          
          return(output)
          
        }else {
          
          output = list(data = list(X = X.raw,
                                    standardize = standardize,
                                    X.mean = X.mean, 
                                    X.sd = X.sd,
                                    X.bal = initialise$X.bal,
                                    eigen = initialise$eigen,
                                    predict = fore_X),
                        params = list(A = A.tilde,
                                      Lambda = Lambda.tilde,
                                      Sigma_u = Sigma.u.tilde,
                                      Sigma_epsilon = Sigma.eta),
                        state = list(factors = state.EM,
                                     factors.cov = covariance.EM),
                        converged = list(converged = converged,
                                         alpha_grid = alphas.used,
                                         alpha_opt = best.alpha,
                                         bic = bic,
                                         loglik = loglik.store,
                                         num_iter = num_iter,
                                         tol = threshold,
                                         max_iter = max_iter))
          
          return(output)
          
        }
        
    
  }
  
}
  
