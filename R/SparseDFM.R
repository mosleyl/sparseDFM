#' Estimate a Sparse Dynamic Factor Model
#'
#' @description
#' Main function to allow estimation of a DFM or a sparse DFM (with sparse loadings) on stationary data that may have arbitrary patterns of missing data. We allow the user:
#' \itemize{
#' \item an option for estimation method - \code{"PCA"}, \code{"2Stage"}, \code{"EM"} or \code{"EM-sparse"}
#' \item an option for \code{IID} or \code{AR1} idiosyncratic errors
#' \item an option for Kalman Filter/Smoother estimation using standard \code{multivariate} equations or fast \code{univariate} filtering equations
#' }  
#'
#' @param X \code{n x p} numeric data matrix or data frame of (stationary) time series.
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
#' @param store.parameters Logical. Store outputs for every alpha L1 penalty parameter. Default is FALSE.
#' @param standardize Logical. Standardize the data before estimating the model. Default is \code{TRUE}.
#' @param max_iter Integer. Maximum number of EM iterations. Default is 100. 
#' @param threshold Numeric. Tolerance on EM iterates. Default is 1e-4. 
#'
#' @details 
#' For full details of the model please refer to Mosley et al. (2023).
#'  
#' @returns A list-of-lists-like object of class 'sparseDFM' with the following elements:
#'  \item{\code{data}}{A list containing information about the data with the following elements:
#'  \tabular{llll}{
#'      \code{X} \tab\tab is the original \eqn{n \times p}{n x p} numeric data matrix of (stationary) time series. \cr\cr
#'      \code{standardize} \tab\tab is a logical value indicating whether the original data was standardized.\cr\cr
#'      \code{X.mean} \tab\tab is a p-dimensional numeric vector of column means of \eqn{X}.  \cr\cr
#'      \code{X.sd} \tab\tab is a p-dimensional numeric vector of column standard deviations of \eqn{X}.  \cr\cr
#'      \code{X.bal} \tab\tab is a \eqn{n \times p}{n x p} numeric data matrix of the original \eqn{X} with missing data interpolated using \code{fillNA()}. \cr\cr
#'      \code{eigen} \tab\tab is the eigen decomposition of \code{X.bal}. \cr\cr 
#'      \code{fitted} \tab\tab is the \eqn{n \times p}{n x p} predicted data matrix using the estimated parameters: \eqn{\hat{\Lambda}\hat{F}}{\hat{\Lambda}\hat{F}} for IID errors and \eqn{\hat{\Lambda}\hat{F}+\hat{\epsilon}}{\hat{\Lambda}\hat{F}+\hat{\epsilon}} for AR(1) errors. \cr\cr
#'      \code{fitted.unscaled} \tab\tab is the \eqn{n \times p}{n x p} predicted data matrix using the estimated parameters: \eqn{\hat{\Lambda}\hat{F}}{\hat{\Lambda}\hat{F}} for IID errors and \eqn{\hat{\Lambda}\hat{F}+\hat{\epsilon}}{\hat{\Lambda}\hat{F}+\hat{\epsilon}} for AR(1) errors. This has been unscaled back to original data scale if \code{standardize} is \code{TRUE}. \cr\cr
#'      \code{method} \tab\tab the estimation algorithm used (\code{alg}). \cr\cr
#'      \code{err} \tab\tab the type of idiosyncratic errors assumed. Either \code{IID} or \code{AR1}. \cr\cr
#'      \code{call} \tab\tab call object obtained from \code{match.call()}. \cr\cr
#'      }
#'  }
#'  \item{\code{params}}{A list containing the estimated parameters of the model with the following elements:
#'  \tabular{llll}{
#'      \code{A} \tab\tab the \eqn{r \times r}{r x r} factor transition matrix. \cr\cr
#'      \code{Phi} \tab\tab the p-dimensional vector of AR(1) coefficients for the idiosyncratic errors. \cr\cr
#'      \code{Lambda} \tab\tab the \eqn{p \times r}{p x r} loadings matrix. \cr\cr
#'      \code{Sigma_u} \tab\tab the \eqn{r \times r}{r x r} factor transition error covariance matrix. \cr\cr
#'      \code{Sigma_epsilon} \tab\tab the p-dimensional vector of idiosyncratic error variances. As \eqn{\bm{\Sigma}_{\epsilon}}{\Sigma(\epsilon)} is assumed to be diagonal. \cr\cr  
#'      }
#'  }
#'  \item{\code{state}}{A list containing the estimated states and state covariances with the following elements:
#'  \tabular{llll}{
#'      \code{factors} \tab\tab the \eqn{n \times r}{n x r} matrix of factor estimates. \cr\cr
#'      \code{errors} \tab\tab the \eqn{n \times p}{n x p} matrix of AR(1) idiosyncratic error estimates. For err = AR1 only. \cr\cr
#'      \code{factors.cov} \tab\tab the \eqn{r \times r \times n}{r x r x n} covariance matrices of the factor estimates. \cr\cr
#'      \code{errors.cov} \tab\tab the \eqn{p \times p \times n}{p x p x n} covariance matrices of the AR(1) idiosyncratic error estimates. For err = AR1 only. \cr\cr
#'      }
#'  }
#'  \item{\code{em}}{A list containing information about the EM algorithm with the following elements:
#'  \tabular{llll}{
#'      \code{converged} \tab\tab a logical value indicating whether the EM algorithm converged. \cr\cr
#'      \code{alpha_grid} \tab\tab a numerical vector containing the LASSO tuning parameters considered in BIC evaluation before stopping. \cr\cr
#'      \code{alpha_opt} \tab\tab the optimal LASSO tuning parameter used. \cr\cr
#'      \code{bic} \tab\tab a numerical vector containing BIC values for the corresponding LASSO tuning parameter in \code{alpha_grid}. \cr\cr
#'      \code{loglik} \tab\tab the log-likelihood of the innovations from the Kalman filter in the final model. \cr\cr
#'      \code{num_iter} \tab\tab number of iterations taken by the EM algorithm. \cr\cr
#'      \code{tol} \tab\tab tolerance for EM convergence. Matches \code{threshold} in the input. \cr\cr
#'      \code{max_iter} \tab\tab maximum number of iterations allowed for the EM algorithm. Matches \code{max_iter} in the input. \cr\cr
#'      \code{em_time} \tab\tab time taken for EM convergence \cr\cr
#'      }
#'  \item{\code{alpha.output}}{Parameter and state outputs for each L1-norm penalty parameter in \code{alphas} if \code{store.parameters = TRUE}. \cr\cr}
#'  }
#'  
#' @references 
#' Banbura, M., & Modugno, M. (2014). Maximum likelihood estimation of factor models on datasets with arbitrary pattern of missing data. \emph{Journal of Applied Econometrics, 29}(1), 133-160.
#'
#' Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for large approximate dynamic factor models based on Kalman filtering. \emph{Journal of Econometrics, 164}(1), 188-205.
#' 
#' Giannone, D., Reichlin, L., & Small, D. (2008). Nowcasting: The real-time informational content of macroeconomic data. \emph{Journal of monetary economics, 55}(4), 665-676.
#' 
#' Koopman, S. J., & Durbin, J. (2000). Fast filtering and smoothing for multivariate state space models. \emph{Journal of Time Series Analysis, 21}(3), 281-296. 
#' 
#' Mosley, L., Chan, TS., & Gibberd, A. (2023). sparseDFM: An R Package to Estimate Dynamic Factor Models with Sparse Loadings.  
#' 
#' Shumway, R. H., & Stoffer, D. S. (1982). An approach to time series smoothing and forecasting using the EM algorithm. \emph{Journal of time series analysis, 3}(4), 253-264.
#' 
#' Stock, J. H., & Watson, M. W. (2002). Forecasting using principal components from a large number of predictors. \emph{Journal of the American statistical association, 97}(460), 1167-1179.  
#'  
#' @examples 
#' # load inflation data set 
#' data = inflation
#' 
#' # make stationary by taking first differences 
#' new_data = transformData(data, rep(2,ncol(data)))
#' 
#' # tune for the number of factors to use 
#' tuneFactors(new_data, type = 2)
#' 
#' # fit a PCA using 3 PC's
#' fit.pca <- sparseDFM(new_data, r = 3, alg = 'PCA')
#' 
#' # fit a DFM using the two-stage approach 
#' fit.2stage <- sparseDFM(new_data, r = 3, alg = '2Stage')
#' 
#' # fit a DFM using EM algorithm with 3 factors 
#' fit.dfm <- sparseDFM(new_data, r = 3, alg = 'EM')
#' 
#' # fit a Sparse DFM with 3 factors 
#' fit.sdfm <- sparseDFM(new_data, r = 3, alg = 'EM-sparse')
#' 
#' # observe the factor loadings of the sparse DFM
#' plot(fit.sdfm, type = 'loading.heatmap')
#' 
#' # observe the factors 
#' plot(fit.sdfm, type = 'factor')
#' 
#' # observe the residuals 
#' plot(fit.sdfm, type = 'residual')
#' 
#' # observe the LASSO parameter selected and BIC values 
#' plot(fit.sdfm, type = 'lasso.bic')
#' 
#' # predict 3 steps ahead 
#' predict(fit.sdfm, h = 3)
#' 
#'  
#' @useDynLib sparseDFM, .registration = TRUE
#'  
#' @export

sparseDFM <- function(X, r, q = 0, alphas = logspace(-2,3,100), alg = 'EM-sparse', err = 'IID', kalman = 'univariate', store.parameters = FALSE, standardize = TRUE, max_iter=100, threshold=1e-4) {
  
  
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
  
  # return original X input 
  X.input = X
  
  # Unclass X - make sure it is a numeric matrix 
  X = unclass(X)
  
  # make sure lasso parameter grid is low to high 
  alphas = sort(alphas)
  
  # dimensions 
  n = dim(X)[1]
  p = dim(X)[2]
  k = r + p 
  
  # standardize if TRUE 
  X.scale = scale(X)
  X.mean = attr(X.scale, "scaled:center")
  X.sd = attr(X.scale, "scaled:scale")
  
  if(standardize){
    X = X.scale
  }
  
  # Remove any columns that are entirely NA
  numberNA = as.numeric(colSums(is.na(X)))
  allNA = which(numberNA == n)
  if(length(allNA) > 0){
    X = X[,-allNA]
    X.input = X.input[,-allNA]
    p = dim(X)[2]
    X.mean = X.mean[-allNA]
    X.sd = X.sd[-allNA]
    message('Columns: ',allNA, ' are entirely missing. Removing these columns.')
  }
  
                        
  
  # label variables 
  obs.names = unlist(dimnames(X)[1])
  series.names = unlist(dimnames(X)[2])
  factor.names = paste0('F',1:r)
  f.error.names = paste0('u',1:r)
  
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
        
        fit_x = factors.PCA %*% t(loadings.PCA)
        if(standardize){
          fit_X = kronecker(t(X.sd),rep(1,n))*fit_x + kronecker(t(X.mean),rep(1,n))
        }else{
          fit_X = fit_x
        }
        
        dimnames(fit_x) = dimnames(X)
        dimnames(fit_X) = dimnames(X)
        
        dimnames(factors.PCA) = list(obs.names, factor.names)
      
        
    ## Output for PCA - depends on if err = 'AR1' or 'IID'
        
      if(err == 'AR1'){
        
        A = A.tilde[1:r,1:r]
        Phi = A.tilde[(r+1):k,(r+1):k]
        Lambda = Lambda.tilde[,1:r]
        Sigma_u = Sigma.u.tilde[1:r,1:r]
        Sigma_epsilon = Sigma.u.tilde[(r+1):k,(r+1):k]
        factors.cov = P0_0[1:r,1:r]

        dimnames(A) = list(factor.names,factor.names)
        dimnames(Phi) = list(series.names,series.names)
        dimnames(Lambda) = list(series.names, factor.names)
        dimnames(Sigma_u) = list(f.error.names, f.error.names)
        dimnames(Sigma_epsilon) = list(series.names, series.names)
        dimnames(factors.cov) = list(factor.names, factor.names)

        
        output = list(data = list(X = X.input,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  fitted = fit_x,
                                  fitted.unscaled = fit_X,
                                  method = alg,
                                  err = err,
                                  call = match.call()),
                      params = list(a0_0 = a0_0,
                                    P0_0 = P0_0,
                                    A = A,
                                    Phi = diag(Phi),
                                    Lambda = Lambda,
                                    Sigma_u = Sigma_u,
                                    Sigma_epsilon = diag(Sigma_epsilon)),
                      state = list(factors = factors.PCA,
                                   factors.cov = factors.cov))
        
        
      }else {
        
        A = A.tilde
        Lambda = Lambda.tilde
        Sigma_u = Sigma.u.tilde
        Sigma_epsilon = Sigma.eta
        factors.cov = P0_0[1:r,1:r]

        dimnames(A) = list(factor.names,factor.names)
        dimnames(Lambda) = list(series.names, factor.names)
        dimnames(Sigma_u) = list(f.error.names, f.error.names)
        dimnames(Sigma_epsilon) = list(series.names, series.names)
        dimnames(factors.cov) = list(factor.names, factor.names)

        output = list(data = list(X = X.input,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  fitted = fit_x,
                                  fitted.unscaled = fit_X,
                                  method = alg,
                                  err = err,
                                  call = match.call()),
                      params = list(a0_0 = a0_0,
                                    P0_0 = P0_0,
                                    A = A,
                                    Lambda = Lambda,
                                    Sigma_u = Sigma_u,
                                    Sigma_epsilon = diag(Sigma_epsilon)),
                      state = list(factors = factors.PCA,
                                   factors.cov = factors.cov))
        
        
      }
        
      class(output) <- 'sparseDFM'
      return(output)
        
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
        KFS <- kalmanMultivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
      }
        
        state.KS = t(as.matrix(KFS$at_n))
        covariance.KS = KFS$Pt_n
        
    ## Fill in missing data in X
        
        fit_x = state.KS[,1:r] %*% t(Lambda.tilde[,1:r])
        if(standardize){
          fit_X = kronecker(t(X.sd),rep(1,n))*fit_x + kronecker(t(X.mean),rep(1,n))
        }else{
          fit_X = fit_x
        }
        
        dimnames(fit_x) = dimnames(X)
        dimnames(fit_X) = dimnames(X)
        
    ## Output for 2Stage - depends on if err = 'AR1' or 'IID'
    
      if(err == 'AR1'){
        
        factors = state.KS[,1:r]
        errors = state.KS[,(r+1):k]
        
        dimnames(factors) = list(obs.names, factor.names)
        dimnames(errors) = list(obs.names, series.names)
        
        
        
        A = A.tilde[1:r,1:r]
        Phi = A.tilde[(r+1):k,(r+1):k]
        Lambda = Lambda.tilde[,1:r]
        Sigma_u = Sigma.u.tilde[1:r,1:r]
        Sigma_epsilon = Sigma.u.tilde[(r+1):k,(r+1):k]
        factors.cov = covariance.KS[1:r,1:r,]
        errors.cov = covariance.KS[(r+1):k,(r+1):k,]
        
        dimnames(A) = list(factor.names,factor.names)
        dimnames(Phi) = list(series.names,series.names)
        dimnames(Lambda) = list(series.names, factor.names)
        dimnames(Sigma_u) = list(f.error.names, f.error.names)
        dimnames(Sigma_epsilon) = list(series.names, series.names)
        dimnames(factors.cov) = list(factor.names, factor.names)
        dimnames(errors.cov) = list(series.names, series.names)
        
        output = list(data = list(X = X.input,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  fitted = fit_x,
                                  fitted.unscaled = fit_X,
                                  method = alg,
                                  err = err,
                                  call = match.call()),
                      params = list(a0_0 = a0_0,
                                    P0_0 = P0_0,
                                    A = A,
                                    Phi = diag(Phi),
                                    Lambda = Lambda,
                                    Sigma_u = Sigma_u,
                                    Sigma_epsilon = diag(Sigma_epsilon)),
                      state = list(factors = factors,
                                   errors = errors,
                                   factors.cov = factors.cov,
                                   errors.cov = errors.cov))
                                   
        
        
      }else {
        
        factors = state.KS

        dimnames(factors) = list(obs.names, factor.names)

        A = A.tilde
        Lambda = Lambda.tilde
        Sigma_u = Sigma.u.tilde
        Sigma_epsilon = Sigma.eta
        factors.cov = covariance.KS
        
        dimnames(A) = list(factor.names,factor.names)
        dimnames(Lambda) = list(series.names, factor.names)
        dimnames(Sigma_u) = list(f.error.names, f.error.names)
        dimnames(Sigma_epsilon) = list(series.names, series.names)
        dimnames(factors.cov) = list(factor.names, factor.names)
        
        output = list(data = list(X = X.input,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  fitted = fit_x,
                                  fitted.unscaled = fit_X,
                                  method = alg,
                                  err = err,
                                  call = match.call()),
                      params = list(a0_0 = a0_0,
                                    P0_0 = P0_0,
                                    A = A,
                                    Lambda = Lambda,
                                    Sigma_u = Sigma_u,
                                    Sigma_epsilon = diag(Sigma_epsilon)),
                      state = list(factors = factors,
                                   factors.cov = factors.cov))
        
        
      }
        
      class(output) <- 'sparseDFM'
      return(output)
      
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
        st.em = Sys.time()
        EM.fit <- EM(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde, 
                   err = err, kalman = kalman, sparse = FALSE, max_iter = max_iter, threshold = threshold)
        et.em = Sys.time()
        time.em = et.em - st.em
                   
      
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
          best.KFS <- kalmanMultivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
        }
      
      state.EM = t(best.KFS$at_n)
      covariance.EM = best.KFS$Pt_n
      
      
      # fill in missing data in X
      
      fit_x = state.EM[,1:r] %*% t(Lambda.tilde[,1:r])
      if(standardize){
        fit_X = kronecker(t(X.sd),rep(1,n))*fit_x + kronecker(t(X.mean),rep(1,n))
      }else{
        fit_X = fit_x
      }
      
      dimnames(fit_x) = dimnames(X)
      dimnames(fit_X) = dimnames(X)
      
      ## Output for EM - depends on if err = 'AR1' or 'IID'
      
      if(err == 'AR1'){
        
        factors = state.EM[,1:r]
        errors = state.EM[,(r+1):k]
        
        dimnames(factors) = list(obs.names, factor.names)
        dimnames(errors) = list(obs.names, series.names)
        
        A = A.tilde[1:r,1:r]
        Phi = A.tilde[(r+1):k,(r+1):k]
        Lambda = Lambda.tilde[,1:r]
        Sigma_u = Sigma.u.tilde[1:r,1:r]
        Sigma_epsilon = Sigma.u.tilde[(r+1):k,(r+1):k]
        factors.cov = covariance.EM[1:r,1:r,]
        errors.cov = covariance.EM[(r+1):k,(r+1):k,]
        
        dimnames(A) = list(factor.names,factor.names)
        dimnames(Phi) = list(series.names,series.names)
        dimnames(Lambda) = list(series.names, factor.names)
        dimnames(Sigma_u) = list(f.error.names, f.error.names)
        dimnames(Sigma_epsilon) = list(series.names, series.names)
        dimnames(factors.cov) = list(factor.names, factor.names)
        dimnames(errors.cov) = list(series.names, series.names)
        
        output = list(data = list(X = X.input,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  fitted = fit_x,
                                  fitted.unscaled = fit_X,
                                  method = alg,
                                  err = err,
                                  call = match.call()),
                      params = list(a0_0 = a0_0,
                                    P0_0 = P0_0,
                                    A = A,
                                    Phi = diag(Phi),
                                    Lambda = Lambda,
                                    Sigma_u = Sigma_u,
                                    Sigma_epsilon = diag(Sigma_epsilon)),
                      state = list(factors = factors,
                                   errors = errors,
                                   factors.cov = factors.cov,
                                   errors.cov = errors.cov),
                      em = list(converged = converged,
                                loglik = loglik.store,
                                num_iter = num_iter,
                                tol = threshold,
                                max_iter = max_iter,
                                em_time = time.em))
        
        
      }else {
        
        factors = state.EM
        
        dimnames(factors) = list(obs.names, factor.names)
        
        A = A.tilde
        Lambda = Lambda.tilde
        Sigma_u = Sigma.u.tilde
        Sigma_epsilon = Sigma.eta
        factors.cov = covariance.EM
        
        dimnames(A) = list(factor.names,factor.names)
        dimnames(Lambda) = list(series.names, factor.names)
        dimnames(Sigma_u) = list(f.error.names, f.error.names)
        dimnames(Sigma_epsilon) = list(series.names, series.names)
        dimnames(factors.cov) = list(factor.names, factor.names)
        
        output = list(data = list(X = X.input,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  fitted = fit_x,
                                  fitted.unscaled = fit_X,
                                  method = alg,
                                  err = err,
                                  call = match.call()),
                      params = list(a0_0 = a0_0,
                                    P0_0 = P0_0,
                                    A = A,
                                    Lambda = Lambda,
                                    Sigma_u = Sigma_u,
                                    Sigma_epsilon = diag(Sigma_epsilon)),
                      state = list(factors = factors,
                                   factors.cov = factors.cov),
                      em = list(converged = converged,
                                loglik = loglik.store,
                                num_iter = num_iter,
                                tol = threshold,
                                max_iter = max_iter,
                                em_time = time.em))
        
        
      }

      class(output) <- 'sparseDFM'
      return(output)
    
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
        alpha.output = list()
        time.emsparse = c()
        
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
            
            st.emsparse = Sys.time()
            EM.fit <- EM(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde, 
                       alpha.lasso = alpha.value, err = err, kalman = kalman, sparse = TRUE, 
                       max_iter = max_iter, threshold = threshold)
            et.emsparse = Sys.time()
            time.emsparse[alphas.index] = et.emsparse - st.emsparse
          
          # store number of iterations for each alpha
            
            num_iter[alphas.index] = EM.fit$num_iter
            
          
          # check if a column of Lambda has been set entirely to 0
            
            if(any(colSums(EM.fit$Lambda.tilde[(q+1):p,1:r]) == 0) && alphas.index == 1){
              stop('First alpha value to large. All variables set to 0 in a factor. Make alpha value smaller.')
            }
            
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
              KFS <- kalmanMultivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
            }
          

          # calculate BIC 
            
            bic[alphas.index] = bicFunction(X, t(KFS$at_n[1:r,]), Lambda.tilde[,1:r])
            
            
          # Store parameters for each alpha if store.parameters = TRUE
          
          if(store.parameters){
            
            store.state = t(KFS$at_n)
            store.cov = KFS$Pt_n
            
            if(err == 'AR1'){
              
              store.factors = store.state[,1:r]
              store.errors = store.state[,(r+1):k]
              
              dimnames(store.factors) = list(obs.names, factor.names)
              dimnames(store.errors) = list(obs.names, series.names)
              
              store.A = A.tilde[1:r,1:r]
              store.Phi = A.tilde[(r+1):k,(r+1):k]
              store.Lambda = Lambda.tilde[,1:r]
              store.Sigma_u = Sigma.u.tilde[1:r,1:r]
              store.Sigma_epsilon = Sigma.u.tilde[(r+1):k,(r+1):k]
              store.factors.cov = store.cov[1:r,1:r,]
              store.errors.cov = store.cov[(r+1):k,(r+1):k,]
              
              dimnames(store.A) = list(factor.names,factor.names)
              dimnames(store.Phi) = list(series.names,series.names)
              dimnames(store.Lambda) = list(series.names, factor.names)
              dimnames(store.Sigma_u) = list(f.error.names, f.error.names)
              dimnames(store.Sigma_epsilon) = list(series.names, series.names)
              dimnames(store.factors.cov) = list(factor.names, factor.names)
              dimnames(store.errors.cov) = list(series.names, series.names)
              
              alpha.output[[alphas.index]] = list(params = list(a0_0 = a0_0,
                                                               P0_0 = P0_0,
                                                               A = store.A,
                                                               Phi = diag(store.Phi),
                                                               Lambda = store.Lambda,
                                                               Sigma_u = store.Sigma_u,
                                                               Sigma_epsilon = diag(store.Sigma_epsilon)),
                                                  state = list(factors = store.factors,
                                                               errors = store.errors,
                                                               factors.cov = store.factors.cov,
                                                               errors.cov = store.errors.cov))
              
                                                  
                                        
            }else{
              
              
              store.factors = store.state
              
              dimnames(store.factors) = list(obs.names, factor.names)
              
              store.A = A.tilde
              store.Lambda = Lambda.tilde
              store.Sigma_u = Sigma.u.tilde
              store.Sigma_epsilon = Sigma.eta
              store.factors.cov = store.cov
              
              dimnames(store.A) = list(factor.names,factor.names)
              dimnames(store.Lambda) = list(series.names, factor.names)
              dimnames(store.Sigma_u) = list(f.error.names, f.error.names)
              dimnames(store.Sigma_epsilon) = list(series.names, series.names)
              dimnames(store.factors.cov) = list(factor.names, factor.names)
              
              alpha.output[[alphas.index]] = list(params = list(a0_0 = a0_0,
                                                                P0_0 = P0_0,
                                                                A = store.A,
                                                                Lambda = store.Lambda,
                                                                Sigma_u = store.Sigma_u,
                                                                Sigma_epsilon = diag(store.Sigma_epsilon)),
                                                  state = list(factors = store.factors,
                                                               factors.cov = store.factors.cov))
              
              
            }
            
            
            
          }
          
          
          # store estimates if BIC improved 
          
            if(bic[alphas.index] < best.bic){
            
            best.EM = EM.fit
            best.KFS = KFS
            best.bic = bic[alphas.index]
            
          }
          
          
        } 
        
        
      ## Store the optimal outputs  
        time.emsparse = time.emsparse[1:length(bic)]
        num_iter = num_iter[1:length(bic)]
        alphas.used = alphas[1:length(bic)]
        best.alpha = alphas[which.min(bic)]
        loglik.store = best.EM$loglik.store 
        converged = best.EM$converged 

        a0_0 = best.EM$a0_0
        P0_0 = best.EM$P0_0
        A.tilde = best.EM$A.tilde
        Lambda.tilde = best.EM$Lambda.tilde
        Sigma.u.tilde = best.EM$Sigma.u.tilde
        Sigma.eta = best.EM$Sigma.eta
        
        state.EM = t(best.KFS$at_n)
        covariance.EM = best.KFS$Pt_n
        
        
      ## Fill in missing data in X
        
        fit_x = state.EM[,1:r] %*% t(Lambda.tilde[,1:r])
        if(standardize){
          fit_X = kronecker(t(X.sd),rep(1,n))*fit_x + kronecker(t(X.mean),rep(1,n))
        }else{
          fit_X = fit_x
        }
        
        dimnames(fit_x) = dimnames(X)
        dimnames(fit_X) = dimnames(X)
        
        if(store.parameters){
          names(alpha.output) = paste0('alpha=',round(alphas.used,4))
        }else{
          alpha.output = NULL
        }
        
        
      ## Output for EM-sparse - depends on if err = 'AR1' or 'IID'
        
        if(err == 'AR1'){
          
          factors = state.EM[,1:r]
          errors = state.EM[,(r+1):k]
          
          dimnames(factors) = list(obs.names, factor.names)
          dimnames(errors) = list(obs.names, series.names)
          
          A = A.tilde[1:r,1:r]
          Phi = A.tilde[(r+1):k,(r+1):k]
          Lambda = Lambda.tilde[,1:r]
          Sigma_u = Sigma.u.tilde[1:r,1:r]
          Sigma_epsilon = Sigma.u.tilde[(r+1):k,(r+1):k]
          factors.cov = covariance.EM[1:r,1:r,]
          errors.cov = covariance.EM[(r+1):k,(r+1):k,]
          
          dimnames(A) = list(factor.names,factor.names)
          dimnames(Phi) = list(series.names,series.names)
          dimnames(Lambda) = list(series.names, factor.names)
          dimnames(Sigma_u) = list(f.error.names, f.error.names)
          dimnames(Sigma_epsilon) = list(series.names, series.names)
          dimnames(factors.cov) = list(factor.names, factor.names)
          dimnames(errors.cov) = list(series.names, series.names)
          
          output = list(data = list(X = X.input,
                                    standardize = standardize,
                                    X.mean = X.mean, 
                                    X.sd = X.sd,
                                    X.bal = initialise$X.bal,
                                    eigen = initialise$eigen,
                                    fitted = fit_x,
                                    fitted.unscaled = fit_X,
                                    method = alg,
                                    err = err,
                                    call = match.call()),
                        params = list(a0_0 = a0_0,
                                      P0_0 = P0_0,
                                      A = A,
                                      Phi = diag(Phi),
                                      Lambda = Lambda,
                                      Sigma_u = Sigma_u,
                                      Sigma_epsilon = diag(Sigma_epsilon)),
                        state = list(factors = factors,
                                     errors = errors,
                                     factors.cov = factors.cov,
                                     errors.cov = errors.cov),
                        em = list(converged = converged,
                                  alpha_grid = alphas.used,
                                  alpha_opt = best.alpha,
                                  bic = bic,
                                  loglik = loglik.store,
                                  num_iter = num_iter,
                                  tol = threshold,
                                  max_iter = max_iter,
                                  em_time = time.emsparse),
                        alpha.output = alpha.output)

          
        }else {
          
          factors = state.EM
          
          dimnames(factors) = list(obs.names, factor.names)
          
          A = A.tilde
          Lambda = Lambda.tilde
          Sigma_u = Sigma.u.tilde
          Sigma_epsilon = Sigma.eta
          factors.cov = covariance.EM
          
          dimnames(A) = list(factor.names,factor.names)
          dimnames(Lambda) = list(series.names, factor.names)
          dimnames(Sigma_u) = list(f.error.names, f.error.names)
          dimnames(Sigma_epsilon) = list(series.names, series.names)
          dimnames(factors.cov) = list(factor.names, factor.names)
          
          output = list(data = list(X = X.input,
                                    standardize = standardize,
                                    X.mean = X.mean, 
                                    X.sd = X.sd,
                                    X.bal = initialise$X.bal,
                                    eigen = initialise$eigen,
                                    fitted = fit_x,
                                    fitted.unscaled = fit_X,
                                    method = alg,
                                    err = err,
                                    call = match.call()),
                        params = list(a0_0 = a0_0,
                                      P0_0 = P0_0,
                                      A = A,
                                      Lambda = Lambda,
                                      Sigma_u = Sigma_u,
                                      Sigma_epsilon = diag(Sigma_epsilon)),
                        state = list(factors = factors,
                                     factors.cov = factors.cov),
                        em = list(converged = converged,
                                  alpha_grid = alphas.used,
                                  alpha_opt = best.alpha,
                                  bic = bic,
                                  loglik = loglik.store,
                                  num_iter = num_iter,
                                  tol = threshold,
                                  max_iter = max_iter,
                                  em_time = time.emsparse),
                        alpha.output = alpha.output)
          
          
        }
        
        class(output) <- 'sparseDFM'
        return(output)
    
  }
  
}
  
