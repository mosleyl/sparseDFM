#' New solver for Lambda M-step in Sparsified EM algorithm that handles missing data 
#'  
#' @param D_cube r x r x p 
#' @param C_n p x r 
#' @param nu scaling term for the augmentation in ADMM 
#' @param alpha LASSO regularisation parameter 
#' @noRd

## Outputs
##
##  Lambda: primal variable 
##  Z: auxiliary variable (should be sparse)
##  U: dual variable 
##  opt: difference in estimates along optimization path 
##  iter: number of iterations till convergence/end



solveLambda <- function(D_cube, C_n, nu = 1, alpha){
  
  maxiter = 100
  eps = 1e-4
  
  p = dim(C_n)[1]
  r = dim(C_n)[2]
  
  # initialise at 0 
  Lambda = matrix(data = 0, nrow=p, ncol=r)
  Z = Lambda
  U = Lambda
  
  iter = 1
  converged = FALSE
  opt = 0; opt2=0
  
  
  # Main ADMM Loop
  while (converged == FALSE && iter < maxiter){
    
    Lambdaold = Lambda
    
    # Solve for primal variable
    CC_n = C_n + nu*(Z-U)
    Lambda = fastLambda(D_cube, CC_n)
    
    # Shrinkage Operator - nu = 1
    Z = softThreshMatrix(Lambda+U,alpha)
    
    # Dual variable 
    U = U + Lambda - Z
    
    # Track convergence
    opt[iter] = norm(Lambda-Lambdaold,'F')
    if(opt[iter] < eps){
      converged = TRUE
    }else{
      iter = iter+1
    }
  }
  
  return(list('Lambda'=Lambda,'Z'=Z,'U'=U,'iter'=iter,'opt'=opt))
} 
