solveLambda <- function(D_cube, C_n, nu = 1, alpha, maxiter = 100, eps = 1e-4){

  ## New solver for Lambda M-step in Sparsified EM algorithm that handles missing data

  ## Inputs:
  ##
  ##  D_cube: r x r x N - solve(\sum_{t=1}^n [E(F_t F_t^T | \Omega_n) \otimes (W_t \bsig_e^{-1} W_t) + \nu*I])
  ##  C_n: N x r - \sum_{t=1}^n [W_t \bsig_e^{-1} W_t X_t E(F_t | \Omega_n)]
  ##  nu: scaling term for the augmentation in ADMM
  ##  alpha: LASSO regularisation parameter
  ##  maxiter: maximum number of iterations
  ##  eps: convergence threshold

  ## Outputs:
  ##
  ##  Lambda: primal variable
  ##  Z: auxiliary variable (should be sparse)
  ##  U: dual variable
  ##  opt: difference in estimates along optimization path
  ##  iter: number of iterations till convergence/end

  N = dim(C_n)[1]
  r = dim(C_n)[2]

  # initialise at 0
  Lambda = matrix(data = 0, nrow=N, ncol=r)
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
    Lambda = SparseDFM::fastLambda(D_cube, CC_n)

    # Shrinkage Operator - nu = 1
    #Z = softThresh(Lambda+U,alpha)
    Z = SparseDFM::softThreshMatrix(Lambda+U,alpha)

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
