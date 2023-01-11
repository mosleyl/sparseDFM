postProcess <- function(lambda_estimate, lambda_true){
  # Evaluation 
  #n = nrow(fact_est)
  N = nrow(lambda_estimate)
  r = ncol(lambda_estimate)
  # 1. correct scale 
  
      Lam_est = lambda_estimate * (norm(lambda_true[,1],'2'))
      #Fact_est = fact_est * (norm(lambda_true[,1],'2'))
  
  # 2. column permutations
  dist = list()
  lam_i = c()
  for(i in 1:r){
    for(j in 1:r){
      lam_i[j] = norm(abs(lambda_true[,j])-abs(Lam_est[,i]), type = '2')
    }
    dist[[i]] = lam_i
  }
  dist_unlist = unlist(dist)
  dist_matrix = matrix(dist_unlist, ncol = r, byrow = TRUE)
  
  new_dist_matrix = dist_matrix
  lam_hat_new = matrix(NA, nrow = N, ncol = r)
  #fact_est_new = matrix(NA, nrow = n, ncol = r)
  
  while(!all(new_dist_matrix == 1e23)) {
    
    min_idx = which(new_dist_matrix == min(new_dist_matrix), arr.ind = TRUE)
    lam_hat_new[,min_idx[2]] = Lam_est[,min_idx[1]]
    new_dist_matrix[min_idx[1],] = 1e23
    new_dist_matrix[,min_idx[2]] = 1e23
    
    #fact_est_new[,min_idx[2]] = fact_est[,min_idx[1]]
    
  }
  
  # 3. make same sign as original Lambda 
  Lam_est_final = abs(lam_hat_new)
  return(Lam_est_final)
  #return(list('lam' = Lam_est_final, 'fact' = fact_est_new))
  
}
