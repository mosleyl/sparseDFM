#' Tune for the number of factors to use 
#' 
#' @description 
#' Uses Bai and Ng (2002) information criteria approach. Missing data is interpolated using the \code{fillNA} function. 
#' 
#' @param X \code{n x p} numeric data matrix or data frame of (stationary) time series.
#' @param type Character. Option for which information criteria to use. Default is 2. 
#' @param standardize Logical. Standardize the data before estimating the model. Default is \code{TRUE}.
#' @param r.max Integer. Maximum number of factors to search for. Default is min(15,ncol(X)-1). 
#' @param plot Logical. Make a plot showing the IC value for each of the number of factors considered. Default is \code{TRUE}.
#' 
#' @returns 
#' The number of factors to use according to Bai and Ng (2002) information criteria.
#'
#' @details 
#' To calculate the number of factors to use in the model, the information criteria approach of Bai and Ng (2002) is used. This can be done before \code{sparseDFM} is fitted to the data to determine \code{r}. Bai and Ng (2002) consider 3 types of information criteria with different penalties of the form:
#'
#' \deqn{IC_1(r) = log\left(V_r(\hat{\bm{F}},\hat{\bm{\Lambda}})\right) + r \left( \frac{n+p}{np}\right)log\left( \frac{np}{n+p}\right)}{IC_1(r) = log(V_r(\hat{\bm{F}},\hat{\Lambda})) + r ((n+p)/(np))log((np)/(n+p))}
#' \deqn{IC_2(r) = log\left(V_r(\hat{\bm{F}},\hat{\bm{\Lambda}})\right) + r \left( \frac{n+p}{np} \right)log\left( min\{n,p\}\right)}{IC_2(r) = log(V_r(\hat{F},\hat{\Lambda})) + r ((n+p)/(np))log(min(n,p))}
#' \deqn{IC_3(r) = log\left(V_r(\hat{\bm{F}},\hat{\bm{\Lambda}})\right) + r \frac{log\left( min\{n,p\}\right)}{min\{n,p\}}}{IC_3(r) = log(V_r(\hat{F},\hat{\Lambda})) + r log (min(n,p)))/(min(n,p))}
#' 
#' The sum of squared residuals for \eqn{r}{r} factors \eqn{V_r(\hat{\bm{F}},\hat{\bm{\Lambda}}) = \sum_{i=1}^p\sum_{t=1}^n E[\hat{\epsilon}_{i,t}^2]/np}{V_r(\hat{F},\hat{\Lambda}) = \sum_{i=1}^p\sum_{t=1}^n E[\hat{\epsilon}_{i,t}^2]/np} with \eqn{\hat{\epsilon}_{i,t} = X_{t,i}-\hat{\bm{F}}_t\hat{\bm{\Lambda}}_i}{\hat{\epsilon}_{i,t} = X_{t,i}-\hat{F}_t\hat{\Lambda}_i} is found using PCA on the standardized data set \eqn{\bm{X}}{X}. The estimated factors \eqn{\hat{\bm{F}}}{\hat{F}} corresponding to the principle components and the estimated loadings \eqn{\hat{\bm{\Lambda}}}{\hat{\Lambda}} corresponding to the eigenvectors. Should the data contain missing values, then the missing data is interpolated using \code{fillNA}.
#'
#' The number of factors to use will correspond to \eqn{\argmin_r IC_i(r)}{argmin(r) IC(i)(r)} for \eqn{i=1,2}{i=1,2} or \eqn{3}{3}. Type 2 is the highest when working in finite samples and therefore is set to default. 
#'
#' @references 
#' Bai, J., & Ng, S. (2002). Determining the number of factors in approximate factor models. \emph{Econometrica, 70}(1), 191-221.
#'
#' @importFrom stats cov 
#' @importFrom graphics plot points 
#' 
#' @export


tuneFactors <- function(X, type = 2, standardize = TRUE, r.max = min(15,ncol(X)-1), plot = TRUE){
  
  X = as.matrix(X)
  
  if(standardize){
    X = scale(X)
  }
  
  if(anyNA(X)){
    message("Data contains missing values: imputing data with fillNA()")
    X = fillNA(X)
    X = X$X
  }
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  evd = eigen(cov(X))
  loadings = evd$vectors 
  factors = X %*% loadings 
  
  IC = c()
  
  for(r in 1:r.max){
    
    eps = X - tcrossprod(factors[,1:r], loadings[,1:r])
    logV = log(sum(colSums(eps^2)/(n*p)))
    
    if(type == 1){
      IC[r] = logV + r*((n+p)/(n*p))*log((n*p)/(n+p))
    }else if(type == 2){
      IC[r] = logV + r*((n+p)/(n*p))*log(min(n,p))
    }else{
      IC[r] = logV + r*log(min(n,p))/(min(n,p))
    }
    
  }
  
  if(plot){
    
    par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    par(mfrow = c(2,1))
    
    plot(IC, main = paste("IC",type), xlab = "Number of factors", ylab = "Index")
    points(which.min(IC), IC[which.min(IC)], pch = 16, col ="red")
  
    ev = evd$values 
    pve = (ev/sum(ev))*100
    if(length(ev) > r.max){
      ev = ev[1:r.max]
      pve = pve[1:r.max]
    }
    plot(pve, type = 'o', ylab = '% Variance explained',xlab = "Number of factors")
    grid()
    points(which.min(IC), pve[which.min(IC)], pch = 16, col='red')
    
    par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    par(mfrow = c(1,1))
    
  }
  
  return(paste("The chosen number of factors using criteria type ", type, " is ", which.min(IC)))
  
  
}
