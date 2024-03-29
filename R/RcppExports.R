# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @import Rcpp
fastLambda <- function(D, C) {
    .Call(`_sparseDFM_fastLambda`, D, C)
}

#' Classic Multivariate KFS Equations
#'
#' @description 
#' Implementation of the classic multivariate Kalman filter and smoother equations of Shumway and Stoffer (1982).
#'
#' @param X n x p, numeric matrix of (stationary) time series 
#' @param a0_0 k x 1, initial state mean vector 
#' @param P0_0 k x k, initial state covariance matrix
#' @param A k x k, state transition matrix
#' @param Lambda p x k, measurement matrix 
#' @param Sig_e p x p, measurement equation residuals covariance matrix (diagonal)
#' @param Sig_u k x k, state equation residuals covariance matrix
#'
#' @details 
#' For full details of the classic multivariate KFS approach, please refer to Mosley et al. (2023). Note that \eqn{n}{n} is the number of observations, \eqn{p}{p} is the number of time series, and \eqn{k}{k} is the number of states.
#'
#' @return logl log-likelihood of the innovations from the Kalman filter 
#' @return at_t \eqn{k \times n}{k x n}, filtered state mean vectors
#' @return Pt_t \eqn{k \times k \times n}{k x k x n}, filtered state covariance matrices
#' @return at_n \eqn{k \times n}{k x n}, smoothed state mean vectors
#' @return Pt_n \eqn{k \times k \times n}{k x k x n}, smoothed state covariance matrices
#' @return Pt_tlag_n \eqn{k \times k \times n}{k x k x n}, smoothed state covariance with lag
#'
#' @references 
#' Mosley, L., Chan, TS., & Gibberd, A. (2023). sparseDFM: An R Package to Estimate Dynamic Factor Models with Sparse Loadings.
#' 
#' Shumway, R. H., & Stoffer, D. S. (1982). An approach to time series smoothing and forecasting using the EM algorithm. \emph{Journal of time series analysis, 3}(4), 253-264.
#'
#' @import Rcpp
#' @export
kalmanMultivariate <- function(X, a0_0, P0_0, A, Lambda, Sig_e, Sig_u) {
    .Call(`_sparseDFM_kalmanMultivariate`, X, a0_0, P0_0, A, Lambda, Sig_e, Sig_u)
}

#' Univariate filtering (sequential processing) for fast KFS
#'
#' @description
#' Univariate treatment (sequential processing) of the multivariate Kalman filter and smoother equations for fast implementation. Refer to Koopman and Durbin (2000).
#'
#' @param X n x p, numeric matrix of (stationary) time series 
#' @param a0_0 k x 1, initial state mean vector 
#' @param P0_0 k x k, initial state covariance matrix
#' @param A k x k, state transition matrix
#' @param Lambda p x k, measurement matrix 
#' @param Sig_e p x p, measurement equation residuals covariance matrix (diagonal)
#' @param Sig_u k x k, state equation residuals covariance matrix
#' 
#' @details 
#' For full details of the univariate filtering approach, please refer to Mosley et al. (2023). Note that \eqn{n}{n} is the number of observations, \eqn{p}{p} is the number of time series, and \eqn{k}{k} is the number of states.
#'
#' @return logl log-likelihood of the innovations from the Kalman filter 
#' @return at_t \eqn{k \times n}{k x n}, filtered state mean vectors
#' @return Pt_t \eqn{k \times k \times n}{k x k x n}, filtered state covariance matrices
#' @return at_n \eqn{k \times n}{k x n}, smoothed state mean vectors
#' @return Pt_n \eqn{k \times k \times n}{k x k x n}, smoothed state covariance matrices
#' @return Pt_tlag_n \eqn{k \times k \times n}{k x k x n}, smoothed state covariance with lag
#'
#' @import Rcpp
#'
#' @references 
#' Koopman, S. J., & Durbin, J. (2000). Fast filtering and smoothing for multivariate state space models. \emph{Journal of Time Series Analysis, 21}(3), 281-296.
#'
#' Mosley, L., Chan, TS., & Gibberd, A. (2023). sparseDFM: An R Package to Estimate Dynamic Factor Models with Sparse Loadings.
#'
#' @export
kalmanUnivariate <- function(X, a0_0, P0_0, A, Lambda, Sig_e, Sig_u) {
    .Call(`_sparseDFM_kalmanUnivariate`, X, a0_0, P0_0, A, Lambda, Sig_e, Sig_u)
}

#' @import Rcpp
softThreshScalar <- function(X, thresh) {
    .Call(`_sparseDFM_softThreshScalar`, X, thresh)
}

#' @import Rcpp
softThreshMatrix <- function(X, thresh) {
    .Call(`_sparseDFM_softThreshMatrix`, X, thresh)
}

#' @import Rcpp
solveCube <- function(A, B, nu = 0.0) {
    .Call(`_sparseDFM_solveCube`, A, B, nu)
}

