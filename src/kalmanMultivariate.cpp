#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Symmetrise the matrix
inline arma::mat symmat(const arma::mat& P) {
  return 0.5 * (P + P.t());
}

//' Classic Multivariate KFS Equations
//'
//' @description 
//' Implementation of the classic multivariate Kalman filter and smoother equations of Shumway and Stoffer (1982).
//'
//' @param X n x p, numeric matrix of (stationary) time series 
//' @param a0_0 k x 1, initial state mean vector 
//' @param P0_0 k x k, initial state covariance matrix
//' @param A k x k, state transition matrix
//' @param Lambda p x k, measurement matrix 
//' @param Sig_e p x p, measurement equation residuals covariance matrix (diagonal)
//' @param Sig_u k x k, state equation residuals covariance matrix
//'
//' @details 
//' For full details of the classic multivariate KFS approach, please refer to Mosley et al. (2023). Note that \eqn{n}{n} is the number of observations, \eqn{p}{p} is the number of time series, and \eqn{k}{k} is the number of states.
//'
//' @return logl log-likelihood of the innovations from the Kalman filter 
//' @return at_t \eqn{k \times n}{k x n}, filtered state mean vectors
//' @return Pt_t \eqn{k \times k \times n}{k x k x n}, filtered state covariance matrices
//' @return at_n \eqn{k \times n}{k x n}, smoothed state mean vectors
//' @return Pt_n \eqn{k \times k \times n}{k x k x n}, smoothed state covariance matrices
//' @return Pt_tlag_n \eqn{k \times k \times n}{k x k x n}, smoothed state covariance with lag
//'
//' @references 
//' Mosley, L., Chan, TS., & Gibberd, A. (2023). sparseDFM: An R Package to Estimate Dynamic Factor Models with Sparse Loadings.
//' 
//' Shumway, R. H., & Stoffer, D. S. (1982). An approach to time series smoothing and forecasting using the EM algorithm. \emph{Journal of time series analysis, 3}(4), 253-264.
//'
//' @import Rcpp
//' @export
// [[Rcpp::export]]
List kalmanMultivariate(const arma::mat& X, const arma::mat& a0_0, const arma::mat& P0_0, const arma::mat& A,
						const arma::mat& Lambda, const arma::mat& Sig_e, const arma::mat& Sig_u) {

  // Initialise
  const arma::uword n = X.n_rows;
  const arma::uword k = A.n_rows;

  arma::mat at_tlag(k, n, arma::fill::zeros);
  arma::cube Pt_tlag(k, k, n, arma::fill::zeros);
  arma::mat at_t(k, n + 1, arma::fill::zeros);
  arma::cube Pt_t(k, k, n + 1, arma::fill::zeros);
  arma::mat at_n(k, n + 1, arma::fill::zeros);
  arma::cube Pt_n(k, k, n + 1, arma::fill::zeros);
  arma::cube Pt_tlag_n(k, k, n, arma::fill::zeros);

  at_t.col(0) = vectorise(a0_0);  // initial state mean at t=0
  Pt_t.slice(0) = P0_0;           // initial state covariance at t=0

  const arma::mat y = X.t();      // work with p x n matrix y

  double logl = 0;        // log-likelihood required for convergence check in EM

  // Kalman filter loop for t = 1,...,n
  arma::mat Lambda_t;
  const arma::mat Ik = arma::eye(k, k);
  arma::mat KG;
  const double log2pi = log(2.0 * arma::datum::pi);

  for (arma::uword t = 0; t < n; ++t) {
    // Remove row of Lambda plus row and column of Sig_e if X is missing at t
    const arma::mat y_t = y.col(t);
    const arma::uvec na_omit_t = find_finite(y_t);
    Lambda_t = Lambda.rows(na_omit_t);
    const arma::mat Sig_e_t = diagmat(Sig_e.submat(na_omit_t, na_omit_t));
    const arma::mat inv_Sig_e_t = inv(Sig_e_t);
    const arma::vec diag_Sig_e_t = diagvec(Sig_e_t);

    // Prediction equations
    at_tlag.col(t) = A * at_t.col(t);
    Pt_tlag.slice(t) = symmat(A * Pt_t.slice(t) * A.t() + Sig_u);

    // Update equations (omitting missing values at t)
    const arma::mat Pt_tlag_Lambda = Pt_tlag.slice(t) * Lambda_t.t();
    const arma::mat inov_res = y_t.elem(na_omit_t) - Lambda_t * at_tlag.col(t);
    const arma::mat inov_cov = Lambda_t * Pt_tlag_Lambda + Sig_e_t;
    const arma::mat GG = Ik + Pt_tlag_Lambda * inv_Sig_e_t * Lambda_t;
    const arma::mat inov_cov_inv = inv_Sig_e_t - inv_Sig_e_t * Lambda_t
        * solve(GG, Pt_tlag_Lambda) * inv_Sig_e_t;        // see Harvey (1990)
    KG = Pt_tlag_Lambda * inov_cov_inv;                   // Kalman gain
    at_t.col(t + 1) = at_tlag.col(t) + KG * inov_res;
    Pt_t.slice(t + 1) = symmat(Pt_tlag.slice(t) - KG * Lambda_t * Pt_tlag.slice(t));

    // Calculate log-likelihood
    if (all(diag_Sig_e_t > 0)) {
      const double detGG = det(GG);
      if (detGG > 0) {
        logl -= 0.5 * (inov_res.n_elem * log2pi + sum(log(diag_Sig_e_t))
            + log(detGG) + dot(inov_res, inov_cov_inv * inov_res));
      }
    }
  }

  // Output for Kalman filter
  const arma::mat factors_KF = at_t.cols(1, n);
  const arma::cube covariance_KF = Pt_t.slices(1, n);

  // Initialise Kalman smoother with t=n of filtered mean and covariance
  at_n.col(n) = at_t.col(n);
  Pt_n.slice(n) = Pt_t.slice(n);

  // Initialise Kalman smooth covariance with lag and smoother gain
  Pt_tlag_n.slice(n - 1) = (Ik - KG * Lambda_t) * A * Pt_t.slice(n - 1);
  arma::mat J_2 = Pt_t.slice(n - 1) * A.t() * inv_sympd(Pt_tlag.slice(n - 1));

  // Kalman smoother loop for t = n,...,1
  for (arma::uword t = n - 1; t != arma::uword(-1); --t) {
    const arma::mat J_1 = J_2;
    at_n.col(t) = at_t.col(t) + J_1 * (at_n.col(t + 1) - at_tlag.col(t));
    Pt_n.slice(t) = Pt_t.slice(t) + J_1 * (Pt_n.slice(t + 1) - Pt_tlag.slice(t))
        * J_1.t();
    if (t > 0) {
      J_2 = Pt_t.slice(t - 1) * A.t() * inv_sympd(Pt_tlag.slice(t - 1));
      Pt_tlag_n.slice(t - 1) = Pt_t.slice(t) * J_2.t() + J_1
          * (Pt_tlag_n.slice(t) - A * Pt_t.slice(t)) * J_2.t();
    }
  }

  // Output for Kalman smoother
  const arma::mat factors_KS = at_n.cols(1, n);
  const arma::cube covariance_KS = Pt_n.slices(1, n);

  // Output
  return List::create(Named("logl")       = logl,
                      Named("at_t")       = factors_KF,
                      Named("Pt_t")       = covariance_KF,
                      Named("at_n")       = factors_KS,
                      Named("Pt_n")       = covariance_KS,
                      Named("Pt_tlag_n")  = Pt_tlag_n);
}
