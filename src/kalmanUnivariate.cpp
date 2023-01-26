#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Univariate filtering (sequential processing) for fast KFS
//'
//' @description
//' Univariate treatment (sequential processing) of the multivariate Kalman filter and smoother equations for fast implementation. Refer to Koopman and Durbin (2000).
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
//' For full details of the univariate filtering approach, please refer to Mosley et al. (2023). Note that \eqn{n}{n} is the number of observations, \eqn{p}{p} is the number of time series, and \eqn{k}{k} is the number of states.
//'
//' @return logl log-likelihood of the innovations from the Kalman filter 
//' @return at_t \eqn{k \times n}{k x n}, filtered state mean vectors
//' @return Pt_t \eqn{k \times k \times n}{k x k x n}, filtered state covariance matrices
//' @return at_n \eqn{k \times n}{k x n}, smoothed state mean vectors
//' @return Pt_n \eqn{k \times k \times n}{k x k x n}, smoothed state covariance matrices
//' @return Pt_tlag_n \eqn{k \times k \times n}{k x k x n}, smoothed state covariance with lag
//'
//' @import Rcpp
//'
//' @references 
//' Koopman, S. J., & Durbin, J. (2000). Fast filtering and smoothing for multivariate state space models. \emph{Journal of Time Series Analysis, 21}(3), 281-296.
//'
//' Mosley, L., Chan, TS., & Gibberd, A. (2023). sparseDFM: An R Package to Estimate Dynamic Factor Models with Sparse Loadings.
//'
//' @export
// [[Rcpp::export]]
List kalmanUnivariate(const arma::mat& X, const arma::mat& a0_0, const arma::mat& P0_0,
                      const arma::mat& A, const arma::mat& Lambda, const arma::mat& Sig_e,
                      const arma::mat& Sig_u) {

  // Initialise
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  const arma::uword k = A.n_rows;

  arma::mat at_tlag(k, n, arma::fill::zeros);        // predicted state mean
  arma::cube Pt_tlag(k, k, n, arma::fill::zeros);    // predicted state covariance
  arma::mat at_t(k, n, arma::fill::zeros);           // filtered state mean
  arma::cube Pt_t(k, k, n, arma::fill::zeros);       // filtered state covariance
  arma::mat at_n(k, n, arma::fill::zeros);           // smoothed state mean
  arma::cube Pt_n(k, k, n, arma::fill::zeros);       // smoothed state covariance
  arma::cube Pt_tlag_n(k, k, n, arma::fill::zeros);  // smoothed state covariance with lag

  arma::vec at_i = vectorise(a0_0);           // initial state mean
  arma::mat Pt_i = P0_0;                      // initial state covariance

  double logl = 0.0;                    // log-likelihood

  arma::mat vt(p, n, arma::fill::zeros);             // innovation
  arma::mat inv_Ft(p, n, arma::fill::zeros);        // inverse innovation variance
  arma::cube Kt(k, p, n, arma::fill::zeros);         // Kalman gain

  // Kalman filter loop for t = 1,...,n
  const arma::vec diag_Sig_e = diagvec(Sig_e);
  const double log2pi = log(2.0 * arma::datum::pi);

  for (arma::uword t = 0; t < n; ++t) {
    // Prediction equations
    at_i = A * at_i;
    Pt_i = A * Pt_i * A.t() + Sig_u;
    at_tlag.col(t) = at_i;
    Pt_tlag.slice(t) = Pt_i;

    // Update equations
    for (arma::uword i = 0; i < p; ++i) {
      const double yt_i = X(t, i);      // y = t(X)
      if (std::isnan(yt_i)) continue;   // omit NA

      const arma::rowvec Zt_i = Lambda.row(i);
      const double vt_i = yt_i - dot(Zt_i, at_i);
      const arma::vec PZt_i = Pt_i * Zt_i.t();
      const double Ft_i = dot(Zt_i, PZt_i) + diag_Sig_e(i);

      // Skip yt_i if Ft_i <= 0
      if (Ft_i <= 0) continue;

      const double inv_Ft_i = 1.0 / Ft_i;
      const arma::vec Kt_i = PZt_i * inv_Ft_i;
      at_i += Kt_i * vt_i;
      Pt_i -= Kt_i * PZt_i.t();

      // Enforce symmetry
      Pt_i = 0.5 * (Pt_i + Pt_i.t());

      // Calculate log-likelihood
      logl -= 0.5 * (log2pi + log(Ft_i) + vt_i * inv_Ft_i * vt_i);

      vt(i, t) = vt_i;
      inv_Ft(i, t) = inv_Ft_i;
      Kt.slice(t).col(i) = Kt_i;
    }
    at_t.col(t) = at_i;
    Pt_t.slice(t) = Pt_i;
  }

  // Initialise Kalman smoother with zeros
  arma::vec rt_i(k, arma::fill::zeros);
  arma::mat Nt_i(k, k, arma::fill::zeros);

  // Kalman smoother loop for t = n,...,1
  const arma::mat Ik = arma::eye(k, k);
  for (arma::uword t = n - 1; t != arma::uword(-1); --t) {
    for (arma::uword i = p - 1; i != arma::uword(-1); --i) {
      const arma::rowvec Zt_i = Lambda.row(i);
      const double inv_Ft_i = inv_Ft(i, t);
      const arma::vec Kt_i = Kt.slice(t).col(i);
      const arma::mat Lt_i = Ik - Kt_i * Zt_i;
      if (inv_Ft_i) {
        // yt_i is not missing and Ft_i is not zero
        const arma::vec ZFt_i = Zt_i.t() * inv_Ft_i;
        rt_i = ZFt_i * vt(i, t) + Lt_i.t() * rt_i;
        Nt_i = ZFt_i * Zt_i + Lt_i.t() * Nt_i * Lt_i;
      } else {
        rt_i = Lt_i.t() * rt_i;
        Nt_i = Lt_i.t() * Nt_i * Lt_i;
      }
    }
    const arma::mat Pt_1 = Pt_tlag.slice(t);
    at_n.col(t) = at_tlag.col(t) + Pt_1 * rt_i;
    Pt_n.slice(t) = Pt_1 - Pt_1 * Nt_i * Pt_1;
    if (t > 0) {
      rt_i = A.t() * rt_i;
      Nt_i = A.t() * Nt_i * A;
    }
  }

  // Calculate smoothed state covariance with lag
  Pt_tlag_n.slice(0) = Pt_n.slice(0) * solve(Pt_tlag.slice(0), A) * P0_0;
  for (arma::uword t = 1; t < n; ++t) {
    Pt_tlag_n.slice(t) = Pt_n.slice(t) * solve(Pt_tlag.slice(t), A) * Pt_t.slice(t - 1);
  }

  // Output
  return List::create(Named("logl")       = logl,
                      Named("at_t")       = at_t,
                      Named("Pt_t")       = Pt_t,
                      Named("at_n")       = at_n,
                      Named("Pt_n")       = Pt_n,
                      Named("Pt_tlag_n")  = Pt_tlag_n);
}
