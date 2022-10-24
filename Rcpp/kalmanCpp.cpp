#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List kalmanCpp(const mat& X, const vec& x_init, const mat& V_init, const mat& A,
               const mat& Lambda, const mat& Sig_e, const mat& Sig_u) {
  // Function for the kalman filter and smoother
  // Initial parameters obtained from PCA and OLS

  const mat y = X.t();
  const uword n = y.n_cols;
  const uword r = A.n_rows;
  mat xt_t(r, n + 1, fill::zeros);
  cube Vt_t(r, r, n + 1, fill::zeros);
  mat xt_tlag(r, n, fill::zeros);
  cube Vt_tlag(r, r, n, fill::zeros);

  // Kalman filter
  xt_t.col(0) = x_init;
  Vt_t.slice(0) = V_init;
  double logl = 0;
  mat Lambda_j;
  mat KG;

  for (uword j = 0; j < n; ++j) {
    // Missing data
    const mat y_j = y.col(j);
    const uvec na_omit_j = find_finite(y_j);
    Lambda_j = Lambda.rows(na_omit_j);
    const mat Sig_e_j = Sig_e.submat(na_omit_j, na_omit_j);
    const mat inv_Sig_e_j = inv(diagmat(Sig_e_j));

    // Prediction equations
    xt_tlag.col(j) = A * xt_t.col(j);
    Vt_tlag.slice(j) = A * Vt_t.slice(j) * A.t() + Sig_u;

    // Update equations
    const mat inov_res = y_j.elem(na_omit_j) - Lambda_j * xt_tlag.col(j);
    const mat inov_cov = Lambda_j * Vt_tlag.slice(j) * Lambda_j.t() + Sig_e_j;
    const mat GG = eye(r, r) + Vt_tlag.slice(j) * Lambda_j.t() * inv_Sig_e_j
        * Lambda_j;
    const mat inov_cov_inv = inv_Sig_e_j - inv_Sig_e_j * Lambda_j * pinv(GG)
        * Vt_tlag.slice(j) * Lambda_j.t() * inv_Sig_e_j;
    KG = Vt_tlag.slice(j) * Lambda_j.t() * inov_cov_inv;
    xt_t.col(j + 1) = xt_tlag.col(j) + KG * inov_res;
    Vt_t.slice(j + 1) = Vt_tlag.slice(j) - KG * Lambda_j * Vt_tlag.slice(j);

    // Log-likelihood
    const double inov_cov_det = prod(diagvec(Sig_e_j)) * det(GG);
    logl -= 0.5 * (inov_res.n_elem * log(2 * datum::pi) + log(abs(inov_cov_det))
        + as_scalar(inov_res.t() * inov_cov_inv * inov_res));
  }

  // Rauch-Tung-Striebel smoother
  mat xt_n(r, n + 1, fill::zeros);
  cube Vt_n(r, r, n + 1, fill::zeros);
  cube Vt_tlag_n(r, r, n, fill::zeros);

  xt_n.col(n) = xt_t.col(n);
  Vt_n.slice(n) = Vt_t.slice(n);
  Vt_tlag_n.slice(n - 1) = (eye(r, r) - KG * Lambda_j) * A * Vt_t.slice(n - 1);
  mat J_2 = Vt_t.slice(n - 1) * A.t() * pinv(Vt_tlag.slice(n - 1));

  for (uword j = n - 1; j != uword(-1); --j) {
    const mat J_1 = J_2;
    xt_n.col(j) = xt_t.col(j) + J_1 * (xt_n.col(j + 1) - xt_tlag.col(j));
    Vt_n.slice(j) = Vt_t.slice(j) + J_1 * (Vt_n.slice(j + 1) - Vt_tlag.slice(j))
        * J_1.t();
    if (j > 0) {
      J_2 = Vt_t.slice(j - 1) * A.t() * pinv(Vt_tlag.slice(j - 1));
      Vt_tlag_n.slice(j - 1) = Vt_t.slice(j) * J_2.t() + J_1
          * (Vt_tlag_n.slice(j) - A * Vt_t.slice(j)) * J_2.t();
    }
  }

  return List::create(Named("xt_n")           = xt_n,
                      Named("Vt_n")           = Vt_n,
                      Named("Vt_tlag_n")      = Vt_tlag_n,
                      Named("logl")           = logl,
                      Named("factors.KF")     = xt_t.t(),
                      Named("covariance.KF")  = Vt_t,
                      Named("factors.KS")     = xt_n.t(),
                      Named("covariance.KS")  = Vt_n);
}
