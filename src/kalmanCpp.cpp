#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Symmetrise the matrix
inline mat symmat(const mat& P) {
  return 0.5 * (P + P.t());
}

// [[Rcpp::export]]
List kalmanCpp(const mat& X, const rowvec& a0_0, const mat& P0_0, const mat& A,
               const mat& Lambda, const mat& Sig_e, const mat& Sig_u) {

  //  Kalman Filter and Smoother equations from Shumway and Stoffer (1982)
  
  //  Works for a general state space model of the form:
  //
  //  X_t = Lambda*F_t + e_t,  e_t ~ N(0,Sig_e)
  //  F_t = A*F_{t-1} + u_t,   u_t ~ N(0,Sig_u)

  //  Inputs:
  //
  //  X: n x p, matrix of (stationary) time series
  //  a0_0: 1 x k, initial state mean vector
  //  P0_0: k x k, initial state covariance matrix
  //  A: k x k, state matrix
  //  Lambda: p x k, measurement matrix
  //  Sig_e: p x p, measurement equation residuals covariance matrix (diagonal)
  //  Sig_u: k x k, state equation residuals covariance matrix
  //
  //  Outputs:
  //
  //  logl: log-likelihood required for convergence check in EM
  //  at_tlag: k x n, predicted state mean vectors
  //  Pt_tlag: k x k x n, predicted state covariance matrices
  //  at_t: k x n, filtered state mean vectors
  //  Pt_t: k x k x n, filtered state covariance matrices
  //  at_n: k x n, smoothed state mean vectors
  //  Pt_n: k x k x n, smoothed state covariance matrices
  //  Pt_tlag_n: k x k x n, smoothed state covariance with lag
  //
  //  NOTE: For the DFM with AR(1) errors we have that k = r + p where r is
  //        the number of factors and p is the number of variables.

  // Initialise
  const uword n = X.n_rows;
  const uword k = A.n_rows;

  mat at_tlag(k, n, fill::zeros);
  cube Pt_tlag(k, k, n, fill::zeros);
  mat at_t(k, n + 1, fill::zeros);
  cube Pt_t(k, k, n + 1, fill::zeros);
  mat at_n(k, n + 1, fill::zeros);
  cube Pt_n(k, k, n + 1, fill::zeros);
  cube Pt_tlag_n(k, k, n, fill::zeros);

  at_t.col(0) = a0_0.t(); // initial state mean at t=0
  Pt_t.slice(0) = P0_0;   // initial state covariance at t=0

  const mat y = X.t();    // work with p x n matrix y

  double logl = 0;        // log-likelihood required for convergence check in EM

  // Kalman filter loop for t = 1,...,n
  mat Lambda_t;
  const mat Ik = eye(k, k);
  mat KG;
  const double log2pi = log(2.0 * datum::pi);

  for (uword t = 0; t < n; ++t) {
    // Remove row of Lambda plus row and column of Sig_e if X is missing at t
    const mat y_t = y.col(t);
    const uvec na_omit_t = find_finite(y_t);
    Lambda_t = Lambda.rows(na_omit_t);
    const mat Sig_e_t = diagmat(Sig_e.submat(na_omit_t, na_omit_t));
    const mat inv_Sig_e_t = inv(Sig_e_t);
    const vec diag_Sig_e_t = diagvec(Sig_e_t);

    // Prediction equations
    at_tlag.col(t) = A * at_t.col(t);
    Pt_tlag.slice(t) = symmat(A * Pt_t.slice(t) * A.t() + Sig_u);

    // Update equations (omitting missing values at t)
    const mat Pt_tlag_Lambda = Pt_tlag.slice(t) * Lambda_t.t();
    const mat inov_res = y_t.elem(na_omit_t) - Lambda_t * at_tlag.col(t);
    const mat inov_cov = Lambda_t * Pt_tlag_Lambda + Sig_e_t;
    const mat GG = Ik + Pt_tlag_Lambda * inv_Sig_e_t * Lambda_t;
    const mat inov_cov_inv = inv_Sig_e_t - inv_Sig_e_t * Lambda_t
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
  const mat factors_KF = at_t.cols(1, n);
  const cube covariance_KF = Pt_t.slices(1, n);

  // Initialise Kalman smoother with t=n of filtered mean and covariance
  at_n.col(n) = at_t.col(n);
  Pt_n.slice(n) = Pt_t.slice(n);

  // Initialise Kalman smooth covariance with lag and smoother gain
  Pt_tlag_n.slice(n - 1) = (Ik - KG * Lambda_t) * A * Pt_t.slice(n - 1);
  mat J_2 = Pt_t.slice(n - 1) * A.t() * inv_sympd(Pt_tlag.slice(n - 1));

  // Kalman smoother loop for t = n,...,1
  for (uword t = n - 1; t != uword(-1); --t) {
    const mat J_1 = J_2;
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
  const mat factors_KS = at_n.cols(1, n);
  const cube covariance_KS = Pt_n.slices(1, n);

  // Output
  return List::create(Named("logl")       = logl,
                      Named("at_t")       = factors_KF,
                      Named("Pt_t")       = covariance_KF,
                      Named("at_n")       = factors_KS,
                      Named("Pt_n")       = covariance_KS,
                      Named("Pt_tlag_n")  = Pt_tlag_n);
}
