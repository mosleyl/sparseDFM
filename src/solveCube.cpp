#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @import Rcpp
// [[Rcpp::export]]
arma::cube solveCube(const arma::cube& A, const arma::mat& B, const double nu = 0.0) {
  // Validate inputs
  const arma::uword r = A.n_rows;
  const arma::uword n = A.n_slices;
  const arma::uword p = B.n_rows;
  if (A.n_cols != r) stop("solveCube(): A must be r-by-r-by-n");
  if (B.n_cols != n) stop("solveCube(): A and B must have the same n");

  // Invert each block and stack them into a cube
  arma::cube D(r, r, p, arma::fill::none);
  for (arma::uword i = 0; i < p; ++i) {
    arma::mat block = nu * arma::eye(r, r);
    for (arma::uword t = 0; t < n; ++t) {
      block += B(i, t) * A.slice(t);
    }
    D.slice(i) = inv(block);
  }
  return D;
}
