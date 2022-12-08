#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
cube solveCube(const cube& A, const mat& B, const double nu = 0.0) {
  // Validate inputs
  const uword r = A.n_rows;
  const uword n = A.n_slices;
  const uword p = B.n_rows;
  if (A.n_cols != r) stop("solveCube(): A must be r-by-r-by-n");
  if (B.n_cols != n) stop("solveCube(): A and B must have the same n");

  // Invert each block and stack them into a cube
  cube D(r, r, p, fill::none);
  for (uword i = 0; i < p; ++i) {
    mat block = nu * eye(r, r);
    for (uword t = 0; t < n; ++t) {
      block += B(i, t) * A.slice(t);
    }
    D.slice(i) = inv(block);
  }
  return D;
}
