#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat fastLambda(const cube& D, const mat& C) {
  // Validate inputs
  const uword p = C.n_rows;
  const uword r = C.n_cols;
  if (D.n_rows != r || D.n_cols != r || D.n_slices != p)
    stop("fastLambda(): D must be r-by-r-by-p");

  // Multiply by vec(C) then unvec
  mat Lambda(p, r, fill::zeros);
  for (uword i = 0; i < r; ++i) {
    for (uword j = 0; j < r; ++j) {
      vec Dij = D.tube(i, j);
      Lambda.col(i) += Dij % C.col(j);
    }
  }
  return Lambda;
}
