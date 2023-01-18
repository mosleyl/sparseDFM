#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @import Rcpp
// [[Rcpp::export]]
arma::mat fastLambda(const arma::cube& D, const arma::mat& C) {
  // Validate inputs
  const arma::uword p = C.n_rows;
  const arma::uword r = C.n_cols;
  if (D.n_rows != r || D.n_cols != r || D.n_slices != p)
    stop("fastLambda(): D must be r-by-r-by-p");

  // Multiply by vec(C) then unvec
  arma::mat Lambda(p, r, arma::fill::zeros);
  for (arma::uword i = 0; i < r; ++i) {
    for (arma::uword j = 0; j < r; ++j) {
      arma::vec Dij = D.tube(i, j);
      Lambda.col(i) += Dij % C.col(j);
    }
  }
  return Lambda;
}
