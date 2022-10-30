#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
mat softThreshScalar(const mat& X, const double thresh) {
  // Takes matrix X and soft-thresholds entries

  return sign(X) % clamp(abs(X) - thresh, 0, datum::inf);
}

// [[Rcpp::export]]
mat softThreshMatrix(const mat& X, const mat& thresh) {
  // Same as above, but thresh is now a matrix

  return sign(X) % clamp(abs(X) - thresh, 0, datum::inf);
}
