#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @import Rcpp
// [[Rcpp::export]]
arma::mat softThreshScalar(const arma::mat& X, const double thresh) {
  // Takes matrix X and soft-thresholds entries

  return sign(X) % clamp(abs(X) - thresh, 0, arma::datum::inf);
}

//' @import Rcpp
// [[Rcpp::export]]
arma::mat softThreshMatrix(const arma::mat& X, const arma::mat& thresh) {
  // Same as above, but thresh is now a matrix

  return sign(X) % clamp(abs(X) - thresh, 0, arma::datum::inf);
}
