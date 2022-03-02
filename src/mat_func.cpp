// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Eigen/Eigen>
using namespace Rcpp;



// [[Rcpp::export]]
arma::sp_mat mat_add(const arma::sp_mat& a, const arma::sp_mat& b) {
  return (a + b);
}

// [[Rcpp::export]]
arma::sp_mat mat_minus(const arma::sp_mat& a, const arma::sp_mat& b) {
  return (a - b);
}
