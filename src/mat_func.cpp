//' @name mad_add
//' @title matrix addition
//' @param a a matrix
//' @param b a matrix
//' @return a matrix

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS // try removing it

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Eigen/Eigen>
using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat mat_add(const arma::sp_mat& a, const arma::sp_mat& b) {
  return (a + b);
}



//' @name mad_minus
//' @title matrix subtraction
//' @param a a matrix
//' @param b a matrix
//' @return a matrix

// [[Rcpp::export]]
arma::sp_mat mat_minus(const arma::sp_mat& a, const arma::sp_mat& b) {
  return (a - b);
}
