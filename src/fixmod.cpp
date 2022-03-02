
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

Rcpp::List fixmod(Rcpp::List mod, arma::rowvec coeffs, Rcpp::List X){
  
  arma::sp_mat Xmat = Rcpp::as<arma::sp_mat>(X["design"]);
  
  // original stuff
  arma::rowvec y = Rcpp::as<arma::rowvec>(mod["y"]);
  arma::rowvec id = Rcpp::as<arma::rowvec>(mod["id"]);
  int maxid = Rcpp::as<int>(mod["max.id"]);
  
  // recalculate data for output
  arma::mat ilinpred = Xmat * coeffs.t();
  arma::mat ifitted = arma::exp(ilinpred) / (arma::exp(ilinpred) + 1);
  arma::rowvec linpred = arma::vectorise(ilinpred,1);
  arma::rowvec fitted = arma::vectorise(ifitted,1);
  
  // output
  return Rcpp::List::create(Rcpp::Named("y")=y,
                            Rcpp::Named("fitted.values")=fitted,
                            Rcpp::Named("linear.predictors")=linpred,
                            Rcpp::Named("id")=id,
                            Rcpp::Named("max.id")=maxid,
                            Rcpp::Named("coefficients")=coeffs);
  
}
