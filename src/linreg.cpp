#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(.quick_lm_boot_cpp)]]
List quick_lm_boot_cpp(SEXP Xs, SEXP Ys) {
  Rcpp::NumericMatrix Xr(Xs);
  Rcpp::NumericMatrix Yr(Ys);
  int n = Xr.nrow(), k = Xr.ncol(), p = Yr.ncol();

  arma::mat X(Xr.begin(), n, k, false);
  arma::mat Y(Yr.begin(), n, p, false);

  // fit model y ~ X, extract residuals
  arma::mat H = X * inv(X.t() * X) * X.t();
 // arma::mat coef = arma::solve(X.t() * X, X.t() * Y);
  arma::mat fit =  H * Y;
  arma::mat res  = Y - fit;
  arma::vec hatvals = H.diag();

  List ret;
  ret["fitted"] = fit;
  ret["res"] = res ;
  ret["hatvals"] = hatvals ;
  return(ret);
}
