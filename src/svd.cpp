#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

/**
 * SVD
 *
 * @param X Matrix
 * @return eigenvectors and eigenvalues
 */
// [[Rcpp::export(.quick_svd)]]
List quick_svd(arma::mat X) {
  int p = std::min(X.n_cols, X.n_rows);
  List ret ;
  arma::colvec index (p, arma::fill::ones);

  arma::mat U ;
  arma::vec s ;
  arma::mat V ;
  svd_econ(U, s, V, X, "left") ;
  s = pow(s, 2);

for (int i = s.size() - 1; i>=0; --i) {
   if (s(i) > 1e-15) {
   break;
 } else {
   index(i) = 0;
}
}

  arma::uvec index2 = find(index);
  arma::vec s_sub = s.rows(index2);
  arma::mat U_sub = U.cols(index2);
  ret["evalues"] = s_sub;
  ret["evectors"] = U_sub ;
  return(ret) ;
}
