// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pairwise_prod
arma::mat pairwise_prod(arma::mat x);
RcppExport SEXP _kdca_pairwise_prod(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(pairwise_prod(x));
    return rcpp_result_gen;
END_RCPP
}
// dkat_statistic
long double dkat_statistic(arma::mat Ux, arma::mat dx, arma::mat Uy, arma::mat dy);
RcppExport SEXP _kdca_dkat_statistic(SEXP UxSEXP, SEXP dxSEXP, SEXP UySEXP, SEXP dySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Ux(UxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Uy(UySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dy(dySEXP);
    rcpp_result_gen = Rcpp::wrap(dkat_statistic(Ux, dx, Uy, dy));
    return rcpp_result_gen;
END_RCPP
}
// dkat
long double dkat(arma::mat Ux, arma::mat dx, arma::mat Uy, arma::mat dy);
RcppExport SEXP _kdca_dkat(SEXP UxSEXP, SEXP dxSEXP, SEXP UySEXP, SEXP dySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Ux(UxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Uy(UySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dy(dySEXP);
    rcpp_result_gen = Rcpp::wrap(dkat(Ux, dx, Uy, dy));
    return rcpp_result_gen;
END_RCPP
}
// quick_lm_boot_cpp
List quick_lm_boot_cpp(SEXP Xs, SEXP Ys);
RcppExport SEXP _kdca_quick_lm_boot_cpp(SEXP XsSEXP, SEXP YsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Ys(YsSEXP);
    rcpp_result_gen = Rcpp::wrap(quick_lm_boot_cpp(Xs, Ys));
    return rcpp_result_gen;
END_RCPP
}
// quick_svd
List quick_svd(arma::mat X);
RcppExport SEXP _kdca_quick_svd(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(quick_svd(X));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_kdca_pairwise_prod", (DL_FUNC) &_kdca_pairwise_prod, 1},
    {"_kdca_dkat_statistic", (DL_FUNC) &_kdca_dkat_statistic, 4},
    {"_kdca_dkat", (DL_FUNC) &_kdca_dkat, 4},
    {"_kdca_quick_lm_boot_cpp", (DL_FUNC) &_kdca_quick_lm_boot_cpp, 2},
    {"_kdca_quick_svd", (DL_FUNC) &_kdca_quick_svd, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_kdca(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}