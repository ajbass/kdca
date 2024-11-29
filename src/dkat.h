#ifndef DKAT_CPP_H
#define DKAT_CPP_H
#include <RcppArmadillo.h>
using namespace Rcpp;

double dkat_statistic(arma::mat Ux, arma::mat dx,arma::mat Uy,arma::mat dy);
double dkat(arma::mat Ux, arma::mat dx,arma::mat Uy,arma::mat dy);

#endif
