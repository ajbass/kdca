#ifndef CROSSPROD_H
#define CROSSPROD_H
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Rcpp;

arma::mat pairwise_prod(arma::mat x);
Eigen::MatrixXd fastmult(Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);

#endif
