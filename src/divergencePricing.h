#ifndef DIVPRICING_H
#define DIVPRICING_H
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::cube divergenceSwapRateCpp(const arma::vec& p, const arma::cube& coeffs, const arma::mat& stateMat);

arma::cube skewnessSwapRateCpp(const arma::vec& p, const arma::cube& coeffs, const arma::mat& stateMat);

arma::cube quarticitySwapRateCpp(const arma::vec& p, const arma::cube& coeffs, const arma::mat& stateMat);
#endif