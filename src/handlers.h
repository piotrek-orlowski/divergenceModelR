#ifndef HANDLERS_H
#define HANDLERS_H
#include <RcppArmadillo.h>

Rcpp::List affineObservationStateHandler(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int);

Rcpp::List affineTransitionStateHandler(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int);

#endif