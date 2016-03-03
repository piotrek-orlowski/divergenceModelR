#ifndef HANDLERS_H
#define HANDLERS_H
#include <RcppArmadillo.h>

Rcpp::List affineObservationStateHandler(arma::mat stateMat, Rcpp::List modelParameters);

Rcpp::List affineTransitionStateHandler(arma::mat stateMat, Rcpp::List modelParameters);

#endif