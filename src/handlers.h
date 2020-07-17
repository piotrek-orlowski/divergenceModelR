#ifndef HANDLERS_H
#define HANDLERS_H
#include <RcppArmadillo.h>

Rcpp::List affineObservationStateHandler(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int);

Rcpp::List affineObservationStateHandler_D(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int);

Rcpp::List affineObservationStateHandler_DS(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int);

Rcpp::List affineTransitionStateHandler(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int);

Rcpp::List affineObservationStateHandler_optionPortfolios(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);

Rcpp::List affineObservationStateHandler_optionPortfolios_noStock(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);

Rcpp::List affineObservationStateHandler_affineContracts(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);

// raw cumulant filters
Rcpp::List affineObservationStateHandler_cumulant_DS(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);

#endif