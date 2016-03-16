#ifndef COVMAT_H
#define COVMAT_H

#include "RcppArmadillo.h"

arma::mat covMatFun(Rcpp::List covListS, const arma::uvec covDim, const arma::vec currVol);

#endif