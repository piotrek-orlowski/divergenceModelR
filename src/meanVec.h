#ifndef MEANVEC_H
#define MEANVEC_H

#include "RcppArmadillo.h"

// meanList is the stucture list that will be used by linCombMean. The first element should correspond to the stock, the rest are the hedged portfolio returns. 
// U is the number of portfolio returns (so meanList is length U+1)
// currVol is the current state of the volatility factor

arma::vec meanVecFun(Rcpp::List meanListS, const arma::vec& currVol);

#endif