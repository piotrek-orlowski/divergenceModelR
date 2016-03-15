#ifndef STATE_ADMISSIBILITY_CONTROL_H
#define STATE_ADMISSIBILITY_CONTROL_H
#include <RcppArmadillo.h>

using namespace std;

arma::mat affineAdditiveStateController(arma::mat& stateMat);

#endif