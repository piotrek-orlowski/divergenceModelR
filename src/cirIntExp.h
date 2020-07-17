#ifndef _affineCpp_CIRINTEXP
#define _affineCpp_CIRINTEXP

#include "RcppArmadillo.h"

using namespace std;

RcppExport SEXP cirIntExp(SEXP volParams, SEXP t0, SEXP t1, SEXP v0);

#endif