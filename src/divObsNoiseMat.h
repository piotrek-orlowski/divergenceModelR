#ifndef DIVOBSNOISEMAT_H
#define DIVOBSNOISEMAT_H
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;

arma::mat divModelObsNoiseMat(const arma::vec& corrs, arma::vec& bpars, double spotVar, const arma::vec& matVec, const int U);
#endif