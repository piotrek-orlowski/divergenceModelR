#include <RcppArmadillo.h>
#include "meanVec.h"
#include "covMat.h"
#include "linCombMean.h"

using namespace std;

//' @export
// [[Rcpp::export]]
Rcpp::List affineTransitionStateHandler(arma::mat stateMat, Rcpp::List modelParameters){
  
  Rcpp::List res;
  
  Rcpp::List meanList = modelParameters["mean.vec"];
  Rcpp::List covList = modelParameters["cov.array"];
  
  arma::mat meanMat(stateMat);
  
  for(int kcol = 0; kcol < meanMat.n_cols; kcol++){
    meanMat.col(kcol) = meanVecFun(meanList, stateMat.col(kcol));
  }
  
  arma::vec covDim(2);
  covDim(0) = meanList.length();
  covDim(1) = meanList.length();
  
  arma::mat covMat = covMatFun(covList, covDim, stateMat.col(0));
  covMat = covMat - meanMat * meanMat.t();
  
  res = Rcpp::List::create(Rcpp::Named("stateVec") = meanMat, Rcpp::Named("procNoiseMat") = covMat);
  
  return res;
}