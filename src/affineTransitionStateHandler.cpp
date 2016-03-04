#include <RcppArmadillo.h>
#include "meanVec.h"
#include "covMat.h"
#include "linCombMean.h"

using namespace std;

//' @export
// [[Rcpp::export]]
Rcpp::List affineTransitionStateHandler(arma::mat stateMat, Rcpp::List modelParameters){
  // We propagate the past and current state. As such, stateMat.n_rows = 2*N.factors
  int Nf = stateMat.n_rows/2;
  
  Rcpp::List res;
  
  Rcpp::List meanList = modelParameters["mean.vec"];
  Rcpp::List covList = modelParameters["cov.array"];
  
  arma::mat meanMat(stateMat);
  
  for(int kcol = 0; kcol < meanMat.n_cols; kcol++){
    arma::mat locState = stateMat.col(kcol);
    meanMat.submat(0,kcol,Nf-1,kcol) = meanVecFun(meanList, locState.rows(0,Nf-1));
    meanMat.submat(Nf,kcol,2*Nf-1,kcol) = stateMat.submat(0,0,Nf-1,0); // for reasonable increment handling
    // meanMat.col(kcol) = meanVecFun(meanList, stateMat.col(kcol));
  }
  
  arma::vec covDim(2);
  covDim(0) = Nf;
  covDim(1) = Nf;
  
  arma::mat covMat = covMatFun(covList, covDim, stateMat.submat(0,0,Nf-1,0));
  covMat = covMat - meanMat.submat(0,0,Nf-1,0) * meanMat.submat(0,0,Nf-1,0).t();
  
  arma::mat covMatAll(2*Nf,2*Nf,arma::fill::zeros);
  covMatAll.submat(0,0,Nf-1,Nf-1) = covMat;
  
  for(int kcol = 0; kcol < meanMat.n_cols; kcol++){
    meanMat.submat(0,kcol,Nf-1,kcol) += stateMat.submat(0,kcol,Nf-1,kcol);
    // meanMat.col(kcol) += stateMat.col(kcol);
  }
  
  res = Rcpp::List::create(Rcpp::Named("stateVec") = meanMat, Rcpp::Named("procNoiseMat") = covMatAll);
  
  return res;
}