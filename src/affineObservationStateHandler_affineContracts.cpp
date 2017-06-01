#include <RcppArmadillo.h>
#include "meanVec.h"
#include "covMat.h"
#include "linCombMean.h"
#include "divergencePricing.h"
#include "divObsNoiseMat.h"

using namespace std;

//' @export
// [[Rcpp::export]]
Rcpp::List affineObservationStateHandler_affineContracts(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount){

  // In state mat, additionally, to state values, we carry the past filtered state
  // because it is necessary to calculate effects of stock return observation
  int Nf = stateMat.n_rows/2;
  
  // Evaluate divergence prices
  arma::cube cfCoeffs = Rcpp::as<arma::cube>(modelParameters["cfCoeffs"]);
  arma::vec pVecUnique = Rcpp::as<arma::vec>(modelParameters["pVecUnique"]);
  arma::vec tVecUnique = Rcpp::as<arma::vec>(modelParameters["tVecUnique"]);
  arma::vec pVec = Rcpp::as<arma::vec>(modelParameters["pVec"]);
  arma::vec tVec = Rcpp::as<arma::vec>(modelParameters["tVec"]);
  
  arma::cube divPrices = divergenceSwapRateCpp(pVecUnique, cfCoeffs, stateMat.rows(0,Nf-1).t());
  
  arma::cube skewPrices = skewnessSwapRateCpp(pVecUnique, cfCoeffs, stateMat.rows(0,Nf-1).t());
  
  // create return matrix
  arma::mat yhat(pVec.n_elem,stateMat.n_cols);
  
  // Cycle over state variable columns
  for(unsigned int kcol=0; kcol < stateMat.n_cols; kcol++){
    // Cycle through the representation of the contract specification matrix and construct final prices from priced single instruments
    for(unsigned int prow=0; prow < pVec.n_elem; prow++){
      // find p and t values in the unique vectors to know what to pick
      double loc_p = pVec(prow);
      double loc_t = tVec(prow);
      
      arma::vec loc_p_vec(pVecUnique.n_elem);
      loc_p_vec.fill(loc_p);
      arma::vec loc_t_vec(tVecUnique.n_elem);
      loc_t_vec.fill(loc_t);
      arma::uvec which_p = find(loc_p == loc_p_vec);
      arma::uvec which_t = find(loc_t == loc_t_vec);
      
      // pick and calculate necessary affine swap rate
      double loc_div_price = divPrices(which_p(0),which_t(0),kcol);
      double loc_skew_price = skewPrices(which_p(0),which_t(0),kcol);
      if(loc_p == 0.0){
        yhat(prow, kcol) = loc_div_price;
      } else if(loc_p == 1.0){
        // not implemented yet
      } else {
        yhat(prow, kcol) = (loc_skew_price - (1.0 - 2.0*loc_p)/(loc_p*(loc_p-1.0))*loc_div_price) / (loc_div_price + 1.0);
      }
    }
  }
  
  // Handle the noise covariance matrix. The stock ``observation'' noise is
  // uncorrelated with portfolio observation noise.
  arma::mat obsNoiseMat(pVec.n_elem,pVec.n_elem,arma::fill::zeros);
  
  // 
  // // extract observation noise params
  // arma::vec bVec = Rcpp::as<arma::vec>(modelParameters["bVec"]);
  // arma::vec cVec = Rcpp::as<arma::vec>(modelParameters["cVec"]);
  // // double spotVol = arma::accu(stateMat.col(0));
  // double spotVol = arma::accu(stateMat.submat(0,0,Nf-1,0));
  // obsNoiseMat(arma::span(1,obsNoiseMat.n_rows-1),arma::span(1,obsNoiseMat.n_cols-1)) = divModelObsNoiseMat(cVec,bVec,spotVol,tVecUnique,U);
  SEXP divBigMat_SEXP(modelParameters["divNoiseCube"]);
  arma::cube divBigMat = as<arma::cube>(divBigMat_SEXP);
  
  // var-cov matrix scaling
  arma::vec errSdParVec = as<arma::vec>(modelParameters["errSdParVec"]);
  
  obsNoiseMat(arma::span(0,obsNoiseMat.n_rows-1),arma::span(0,obsNoiseMat.n_cols-1)) = errSdParVec(0) * divBigMat.slice(iterCount);
  // Initialize return list
  Rcpp::List res = Rcpp::List::create(Rcpp::Named("yhat") = yhat, Rcpp::Named("obsNoiseMat") = obsNoiseMat);
  
  return res;
}