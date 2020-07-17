#include <RcppArmadillo.h>
#include "meanVec.h"
#include "covMat.h"
#include "linCombMean.h"
#include "divergencePricing.h"
#include "divObsNoiseMat.h"

using namespace std;

//' @export
// [[Rcpp::export]]
Rcpp::List affineObservationStateHandler_cumulant(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount){

// In state mat, additionally, to state values, we carry the past filtered state
// because it is necessary to calculate effects of stock return observation
int Nf = stateMat.n_rows/2;

// Evaluate divergence prices
arma::cube cfCoeffs = Rcpp::as<arma::cube>(modelParameters["cfCoeffs"]);
arma::vec pVec = Rcpp::as<arma::vec>(modelParameters["pVec"]);
arma::vec tVec = Rcpp::as<arma::vec>(modelParameters["tVec"]);

arma::cube divPrices = divergenceSwapRateCpp(pVec, cfCoeffs, stateMat.rows(0,Nf-1).t());

arma::cube skewPrices = skewnessSwapRateCpp(pVec, cfCoeffs, stateMat.rows(0,Nf-1).t());

arma::cube quartPrices = quarticitySwapRateCpp(pVec, cfCoeffs, stateMat.rows(0,Nf-1).t());

// Write into return structures
unsigned int T = divPrices.n_cols;
unsigned int U = divPrices.n_rows;

arma::vec tVecPricing(U*T,1);
for(unsigned int uu = 0; uu < U; uu++){
  for(unsigned int tt = 0; tt < T; tt++){
    tVecPricing(uu + tt*U) = tVec(tt);
  }
}

arma::mat yhat(U*T*3,stateMat.n_cols);

arma::mat tempPrices(U*T,1,arma::fill::zeros);
arma::mat tempDivergencePrices(U*T,1,arma::fill::zeros);
for(unsigned int kcol=0; kcol < stateMat.n_cols; kcol++){
  // Write divergence prices
  tempPrices.reshape(U,T);
  tempPrices = divPrices.slice(kcol);
  tempPrices.reshape(U*T,1);
  // save divergence prices for standardising skewness and kurtosis
  tempDivergencePrices = tempPrices;
  // save divergence prices
  yhat(arma::span(0,U*T-1L),kcol) = tempPrices;

    // Write skewness prices
  tempPrices.reshape(U,T);
  tempPrices = skewPrices.slice(kcol);
  tempPrices.reshape(U*T,1);

  yhat(arma::span(U*T,2*U*T-1),kcol) = tempPrices;
  // Write quarticity prices
  tempPrices.reshape(U,T);
  tempPrices = quartPrices.slice(kcol);
  tempPrices.reshape(U*T,1);

  yhat(arma::span(2*U*T,3*U*T-1),kcol) = tempPrices;
}

// Handle the noise covariance matrix. The stock ``observation'' noise is
// uncorrelated with portfolio observation noise.
arma::mat obsNoiseMat(3*U*T,3*U*T,arma::fill::zeros);

// 
// // extract observation noise params
// arma::vec bVec = Rcpp::as<arma::vec>(modelParameters["bVec"]);
// arma::vec cVec = Rcpp::as<arma::vec>(modelParameters["cVec"]);
// // double spotVol = arma::accu(stateMat.col(0));
// double spotVol = arma::accu(stateMat.submat(0,0,Nf-1,0));
// obsNoiseMat(arma::span(1,obsNoiseMat.n_rows-1),arma::span(1,obsNoiseMat.n_cols-1)) = divModelObsNoiseMat(cVec,bVec,spotVol,tVec,U);
SEXP divBigMat_SEXP(modelParameters["divNoiseCube"]);
arma::cube divBigMat = as<arma::cube>(divBigMat_SEXP);

arma::vec errSdParVec = as<arma::vec>(modelParameters["errSdParVec"]);
double sumStates = arma::accu(stateMat.submat(0,0,Nf-1,0));

obsNoiseMat(arma::span(0,obsNoiseMat.n_rows-1),arma::span(0,obsNoiseMat.n_cols-1)) = divBigMat.slice(iterCount);

// Initialize return list
Rcpp::List res = Rcpp::List::create(Rcpp::Named("yhat") = yhat, Rcpp::Named("obsNoiseMat") = obsNoiseMat);

return res;
}

// [[Rcpp::export]]
Rcpp::List affineObservationStateHandler_cumulant_D(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount){
  
  // In state mat, additionally, to state values, we carry the past filtered state
  // because it is necessary to calculate effects of stock return observation
  int Nf = stateMat.n_rows/2;
  
  
  // Evaluate divergence prices
  arma::cube cfCoeffs = Rcpp::as<arma::cube>(modelParameters["cfCoeffs"]);
  arma::vec pVec = Rcpp::as<arma::vec>(modelParameters["pVec"]);
  arma::vec tVec = Rcpp::as<arma::vec>(modelParameters["tVec"]);
  
  arma::cube divPrices = divergenceSwapRateCpp(pVec, cfCoeffs, stateMat.rows(0,Nf-1).t());
  
  // Write into return structures
  unsigned int T = divPrices.n_cols;
  unsigned int U = divPrices.n_rows;
  
  arma::vec tVecPricing(U*T,1);
  for(unsigned int uu = 0; uu < U; uu++){
    for(unsigned int tt = 0; tt < T; tt++){
      tVecPricing(uu + tt*U) = tVec(tt);
    }
  }
  
  arma::mat yhat(U*T,stateMat.n_cols);
  
  arma::mat tempPrices(U*T,1,arma::fill::zeros);
  
  for(unsigned int kcol=0; kcol < stateMat.n_cols; kcol++){
    // Write divergence prices
    tempPrices.reshape(U,T);
    tempPrices = divPrices.slice(kcol);
    tempPrices.reshape(U*T,1);
    // annualise divergence prices
    yhat(arma::span(0,U*T-1),kcol) = tempPrices % arma::pow(tVecPricing,-1.0);
  }
  
  // Handle the noise covariance matrix. The stock ``observation'' noise is
  // uncorrelated with portfolio observation noise.
  arma::mat obsNoiseMat(U*T,U*T,arma::fill::zeros);
  
  // extract observation noise params
  SEXP divBigMat_SEXP(modelParameters["divNoiseCube"]);
  arma::cube divBigMat = as<arma::cube>(divBigMat_SEXP);
  obsNoiseMat(arma::span(0,obsNoiseMat.n_rows-1),arma::span(0,obsNoiseMat.n_cols-1)) = divBigMat.slice(iterCount);
  
  arma::vec errSdParVec = as<arma::vec>(modelParameters["errSdParVec"]);
  double sumStates = arma::accu(stateMat.submat(0,0,Nf-1,0));
  // errSdParVec *= sumStates;
  
  obsNoiseMat(arma::span(0,obsNoiseMat.n_rows-1),arma::span(0,obsNoiseMat.n_cols-1)) = divBigMat.slice(iterCount);
  
  // Initialize return list
  Rcpp::List res = Rcpp::List::create(Rcpp::Named("yhat") = yhat, Rcpp::Named("obsNoiseMat") = obsNoiseMat);
  
  return res;
}

// [[Rcpp::export]]
Rcpp::List affineObservationStateHandler_cumulant_DS(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount){
  
  // In state mat, additionally, to state values, we carry the past filtered state
  // because it is necessary to calculate effects of stock return observation
  int Nf = stateMat.n_rows/2;
  
  // Evaluate divergence prices
  arma::cube cfCoeffs = Rcpp::as<arma::cube>(modelParameters["cfCoeffs"]);
  arma::vec pVec = Rcpp::as<arma::vec>(modelParameters["pVec"]);
  arma::vec tVec = Rcpp::as<arma::vec>(modelParameters["tVec"]);
  
  arma::cube divPrices = divergenceSwapRateCpp(pVec, cfCoeffs, stateMat.rows(0,Nf-1).t());
  
  arma::cube skewPrices = skewnessSwapRateCpp(pVec, cfCoeffs, stateMat.rows(0,Nf-1).t());
  
  // Write into return structures
  unsigned int T = divPrices.n_cols;
  unsigned int U = divPrices.n_rows;
  
  arma::vec tVecPricing(U*T,1);
  for(unsigned int uu = 0; uu < U; uu++){
    for(unsigned int tt = 0; tt < T; tt++){
      tVecPricing(uu + tt*U) = tVec(tt);
    }
  }
  
  arma::mat yhat(U*T*2,stateMat.n_cols);
  
  arma::mat tempPrices(U*T,1,arma::fill::zeros);
  arma::mat tempDivergencePrices(U*T,1,arma::fill::zeros);
  for(unsigned int kcol=0; kcol < stateMat.n_cols; kcol++){
    
    // Write divergence prices
    tempPrices.reshape(U,T);
    tempPrices = divPrices.slice(kcol);
    tempPrices.reshape(U*T,1);
    yhat(arma::span(0, U*T-1), kcol) = tempPrices;
    
    // Write skewness prices
    tempPrices.reshape(U,T);
    tempPrices = skewPrices.slice(kcol);
    tempPrices.reshape(U*T,1);
    
    yhat(arma::span(U*T,2*U*T-1),kcol) = tempPrices;
  }
  
  // Handle the noise covariance matrix. 
  arma::mat obsNoiseMat(2*U*T,2*U*T,arma::fill::zeros);
  
  // extract observation noise params
  SEXP divBigMat_SEXP(modelParameters["divNoiseCube"]);
  arma::cube divBigMat = as<arma::cube>(divBigMat_SEXP);
  obsNoiseMat(arma::span(0,obsNoiseMat.n_rows-1),arma::span(0,obsNoiseMat.n_cols-1)) = divBigMat.slice(iterCount);
  
  // arma::vec errSdParVec = as<arma::vec>(modelParameters["errSdParVec"]);
  // double sumStates = arma::accu(stateMat.submat(0,0,Nf-1,0));
  // errSdParVec *= sumStates;
  
  // obsNoiseMat(arma::span(0,obsNoiseMat.n_rows-1),arma::span(0,obsNoiseMat.n_cols-1)) = divBigMat.slice(iterCount);
  
  // Initialize return list
  Rcpp::List res = Rcpp::List::create(Rcpp::Named("yhat") = yhat, Rcpp::Named("obsNoiseMat") = obsNoiseMat);
  
  return res;
}