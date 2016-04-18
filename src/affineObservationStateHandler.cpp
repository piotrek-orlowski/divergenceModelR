#include <RcppArmadillo.h>
#include "meanVec.h"
#include "covMat.h"
#include "linCombMean.h"
#include "divergencePricing.h"
#include "divObsNoiseMat.h"

using namespace std;

//' @export
// [[Rcpp::export]]
Rcpp::List affineObservationStateHandler(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount){

// In state mat, additionally, to state values, we carry the past filtered state
// because it is necessary to calculate effects of stock return observation
int Nf = stateMat.n_rows/2;

// Grab stock parameters  
Rcpp::List stockParams = modelParameters["stockParams"];

// First handle the stock. The observation equation uses the beta of the stock 
// on volatility factors. The noise is the residual variance. First compute 
// expected stock returns at all states. Then calculate the betas at the
// original state (stateMat.col(0).)

Rcpp::List stockParamsMeanVec = stockParams["mean.vec"];
Rcpp::List stockMeanIndividual = Rcpp::List::create(stockParamsMeanVec[0]);

// Initialize a vector for stock means, length = number of states in evaluation
arma::vec stockMeans(stateMat.n_cols);
stockMeans.fill(0.0);

// Calculate expected stock returns at previous state
for(unsigned int kcol = 0; kcol < stateMat.n_cols; kcol++){
  stockMeans(kcol) = arma::as_scalar(meanVecFun(stockMeanIndividual , stateMat.submat(Nf,kcol,2*Nf-1,kcol)));
  // stockMeans(kcol) = arma::as_scalar(meanVecFun(stockMeanIndividual , stateMat.col(kcol)));
}

// Start working on the covariance matrix between stock return and states.
arma::vec stockAndVolMeans = meanVecFun(stockParamsMeanVec, stateMat.submat(Nf,0,2*Nf-1,0));
// arma::vec stockAndVolMeans = meanVecFun(stockParamsMeanVec, stateMat.col(0));
arma::uvec covDim(2);
covDim.fill(stockParamsMeanVec.length());

// Covariance matrix and beta of stock on vol factors
arma::mat stockAndVolCov = covMatFun(stockParams["cov.array"], covDim, stateMat.submat(Nf,0,2*Nf-1,0)) - stockAndVolMeans * stockAndVolMeans.t();
arma::mat stockAndVolBeta = stockAndVolCov.submat(0,1,0,stockAndVolCov.n_cols-1) * arma::inv(stockAndVolCov.submat(1,1,stockAndVolCov.n_rows-1,stockAndVolCov.n_cols-1));

// This is the residual volatility of the stock that is not coming from vol
// co-movement. This will go into the observation noise matrix position (1,1)
arma::mat stockNoise = stockAndVolCov(0,0) - stockAndVolBeta * stockAndVolCov.submat(1,1,stockAndVolCov.n_rows-1,stockAndVolCov.n_cols-1) * stockAndVolBeta.t();

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

arma::mat yhat(1+U*T*3,stateMat.n_cols);

arma::mat tempPrices(U*T,1,arma::fill::zeros);
arma::mat tempDivergencePrices(U*T,1,arma::fill::zeros);
for(unsigned int kcol=0; kcol < stateMat.n_cols; kcol++){
  // Write mean returns
  yhat(0,kcol) = stockMeans(kcol);
  yhat(0,kcol) += arma::as_scalar(stockAndVolBeta * (stateMat.submat(0,kcol,Nf-1,kcol) - stateMat.submat(Nf,kcol,2*Nf-1,kcol)));
  // Write divergence prices
  tempPrices.reshape(U,T);
  tempPrices = divPrices.slice(kcol);
  tempPrices.reshape(U*T,1);
  // save divergence prices for standardising skewness and kurtosis
  tempDivergencePrices = tempPrices;
  // annualise divergence prices
  yhat(arma::span(1,U*T),kcol) = tempPrices % arma::pow(tVecPricing,-1.0);

    // Write skewness prices
  tempPrices.reshape(U,T);
  tempPrices = skewPrices.slice(kcol);
  tempPrices.reshape(U*T,1);
  // standardise skewness prices
  tempPrices /= arma::pow(tempDivergencePrices,1.5);

  yhat(arma::span(U*T+1,2*U*T),kcol) = tempPrices;
  // Write quarticity prices
  tempPrices.reshape(U,T);
  tempPrices = quartPrices.slice(kcol);
  tempPrices.reshape(U*T,1);
  // standardise kurtosis prices
  tempPrices /= arma::pow(tempDivergencePrices,2.0);

  yhat(arma::span(2*U*T+1,3*U*T),kcol) = tempPrices;
}

// Handle the noise covariance matrix. The stock ``observation'' noise is
// uncorrelated with portfolio observation noise.
arma::mat obsNoiseMat(1+3*U*T,1+3*U*T,arma::fill::zeros);
obsNoiseMat(0,0) = stockNoise(0,0);
// 
// // extract observation noise params
// arma::vec bVec = Rcpp::as<arma::vec>(modelParameters["bVec"]);
// arma::vec cVec = Rcpp::as<arma::vec>(modelParameters["cVec"]);
// // double spotVol = arma::accu(stateMat.col(0));
// double spotVol = arma::accu(stateMat.submat(0,0,Nf-1,0));
// obsNoiseMat(arma::span(1,obsNoiseMat.n_rows-1),arma::span(1,obsNoiseMat.n_cols-1)) = divModelObsNoiseMat(cVec,bVec,spotVol,tVec,U);
SEXP divBigMat_SEXP(modelParameters["divNoiseCube"]);
arma::cube divBigMat = as<arma::cube>(divBigMat_SEXP);
obsNoiseMat(arma::span(1,obsNoiseMat.n_rows-1),arma::span(1,obsNoiseMat.n_cols-1)) = divBigMat.slice(iterCount);

// Initialize return list
Rcpp::List res = Rcpp::List::create(Rcpp::Named("yhat") = yhat, Rcpp::Named("obsNoiseMat") = obsNoiseMat);

return res;
}