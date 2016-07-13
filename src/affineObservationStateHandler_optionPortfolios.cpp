#include <RcppArmadillo.h>
#include "meanVec.h"
#include "covMat.h"
#include "linCombMean.h"
#include "divergencePricing.h"
#include "divObsNoiseMat.h"
#include "affineCF.h"
#include "transformOptionPricer.h"

using namespace std;

//' @export
// [[Rcpp::export]]
Rcpp::List affineObservationStateHandler_optionPortfolios(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount){

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
  }
  
  // Start working on the covariance matrix between stock return and states.
  arma::vec stockAndVolMeans = meanVecFun(stockParamsMeanVec, stateMat.submat(Nf,0,2*Nf-1,0));
  arma::uvec covDim(2);
  covDim.fill(stockParamsMeanVec.length());
  
  // Covariance matrix and beta of stock on vol factors
  arma::mat stockAndVolCov = covMatFun(stockParams["cov.array"], covDim, stateMat.submat(Nf,0,2*Nf-1,0)) - stockAndVolMeans * stockAndVolMeans.t();
  arma::mat stockAndVolBeta = stockAndVolCov.submat(0,1,0,stockAndVolCov.n_cols-1) * arma::inv(stockAndVolCov.submat(1,1,stockAndVolCov.n_rows-1,stockAndVolCov.n_cols-1));
  
  // This is the residual volatility of the stock that is not coming from vol
  // co-movement. This will go into the observation noise matrix position (1,1)
  arma::mat stockNoise = stockAndVolCov(0,0) - stockAndVolBeta * stockAndVolCov.submat(1,1,stockAndVolCov.n_rows-1,stockAndVolCov.n_cols-1) * stockAndVolBeta.t();
  
  // Evaluate option prices, form and store option portfolios
  // Recover list of cubes of coefficients for affine characteristic function, dimension length(gl nodes) x T x (Nf+1)
  Rcpp::List cfCoeffs_list = modelParameters["cfCoeffs"];
  // Recover list of strikeMats, dimension T x K x (number of columns in stateMat, i.e. for how many states you're evaluating)
  Rcpp::List strikeMat_list = modelParameters["strikeMats"];
  // Recover list of mkt parameters, each is matrix of size T x (3 or 4) with maturity in first column (in years), interest rate in second, dividend yield in third
  Rcpp::List mkt_list = modelParameters["mkts"];
  // Recover list of option portfolio weight parameters, each is matrix of size (T x K x W with W is number of portfolios per maturity T)
  Rcpp::List wts_list = modelParameters["wts"];
  
  // Recover Gauss-Laguerre quadrature weights
  arma::vec glWts = modelParameters["quadWeights"];
  arma::vec glNodes = modelParameters["quadNodes"];
  
  // Pick objects relevant for current iteration from list above
  arma::cx_cube cfCoeffs = Rcpp::as<arma::cx_cube>(cfCoeffs_list[iterCount]);
  arma::mat strikeMat_base = Rcpp::as<arma::cube>(strikeMat_list[iterCount]);
  arma::cube strikeMat(strikeMat_base.n_rows, strikeMat_base.n_cols, stateMat.n_cols);
  for(int kk = 0; kk < stateMat.n_cols; kk++){
    strikeMat.slice(kk) = strikeMat_base;
  }
  arma::mat mkt  = Rcpp::as<arma::mat>(mkt_list[iterCount]);
  arma::cube wts = Rcpp::as<arma::cube>(wts_list[iterCount]); // cube dimension T x K x W
  
  // define 1i imaginary
  std::complex<double> i1(0,1);
  
  // evaluate characteristic function at test states (upper rows of stateMat, transposed)
  arma::cx_cube cfVals = affineCFevalCpp(cfCoeffs, stateMat.rows(0,Nf-1).t(), false);
  
  // price out of the money options with the use of Gauss-Laguerre quadrature
  arma::cube otm_options = glPricer_cpp(strikeMat, mkt, glWts, glNodes, cfVals, Nf, 0.0, 0.75); // dimension T x K x S
  
  unsigned int W,T,K;
  T = mkt.n_rows; // number of maturities
  K = strikeMat.n_cols; // number of strikes
  W = wts.n_slices; // number of option portfolios per maturity
  
  // declare parts of return objects
  arma::mat yhat(1L+W*T,stateMat.n_cols,arma::fill::zeros);
  arma::vec tempPrices(W*T);
  
  for(unsigned int kcol=0; kcol < stateMat.n_cols; kcol++){
    // Write mean returns
    yhat(0,kcol) = stockMeans(kcol);
    // Form stock updates in terms of deviations from the central prediction value (!!!)
    yhat(0,kcol) += arma::as_scalar(stockAndVolBeta * (stateMat.submat(0,kcol,Nf-1,kcol) - stateMat.submat(0,0,Nf-1,0)));
    for(unsigned int tt = 0; tt < T; tt++){
      for(unsigned int ww = 0; ww < W; ww++){
        // Note: the results will be ordered first by maturity, then by ordering of portfolios
        unsigned int locIndex = tt * W + ww;
        arma::mat locWts_cube = wts(arma::span(tt,tt),arma::span(0,K-1),arma::span(ww,ww));
        arma::mat locOtm_cube = otm_options(arma::span(tt,tt),arma::span(0,K-1),arma::span(kcol,kcol));
        arma::vec locWts = arma::conv_to<arma::vec>::from(locWts_cube);
        arma::vec locOtm = arma::conv_to<arma::vec>::from(locOtm_cube);
        // calculate price of portfolio as scalar product
        tempPrices(locIndex) = arma::as_scalar(locWts.t() * locOtm);
      }
    }
    // assign calculated portfolio prices to output
    yhat(arma::span(1,T*W),kcol) = tempPrices;
  }
  
  // Handle the noise covariance matrix. 
  // create observation noise variance-covariance
  arma::mat obsNoiseMat(1+W*T,1+W*T,arma::fill::zeros);
  //The stock ``observation'' noise is uncorrelated with portfolio observation noise.
  obsNoiseMat(0,0) = stockNoise(0,0);
  
  // extract observation noise correlation matrix
  SEXP divBigMat_SEXP(modelParameters["divNoiseCube"]);
  arma::cube divBigMat = as<arma::cube>(divBigMat_SEXP);
  
  // var-cov matrix scaling
  arma::vec errSdParVec = as<arma::vec>(modelParameters["errSdParVec"]);
  
  obsNoiseMat(arma::span(1,obsNoiseMat.n_rows-1),arma::span(1,obsNoiseMat.n_cols-1)) = errSdParVec(0) * divBigMat.slice(iterCount);
  
  // Initialize return list
  Rcpp::List res = Rcpp::List::create(Rcpp::Named("yhat") = yhat, Rcpp::Named("obsNoiseMat") = obsNoiseMat);
  
  return res;
}