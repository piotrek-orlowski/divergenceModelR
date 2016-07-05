#include <RcppArmadillo.h>
#include "ukfClass.h"
#include "handlers.h"
#include "stateAdmissibilityControl.h"
using namespace std;
using namespace Rcpp;

//' @describeIn filterCaller
//' @details \code{portfolio_sqrtFilter} works on whatever portfolios you supply
//' @export
//[[Rcpp::export]]
List portfolio_sqrtFilter(const arma::mat dataMat, const arma::vec initState, const arma::mat initProcCov, const List modelParams){

// Set state handlers  
stateHandler transitionPtr = &affineTransitionStateHandler;
stateHandler observationPtr = &affineObservationStateHandler_optionPortfolios;
stateControl controlPtr = &affineAdditiveStateController;

// Instantiate filter instance
ukfClass filterInstance(dataMat, initState, initProcCov, transitionPtr, observationPtr, controlPtr, modelParams);

// Run filter
filterInstance.filterSqrtAdditiveNoise();

// Set return variable
List res = List::create(Named("estimState") = filterInstance.getStateMat(), Named("stateCovCube") = filterInstance.getCovCube(), Named("logL") = filterInstance.getLogL(), Named("predMat") = filterInstance.getPredMat());

return res;
}