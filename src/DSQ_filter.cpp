#include <RcppArmadillo.h>
#include "ukfClass.h"
#include "handlers.h"
using namespace std;
using namespace Rcpp;

//' @describeIn filterCaller
//' @title C++ filter call
//' @description This function wraps the ukfClass construction and exposes the call to R. Should not be used directly by the user, rather its exported R interface should be given as an argument to \code{\link{modelLikelihood}}
//[[Rcpp::export]]
List DSQ_filter(const arma::mat dataMat, const arma::vec initState, const arma::mat initProcCov, const List modelParams){

// Set state handlers  
stateHandler transitionPtr = &affineTransitionStateHandler;
stateHandler observationPtr = &affineObservationStateHandler;

// Instantiate filter instance
ukfClass filterInstance(dataMat, initState, initProcCov, transitionPtr, observationPtr, modelParams);

// Run filter
filterInstance.filterAdditiveNoise();

// Set return variable
List res = List::create(Named("estimState") = filterInstance.getStateMat(), Named("stateCovCube") = filterInstance.getCovCube(), Named("logL") = filterInstance.getLogL());

return res;
}