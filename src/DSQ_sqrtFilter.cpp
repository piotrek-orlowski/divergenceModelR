#include <RcppArmadillo.h>
#include "ukfClass.h"
#include "handlers.h"
#include "stateAdmissibilityControl.h"
using namespace std;
using namespace Rcpp;

//' @name filterCaller
//' @title C++ filter call
//' @description This function wraps the ukfClass construction and exposes the call to R. Should not be used directly by the user, rather its exported R interface should be given as an argument to \code{\link{modelLikelihood}}
//' @details \code{DSQ_sqrtFilter} works on divergence, skewness and quarticity
//' @export
//[[Rcpp::export]]
List DSQ_sqrtFilter(const arma::mat dataMat, const arma::vec initState, const arma::mat initProcCov, const List modelParams){
  
  // Set state handlers  
  stateHandler transitionPtr = &affineTransitionStateHandler;
  stateHandler observationPtr = &affineObservationStateHandler;
  stateControl controlPtr = &affineAdditiveStateController;
  
  // Instantiate filter instance
  ukfClass filterInstance(dataMat, initState, initProcCov, transitionPtr, observationPtr, controlPtr, modelParams);
  
  // Run filter
  filterInstance.filterSqrtAdditiveNoise();
  
  // Set return variable
  List res = List::create(Named("estimState") = filterInstance.getStateMat(), Named("stateCovCube") = filterInstance.getCovCube(), Named("logL") = filterInstance.getLogL(), Named("predMat") = filterInstance.getPredMat());
  
  return res;
}

//' @describeIn filterCaller
//' @details \code{D_sqrtFilter} works on divergence
//' @export
//[[Rcpp::export]]
List D_sqrtFilter(const arma::mat dataMat, const arma::vec initState, const arma::mat initProcCov, const List modelParams){
  
  // Set state handlers  
  stateHandler transitionPtr = &affineTransitionStateHandler;
  stateHandler observationPtr = &affineObservationStateHandler_D;
  stateControl controlPtr = &affineAdditiveStateController;
  
  // Instantiate filter instance
  ukfClass filterInstance(dataMat, initState, initProcCov, transitionPtr, observationPtr, controlPtr, modelParams);
  
  // Run filter
  filterInstance.filterSqrtAdditiveNoise();
  
  // Set return variable
  List res = List::create(Named("estimState") = filterInstance.getStateMat(), Named("stateCovCube") = filterInstance.getCovCube(), Named("logL") = filterInstance.getLogL());
  
  return res;
}

//' @describeIn filterCaller
//' @details \code{DS_sqrtFilter} works on divergence
//' @export
//[[Rcpp::export]]
List DS_sqrtFilter(const arma::mat dataMat, const arma::vec initState, const arma::mat initProcCov, const List modelParams){
  
  // Set state handlers  
  stateHandler transitionPtr = &affineTransitionStateHandler;
  stateHandler observationPtr = &affineObservationStateHandler_DS;
  stateControl controlPtr = &affineAdditiveStateController;
  
  // Instantiate filter instance
  ukfClass filterInstance(dataMat, initState, initProcCov, transitionPtr, observationPtr, controlPtr, modelParams);
  
  // Run filter
  filterInstance.filterSqrtAdditiveNoise();
  
  // Set return variable
  List res = List::create(Named("estimState") = filterInstance.getStateMat(), Named("stateCovCube") = filterInstance.getCovCube(), Named("logL") = filterInstance.getLogL());
  
  return res;
}