#include <RcppArmadillo.h>
#include "ukfClass.h"
#include "handlers.h"
#include "stateAdmissibilityControl.h"
using namespace std;
using namespace Rcpp;

//' @useDynLib divergenceModelR
//' @export
//[[Rcpp::export]]
List testAffineFilter(const arma::mat testDataMat, const arma::vec testInitState, const arma::mat testInitProcCov, const List testModelParams){
  
  stateHandler transitionPtr = &affineTransitionStateHandler;
  stateHandler observationPtr = &affineObservationStateHandler;
  stateControl controlPtr = &affineAdditiveStateController;
  
  ukfClass filterInstance(testDataMat, testInitState, testInitProcCov, transitionPtr, observationPtr, controlPtr, testModelParams);
  
  filterInstance.filterAdditiveNoise();
  
  List res;
  
  res = List::create(Named("estimState") = filterInstance.getStateMat(), Named("stateCovCube") = filterInstance.getCovCube(), Named("logL") = filterInstance.getLogL());
  
  return res;
}

//' @export
//[[Rcpp::export]]
List testSqrtAffineFilter(const arma::mat testDataMat, const arma::vec testInitState, const arma::mat testInitProcCov, const List testModelParams){
  
  stateHandler transitionPtr = &affineTransitionStateHandler;
  stateHandler observationPtr = &affineObservationStateHandler;
  stateControl controlPtr = &affineAdditiveStateController;
  
  ukfClass filterInstance(testDataMat, testInitState, testInitProcCov, transitionPtr, observationPtr, controlPtr, testModelParams);
  
  filterInstance.filterSqrtAdditiveNoise();
  
  List res;
  
  res = List::create(Named("estimState") = filterInstance.getStateMat(), Named("stateCovCube") = filterInstance.getCovCube(), Named("logL") = filterInstance.getLogL());
  
  return res;
}