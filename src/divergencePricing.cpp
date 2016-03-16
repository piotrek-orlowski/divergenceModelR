#include <RcppArmadillo.h>
//#include "affineCF.h"
#include "affineModelR_RcppExports.h"
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
arma::cube divergenceSwapRateCpp(const arma::vec p, const arma::cube coeffs, const arma::mat stateMat){
  
  arma::cube cfVals(p.n_elem, coeffs.n_rows, stateMat.n_rows, arma::fill::zeros);
  arma::cube cfAndDerivVals(p.n_elem, coeffs.n_rows, stateMat.n_rows * 4, arma::fill::zeros);
  
  // First calculate CF values and less 1
  cfVals = affineModelR::affineCFevalCpp(coeffs.slices(0,stateMat.n_cols),stateMat,false);
  cfVals -= arma::ones(arma::size(cfVals));
  // Rescale for p different from 0 and 1
  for(unsigned int pp = 0; pp < p.n_elem ;pp++){
    double ploc = p(pp);
    if(ploc!=0.0 && ploc!=1.0){
      cfVals.subcube(arma::span(pp,pp), arma::span::all, arma::span::all) *= 1.0/(ploc*(ploc-1.0));
    }
  }
  
  // 0 and 1 require special treatment (derivatives of the CF)
  if(any(p == 0 || p == 1)){
    // locate 0 and 1
    arma::vec pZeros(p.n_elem,arma::fill::zeros);
    arma::vec pOnes(p.n_elem,arma::fill::ones);
    arma::uvec whichZero = find(p == pZeros);
    arma::uvec whichOne = find(p == pOnes);
    // calculate derivatives
    cfAndDerivVals = affineModelR::affineCFderivsEvalCpp(coeffs,stateMat);
    // assign in zeros
    if(whichZero.n_elem > 0){
      arma::cube tempCube = -cfAndDerivVals.subcube(arma::span(whichZero(0),whichZero(0)), arma::span::all, arma::span(stateMat.n_rows,2*stateMat.n_rows-1));
      cfVals.subcube(arma::span(whichZero(0),whichZero(0)), arma::span::all, arma::span::all) = tempCube;  
    }
    // assign in ones
    if(whichOne.n_elem > 0){
      arma::cube tempCube = cfAndDerivVals.subcube(arma::span(whichOne(0),whichOne(0)), arma::span::all, arma::span(stateMat.n_rows,2*stateMat.n_rows-1));
      cfVals.subcube(arma::span(whichOne(0),whichOne(0)), arma::span::all, arma::span::all) = tempCube;
    }
  }
  return cfVals;
}

//' @export
// [[Rcpp::export]]
arma::cube skewnessSwapRateCpp(const arma::vec p, const arma::cube coeffs, const arma::mat stateMat){
  
  arma::cube cfVals(p.n_elem, coeffs.n_rows, stateMat.n_rows, arma::fill::zeros);
  arma::cube cfAndDerivVals(p.n_elem, coeffs.n_rows, stateMat.n_rows * 4, arma::fill::zeros);
  arma::cube swapRates(p.n_elem, coeffs.n_rows, stateMat.n_rows, arma::fill::zeros);
  
  // First calculate CF values and less 1
  cfAndDerivVals = affineModelR::affineCFderivsEvalCpp(coeffs,stateMat);
  cfVals = cfAndDerivVals(arma::span::all, arma::span::all, arma::span(0,stateMat.n_rows-1));
  cfVals -= arma::ones(arma::size(cfVals));
  
  // collect cf derivs in more convenient containers
  arma::cube cfFirstDeriv = cfAndDerivVals(arma::span::all, arma::span::all, arma::span(stateMat.n_rows, 2*stateMat.n_rows-1));
  
  // Rescale for p different from 0 and 1
  for(unsigned int pp = 0; pp < p.n_elem ;pp++){
    double ploc = p(pp);
    if(ploc!=0.0 && ploc!=1.0){
      cfVals.subcube(arma::span(pp,pp), arma::span::all, arma::span::all) *= -(2*ploc-1)/(pow(ploc,2.0)*pow(ploc-1.0,2.0)); // -(2*p0-1)/(p0^2*(p0-1)^2)
      cfFirstDeriv.subcube(arma::span(pp,pp), arma::span::all, arma::span::all) *= 1.0/(ploc*(ploc-1.0));
    }
  }
  
  // Calculate swap rates for regular cases
  swapRates = cfFirstDeriv + cfVals;
  
  // 0 and 1 require special treatment (derivatives of the CF)
  if(any(p == 0 || p == 1)){
    // locate 0 and 1
    arma::vec pZeros(p.n_elem,arma::fill::zeros);
    arma::vec pOnes(p.n_elem,arma::fill::ones);
    arma::uvec whichZero = find(p == pZeros);
    arma::uvec whichOne = find(p == pOnes);
    // assign in zeros
    if(whichZero.n_elem > 0){
      arma::cube tempFirstDerivCube = -cfAndDerivVals.subcube(arma::span(whichZero(0),whichZero(0)), arma::span::all, arma::span(stateMat.n_rows,2*stateMat.n_rows-1));
      arma::cube tempSecondDerivCube = -0.5*cfAndDerivVals.subcube(arma::span(whichZero(0),whichZero(0)), arma::span::all, arma::span(2*stateMat.n_rows,3*stateMat.n_rows-1));
      swapRates.subcube(arma::span(whichZero(0),whichZero(0)), arma::span::all, arma::span::all) = tempFirstDerivCube + tempSecondDerivCube;  
    }
    // assign in ones
    if(whichOne.n_elem > 0){
      arma::cube tempFirstDerivCube = -cfAndDerivVals.subcube(arma::span(whichOne(0),whichOne(0)), arma::span::all, arma::span(stateMat.n_rows,2*stateMat.n_rows-1));
      arma::cube tempSecondDerivCube = 0.5*cfAndDerivVals.subcube(arma::span(whichOne(0),whichOne(0)), arma::span::all, arma::span(2*stateMat.n_rows,3*stateMat.n_rows-1));
      swapRates.subcube(arma::span(whichOne(0),whichOne(0)), arma::span::all, arma::span::all) = tempFirstDerivCube + tempSecondDerivCube;  
    }
  }
  return swapRates;
}

//' @export
// [[Rcpp::export]]
arma::cube quarticitySwapRateCpp(const arma::vec p, const arma::cube coeffs, const arma::mat stateMat){
  
  arma::cube cfVals(p.n_elem, coeffs.n_rows, stateMat.n_rows, arma::fill::zeros);
  arma::cube cfAndDerivVals(p.n_elem, coeffs.n_rows, stateMat.n_rows * 4, arma::fill::zeros);
  arma::cube swapRates(p.n_elem, coeffs.n_rows, stateMat.n_rows, arma::fill::zeros);
  
  // First calculate CF values and less 1
  cfAndDerivVals = affineModelR::affineCFderivsEvalCpp(coeffs,stateMat);
  cfVals = cfAndDerivVals(arma::span::all, arma::span::all, arma::span(0,stateMat.n_rows-1));
  cfVals -= arma::ones(arma::size(cfVals));
  
  // collect cf derivs in more convenient containers
  arma::cube cfFirstDeriv = cfAndDerivVals(arma::span::all, arma::span::all, arma::span(stateMat.n_rows, 2*stateMat.n_rows-1));
  arma::cube cfSecondDeriv = cfAndDerivVals(arma::span::all, arma::span::all, arma::span(2*stateMat.n_rows, 3*stateMat.n_rows-1));
  
  // Rescale for p different from 0 and 1
  for(unsigned int pp = 0; pp < p.n_elem; pp++){
    double ploc = p(pp);
    if((ploc!=0.0) && (ploc!=1.0)){
      cfVals.subcube(arma::span(pp,pp), arma::span::all, arma::span::all) *= (2.0-6.0*ploc+6.0*pow(ploc,2.0))/(pow(ploc,3.0)*pow(ploc-1.0,3.0)); // (2-6*p0+6*p0^2)/(p0^3*(p0-1)^3)
      cfFirstDeriv.subcube(arma::span(pp,pp), arma::span::all, arma::span::all) *= (-2.0*ploc + 6.0*pow(ploc,2.0) - 4.0*pow(ploc,3.0))/(pow(ploc,3.0)*pow(ploc-1.0,3.0)); //(-2*p0+6*p0^2-4*p0^3)/(p0^3*(p0-1)^3)
      cfSecondDeriv.subcube(arma::span(pp,pp), arma::span::all, arma::span::all) *= (pow(ploc,2.0)-2.0*pow(ploc,3.0)+pow(ploc,4.0))/(pow(ploc,3.0)*pow(ploc-1.0,3.0)); //(p0^2-2*p0^3+p0^4)/(p0^3*(p0-1)^3)
    }
  }
  
  // Calculate swap rates for regular cases
  swapRates = cfSecondDeriv + cfFirstDeriv + cfVals;
  
  arma::vec pZeros(p.n_elem,arma::fill::zeros);
  arma::vec pOnes(p.n_elem,arma::fill::ones);
  
  // 0 and 1 require special treatment (derivatives of the CF)
  if(any((p == 0) || (p == 1))){
    // locate 0 and 1
    arma::uvec whichZero = find(p == pZeros);
    arma::uvec whichOne = find(p == pOnes);
    // assign in zeros
    if(whichZero.n_elem > 0){
      arma::cube tempFirstDerivCube = -2.0*cfAndDerivVals.subcube(arma::span(whichZero(0),whichZero(0)), arma::span::all, arma::span(stateMat.n_rows,2*stateMat.n_rows-1));
      arma::cube tempSecondDerivCube = -cfAndDerivVals.subcube(arma::span(whichZero(0),whichZero(0)), arma::span::all, arma::span(2*stateMat.n_rows,3*stateMat.n_rows-1));
      arma::cube tempThirdDerivCube = -1.0/3.0*cfAndDerivVals.subcube(arma::span(whichZero(0),whichZero(0)), arma::span::all, arma::span(3*stateMat.n_rows,4*stateMat.n_rows-1));
      swapRates.subcube(arma::span(whichZero(0),whichZero(0)), arma::span::all, arma::span::all) = tempFirstDerivCube + tempSecondDerivCube + tempThirdDerivCube;  
    }
    // assign in ones
    if(whichOne.n_elem > 0){
      arma::cube tempFirstDerivCube = 2.0*cfAndDerivVals.subcube(arma::span(whichOne(0),whichOne(0)), arma::span::all, arma::span(stateMat.n_rows,2*stateMat.n_rows-1));
      arma::cube tempSecondDerivCube = -cfAndDerivVals.subcube(arma::span(whichOne(0),whichOne(0)), arma::span::all, arma::span(2*stateMat.n_rows,3*stateMat.n_rows-1));
      arma::cube tempThirdDerivCube = 1.0/3.0*cfAndDerivVals.subcube(arma::span(whichOne(0),whichOne(0)), arma::span::all, arma::span(3*stateMat.n_rows,4*stateMat.n_rows-1));
      swapRates.subcube(arma::span(whichOne(0),whichOne(0)), arma::span::all, arma::span::all) = tempFirstDerivCube + tempSecondDerivCube + tempThirdDerivCube;  
    }
  }
  return swapRates;
}