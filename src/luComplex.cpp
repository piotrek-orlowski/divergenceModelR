#include <R_ext/Rdynload.h>
#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

extern "C" SEXP luComplex(SEXP xs) {
  try{
    List matList(xs);
    NumericMatrix rMat = as<NumericMatrix>(matList["rMat"]);
    NumericMatrix iMat = as<NumericMatrix>(matList["iMat"]);
    int n = rMat.nrow();
    int p = rMat.ncol();
    
    mat rM = mat(rMat.begin(),n,p);
    mat iM = mat(iMat.begin(),n,p);
    
    cx_mat M = cx_mat(rM,iM);
    
    cx_mat L,U;
    
    lu(L,U,M);
    
    rM = real(U);
    iM = imag(U);
    
    return List::create(Named("rU")= rM, Named("iU") = iM);
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}