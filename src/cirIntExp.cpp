// quick implementation of the integrated expectation of CIR process (without jumps)
#include "cirIntExp.h"
#include "Rcpp.h"

using namespace arma;
using namespace std;
using namespace Rcpp;

RcppExport SEXP cirIntExp(SEXP volParams, SEXP t0, SEXP t1, SEXP v0)
{
  try{
    NumericVector volP(volParams);
    double kpp = volP[0];
    double eta = volP[1];
    double t0c = as<double>(t0);
    double t1c = as<double>(t1);
    double v0c = as<double>(v0);
    // double lmb = volParams(2);
    
    //  	double res = exp( -t1c * kpp ) * ( t1c * eta * kpp * exp( t1c * kpp) - v0c + eta);
    //  	res -= exp( - t0c * kpp ) * ( t0c * eta * kpp *exp( t0c * kpp ) - v0c + eta);
    //  	res /= kpp;
    
    double res = eta*(t1c-t0c)+1.0/kpp*exp(-kpp*t1c)*(eta - v0c)+1.0/kpp*exp(-kpp*t0c)*(-eta + v0c);
    
    return wrap(res);
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return -999;
}