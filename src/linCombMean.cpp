// Calculate the mean given the linComb structure
// Calculate the mean coming from the linComb structure (note that the linComb calculates accurately the conditional log c.f., but we need to Taylor expand it back)
#include "linCombMean.h"
#include "RcppArmadillo.h"

using namespace arma;
using namespace std;
using namespace Rcpp;
 
double linCombMeancpp(SEXP matList, const arma::vec volFactor) {
  try{
    List matListRcpp(matList);

    NumericVector aRcpp = as<NumericVector>(matListRcpp["a"]);
    vec a = vec(aRcpp.begin(),aRcpp.size());
    
    NumericMatrix bRcpp = as<NumericMatrix>(matListRcpp["b"]);
    mat b = mat(bRcpp.begin(),bRcpp.nrow(),bRcpp.ncol());
    
    NumericVector linCoeffRcpp = as<NumericVector>(matListRcpp["linCoeff"]);
    vec linCoeff = vec(linCoeffRcpp.begin(),linCoeffRcpp.size());

    vec outA = a + b * volFactor;
    
    for (unsigned int ii=0; ii<outA.n_elem;ii++) {
      outA(ii) = expm1(outA(ii));
    }
    
    double mm = dot(linCoeff, outA);
    
    return(mm);
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}
