#include "meanVec.h"
#include "RcppArmadillo.h"
#include "linCombMean.h"

using namespace std;
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
arma::vec meanVecFun(List meanListS, const arma::vec& currVol) {
  try{
      List temp;
     
      // get the number of vol factors
      int N = meanListS.length();
     
      arma::vec meanVec(N);
      
      meanVec.zeros();
      
      List meanList(meanListS);
      
      // calculate predicted return means
      for (int uu = 0; uu<N; uu++) {
            temp = meanList[uu];
            meanVec(uu) = linCombMeancpp(temp,currVol);
      }
     
      return(meanVec);
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}
