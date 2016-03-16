#include "covMat.h"
#include "RcppArmadillo.h"
#include "linCombMean.h"

using namespace std;
using namespace Rcpp;

//' @export
// [[Rcpp::export]] 
arma::mat covMatFun(List covListS, const arma::uvec covListDim, const arma::vec currVol) {
  try{
      List temp;
      // get the output covariance dimension
      // int N = currVol.n_elem;
      unsigned int N = covListDim(0);
      unsigned int M = covListDim(1);
      
      arma::mat covMat(N,M,arma::fill::zeros);

      // calculate covariance matrix, which will be used to calculate C and also the residual R
      covMat.zeros();
      // calculate predicted return innovation matrix
      unsigned int uInd = 0;
      
      List covList(covListS);

      for (unsigned int uu1 = 0; uu1<M; uu1++) {
        for (unsigned int uu2 =uu1; uu2<N; uu2++) {
            temp = covList[uInd];
            covMat(uu2,uu1) = linCombMeancpp(temp,currVol);
            uInd++;
        }
      }
      // make the matrix symmetric if need be
      if(N == M){
        covMat = arma::symmatl(covMat); 
      }
      
    return(covMat);
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}
