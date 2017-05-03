#include <RcppArmadillo.h>

using namespace std;

arma::mat affineAdditiveStateController(arma::mat stateMat){
  
  // This function sets all non-positive state values to a small number
  arma::uvec nonPosElements = arma::find(stateMat <= 0.0);
  
  // set
  if(nonPosElements.n_elem > 0){
    stateMat.elem(nonPosElements).fill(1e-6);
  }
  
  return stateMat;
}