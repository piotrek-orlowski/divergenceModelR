#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;

arma::mat divModelObsNoiseMat(const arma::vec& corrs, arma::vec& bpars, double spotVar, const arma::vec& matVec, const int U){
  
  // Initialize return matrix
  arma::mat res(3*U*matVec.n_elem,3*U*matVec.n_elem);
  
  // corrs is a 3x1 arma::vec and it serves to set up the correlation structure
  // among types of portfolios.
  arma::mat corrSmall = arma::eye(corrs.n_elem,corrs.n_elem);
  corrSmall(1,0) = corrs(0);
  corrSmall(2,0) = corrs(1);
  corrSmall(2,1) = corrs(2);
  corrSmall = arma::symmatl(corrSmall);
  arma::mat corrBig(res);
  
  arma::mat onesBig(U*matVec.n_elem,U*matVec.n_elem);
  onesBig = arma::kron(arma::eye(matVec.n_elem,matVec.n_elem), arma::ones<arma::mat>(U,U));
  // onesBig.fill(0.975);
  // onesBig.diag() = arma::ones(onesBig.n_cols);
  corrBig = arma::kron(corrSmall,onesBig);
  corrBig.diag(1) *= 0.95;
  corrBig.diag(-1) *= 0.95;
  
  // now set up the correlation structure between maturities using spotVar,
  // matVec and bpars. bpars is by type (length 3), we need a smart kronecker
  // product with matVec.
  
  spotVar = sqrt(spotVar);
  bpars *= spotVar;
  arma::mat btimeVar = arma::kron(bpars,arma::kron(arma::pow(matVec,0.25),arma::ones(U)));
  arma::mat btimeVarDiag = arma::diagmat(btimeVar);
  
  // res = corrBig;
  res = btimeVarDiag * corrBig * btimeVarDiag;
  return res;
}