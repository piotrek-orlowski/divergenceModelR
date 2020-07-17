// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/divergenceModelR.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// affineContract_sqrtFilter
List affineContract_sqrtFilter(const arma::mat dataMat, const arma::vec initState, const arma::mat initProcCov, const List modelParams);
RcppExport SEXP _divergenceModelR_affineContract_sqrtFilter(SEXP dataMatSEXP, SEXP initStateSEXP, SEXP initProcCovSEXP, SEXP modelParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dataMat(dataMatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type initState(initStateSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type initProcCov(initProcCovSEXP);
    Rcpp::traits::input_parameter< const List >::type modelParams(modelParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(affineContract_sqrtFilter(dataMat, initState, initProcCov, modelParams));
    return rcpp_result_gen;
END_RCPP
}
// portfolio_sqrtFilter
List portfolio_sqrtFilter(const arma::mat dataMat, const arma::vec initState, const arma::mat initProcCov, const List modelParams);
RcppExport SEXP _divergenceModelR_portfolio_sqrtFilter(SEXP dataMatSEXP, SEXP initStateSEXP, SEXP initProcCovSEXP, SEXP modelParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dataMat(dataMatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type initState(initStateSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type initProcCov(initProcCovSEXP);
    Rcpp::traits::input_parameter< const List >::type modelParams(modelParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(portfolio_sqrtFilter(dataMat, initState, initProcCov, modelParams));
    return rcpp_result_gen;
END_RCPP
}
// portfolio_noStock_sqrtFilter
List portfolio_noStock_sqrtFilter(const arma::mat dataMat, const arma::vec initState, const arma::mat initProcCov, const List modelParams);
RcppExport SEXP _divergenceModelR_portfolio_noStock_sqrtFilter(SEXP dataMatSEXP, SEXP initStateSEXP, SEXP initProcCovSEXP, SEXP modelParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dataMat(dataMatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type initState(initStateSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type initProcCov(initProcCovSEXP);
    Rcpp::traits::input_parameter< const List >::type modelParams(modelParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(portfolio_noStock_sqrtFilter(dataMat, initState, initProcCov, modelParams));
    return rcpp_result_gen;
END_RCPP
}
// portfolio_noStock_simpleFilter
List portfolio_noStock_simpleFilter(const arma::mat dataMat, const arma::vec initState, const arma::mat initProcCov, const List modelParams);
RcppExport SEXP _divergenceModelR_portfolio_noStock_simpleFilter(SEXP dataMatSEXP, SEXP initStateSEXP, SEXP initProcCovSEXP, SEXP modelParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dataMat(dataMatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type initState(initStateSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type initProcCov(initProcCovSEXP);
    Rcpp::traits::input_parameter< const List >::type modelParams(modelParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(portfolio_noStock_simpleFilter(dataMat, initState, initProcCov, modelParams));
    return rcpp_result_gen;
END_RCPP
}
// DSQ_sqrtFilter
List DSQ_sqrtFilter(const arma::mat dataMat, const arma::vec initState, const arma::mat initProcCov, const List modelParams);
RcppExport SEXP _divergenceModelR_DSQ_sqrtFilter(SEXP dataMatSEXP, SEXP initStateSEXP, SEXP initProcCovSEXP, SEXP modelParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dataMat(dataMatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type initState(initStateSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type initProcCov(initProcCovSEXP);
    Rcpp::traits::input_parameter< const List >::type modelParams(modelParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(DSQ_sqrtFilter(dataMat, initState, initProcCov, modelParams));
    return rcpp_result_gen;
END_RCPP
}
// D_sqrtFilter
List D_sqrtFilter(const arma::mat dataMat, const arma::vec initState, const arma::mat initProcCov, const List modelParams);
RcppExport SEXP _divergenceModelR_D_sqrtFilter(SEXP dataMatSEXP, SEXP initStateSEXP, SEXP initProcCovSEXP, SEXP modelParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dataMat(dataMatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type initState(initStateSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type initProcCov(initProcCovSEXP);
    Rcpp::traits::input_parameter< const List >::type modelParams(modelParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(D_sqrtFilter(dataMat, initState, initProcCov, modelParams));
    return rcpp_result_gen;
END_RCPP
}
// DS_sqrtFilter
List DS_sqrtFilter(const arma::mat dataMat, const arma::vec initState, const arma::mat initProcCov, const List modelParams);
RcppExport SEXP _divergenceModelR_DS_sqrtFilter(SEXP dataMatSEXP, SEXP initStateSEXP, SEXP initProcCovSEXP, SEXP modelParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dataMat(dataMatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type initState(initStateSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type initProcCov(initProcCovSEXP);
    Rcpp::traits::input_parameter< const List >::type modelParams(modelParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(DS_sqrtFilter(dataMat, initState, initProcCov, modelParams));
    return rcpp_result_gen;
END_RCPP
}
// affineObservationStateHandler
Rcpp::List affineObservationStateHandler(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);
RcppExport SEXP _divergenceModelR_affineObservationStateHandler(SEXP stateMatSEXP, SEXP modelParametersSEXP, SEXP iterCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type stateMat(stateMatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type modelParameters(modelParametersSEXP);
    Rcpp::traits::input_parameter< const int >::type iterCount(iterCountSEXP);
    rcpp_result_gen = Rcpp::wrap(affineObservationStateHandler(stateMat, modelParameters, iterCount));
    return rcpp_result_gen;
END_RCPP
}
// affineObservationStateHandler_D
Rcpp::List affineObservationStateHandler_D(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);
RcppExport SEXP _divergenceModelR_affineObservationStateHandler_D(SEXP stateMatSEXP, SEXP modelParametersSEXP, SEXP iterCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type stateMat(stateMatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type modelParameters(modelParametersSEXP);
    Rcpp::traits::input_parameter< const int >::type iterCount(iterCountSEXP);
    rcpp_result_gen = Rcpp::wrap(affineObservationStateHandler_D(stateMat, modelParameters, iterCount));
    return rcpp_result_gen;
END_RCPP
}
// affineObservationStateHandler_DS
Rcpp::List affineObservationStateHandler_DS(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);
RcppExport SEXP _divergenceModelR_affineObservationStateHandler_DS(SEXP stateMatSEXP, SEXP modelParametersSEXP, SEXP iterCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type stateMat(stateMatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type modelParameters(modelParametersSEXP);
    Rcpp::traits::input_parameter< const int >::type iterCount(iterCountSEXP);
    rcpp_result_gen = Rcpp::wrap(affineObservationStateHandler_DS(stateMat, modelParameters, iterCount));
    return rcpp_result_gen;
END_RCPP
}
// affineObservationStateHandler_affineContracts
Rcpp::List affineObservationStateHandler_affineContracts(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);
RcppExport SEXP _divergenceModelR_affineObservationStateHandler_affineContracts(SEXP stateMatSEXP, SEXP modelParametersSEXP, SEXP iterCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type stateMat(stateMatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type modelParameters(modelParametersSEXP);
    Rcpp::traits::input_parameter< const int >::type iterCount(iterCountSEXP);
    rcpp_result_gen = Rcpp::wrap(affineObservationStateHandler_affineContracts(stateMat, modelParameters, iterCount));
    return rcpp_result_gen;
END_RCPP
}
// affineObservationStateHandler_cumulant
Rcpp::List affineObservationStateHandler_cumulant(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);
RcppExport SEXP _divergenceModelR_affineObservationStateHandler_cumulant(SEXP stateMatSEXP, SEXP modelParametersSEXP, SEXP iterCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type stateMat(stateMatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type modelParameters(modelParametersSEXP);
    Rcpp::traits::input_parameter< const int >::type iterCount(iterCountSEXP);
    rcpp_result_gen = Rcpp::wrap(affineObservationStateHandler_cumulant(stateMat, modelParameters, iterCount));
    return rcpp_result_gen;
END_RCPP
}
// affineObservationStateHandler_cumulant_D
Rcpp::List affineObservationStateHandler_cumulant_D(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);
RcppExport SEXP _divergenceModelR_affineObservationStateHandler_cumulant_D(SEXP stateMatSEXP, SEXP modelParametersSEXP, SEXP iterCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type stateMat(stateMatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type modelParameters(modelParametersSEXP);
    Rcpp::traits::input_parameter< const int >::type iterCount(iterCountSEXP);
    rcpp_result_gen = Rcpp::wrap(affineObservationStateHandler_cumulant_D(stateMat, modelParameters, iterCount));
    return rcpp_result_gen;
END_RCPP
}
// affineObservationStateHandler_cumulant_DS
Rcpp::List affineObservationStateHandler_cumulant_DS(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);
RcppExport SEXP _divergenceModelR_affineObservationStateHandler_cumulant_DS(SEXP stateMatSEXP, SEXP modelParametersSEXP, SEXP iterCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type stateMat(stateMatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type modelParameters(modelParametersSEXP);
    Rcpp::traits::input_parameter< const int >::type iterCount(iterCountSEXP);
    rcpp_result_gen = Rcpp::wrap(affineObservationStateHandler_cumulant_DS(stateMat, modelParameters, iterCount));
    return rcpp_result_gen;
END_RCPP
}
// affineObservationStateHandler_optionPortfolios
Rcpp::List affineObservationStateHandler_optionPortfolios(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);
RcppExport SEXP _divergenceModelR_affineObservationStateHandler_optionPortfolios(SEXP stateMatSEXP, SEXP modelParametersSEXP, SEXP iterCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type stateMat(stateMatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type modelParameters(modelParametersSEXP);
    Rcpp::traits::input_parameter< const int >::type iterCount(iterCountSEXP);
    rcpp_result_gen = Rcpp::wrap(affineObservationStateHandler_optionPortfolios(stateMat, modelParameters, iterCount));
    return rcpp_result_gen;
END_RCPP
}
// affineObservationStateHandler_optionPortfolios_noStock
Rcpp::List affineObservationStateHandler_optionPortfolios_noStock(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);
RcppExport SEXP _divergenceModelR_affineObservationStateHandler_optionPortfolios_noStock(SEXP stateMatSEXP, SEXP modelParametersSEXP, SEXP iterCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type stateMat(stateMatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type modelParameters(modelParametersSEXP);
    Rcpp::traits::input_parameter< const int >::type iterCount(iterCountSEXP);
    rcpp_result_gen = Rcpp::wrap(affineObservationStateHandler_optionPortfolios_noStock(stateMat, modelParameters, iterCount));
    return rcpp_result_gen;
END_RCPP
}
// affineTransitionStateHandler
Rcpp::List affineTransitionStateHandler(const arma::mat& stateMat, const Rcpp::List& modelParameters, const int iterCount);
RcppExport SEXP _divergenceModelR_affineTransitionStateHandler(SEXP stateMatSEXP, SEXP modelParametersSEXP, SEXP iterCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type stateMat(stateMatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type modelParameters(modelParametersSEXP);
    Rcpp::traits::input_parameter< const int >::type iterCount(iterCountSEXP);
    rcpp_result_gen = Rcpp::wrap(affineTransitionStateHandler(stateMat, modelParameters, iterCount));
    return rcpp_result_gen;
END_RCPP
}
// covMatFun
arma::mat covMatFun(List covListS, const arma::uvec covListDim, const arma::vec currVol);
RcppExport SEXP _divergenceModelR_covMatFun(SEXP covListSSEXP, SEXP covListDimSEXP, SEXP currVolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type covListS(covListSSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type covListDim(covListDimSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type currVol(currVolSEXP);
    rcpp_result_gen = Rcpp::wrap(covMatFun(covListS, covListDim, currVol));
    return rcpp_result_gen;
END_RCPP
}
// cumulant_DS_sqrtFilter
List cumulant_DS_sqrtFilter(const arma::mat dataMat, const arma::vec initState, const arma::mat initProcCov, const List modelParams);
RcppExport SEXP _divergenceModelR_cumulant_DS_sqrtFilter(SEXP dataMatSEXP, SEXP initStateSEXP, SEXP initProcCovSEXP, SEXP modelParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dataMat(dataMatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type initState(initStateSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type initProcCov(initProcCovSEXP);
    Rcpp::traits::input_parameter< const List >::type modelParams(modelParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(cumulant_DS_sqrtFilter(dataMat, initState, initProcCov, modelParams));
    return rcpp_result_gen;
END_RCPP
}
// divModelObsNoiseMat
arma::mat divModelObsNoiseMat(const arma::vec corrs, arma::vec bpars, double spotVar, const arma::vec matVec, const int U);
RcppExport SEXP _divergenceModelR_divModelObsNoiseMat(SEXP corrsSEXP, SEXP bparsSEXP, SEXP spotVarSEXP, SEXP matVecSEXP, SEXP USEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type corrs(corrsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bpars(bparsSEXP);
    Rcpp::traits::input_parameter< double >::type spotVar(spotVarSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type matVec(matVecSEXP);
    Rcpp::traits::input_parameter< const int >::type U(USEXP);
    rcpp_result_gen = Rcpp::wrap(divModelObsNoiseMat(corrs, bpars, spotVar, matVec, U));
    return rcpp_result_gen;
END_RCPP
}
// divergenceSwapRateCpp
arma::cube divergenceSwapRateCpp(const arma::vec p, const arma::cube coeffs, const arma::mat stateMat);
RcppExport SEXP _divergenceModelR_divergenceSwapRateCpp(SEXP pSEXP, SEXP coeffsSEXP, SEXP stateMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type coeffs(coeffsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type stateMat(stateMatSEXP);
    rcpp_result_gen = Rcpp::wrap(divergenceSwapRateCpp(p, coeffs, stateMat));
    return rcpp_result_gen;
END_RCPP
}
// skewnessSwapRateCpp
arma::cube skewnessSwapRateCpp(const arma::vec p, const arma::cube coeffs, const arma::mat stateMat);
RcppExport SEXP _divergenceModelR_skewnessSwapRateCpp(SEXP pSEXP, SEXP coeffsSEXP, SEXP stateMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type coeffs(coeffsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type stateMat(stateMatSEXP);
    rcpp_result_gen = Rcpp::wrap(skewnessSwapRateCpp(p, coeffs, stateMat));
    return rcpp_result_gen;
END_RCPP
}
// quarticitySwapRateCpp
arma::cube quarticitySwapRateCpp(const arma::vec p, const arma::cube coeffs, const arma::mat stateMat);
RcppExport SEXP _divergenceModelR_quarticitySwapRateCpp(SEXP pSEXP, SEXP coeffsSEXP, SEXP stateMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type coeffs(coeffsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type stateMat(stateMatSEXP);
    rcpp_result_gen = Rcpp::wrap(quarticitySwapRateCpp(p, coeffs, stateMat));
    return rcpp_result_gen;
END_RCPP
}
// glPricer_cpp
arma::cube glPricer_cpp(const arma::cube& strikeMat, const arma::mat& mkt, const arma::vec& glWts, const arma::vec& glNodes, arma::cx_cube& cfVals, int Nfactors, double alpha, double sigmaRef);
RcppExport SEXP _divergenceModelR_glPricer_cpp(SEXP strikeMatSEXP, SEXP mktSEXP, SEXP glWtsSEXP, SEXP glNodesSEXP, SEXP cfValsSEXP, SEXP NfactorsSEXP, SEXP alphaSEXP, SEXP sigmaRefSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type strikeMat(strikeMatSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mkt(mktSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type glWts(glWtsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type glNodes(glNodesSEXP);
    Rcpp::traits::input_parameter< arma::cx_cube& >::type cfVals(cfValsSEXP);
    Rcpp::traits::input_parameter< int >::type Nfactors(NfactorsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sigmaRef(sigmaRefSEXP);
    rcpp_result_gen = Rcpp::wrap(glPricer_cpp(strikeMat, mkt, glWts, glNodes, cfVals, Nfactors, alpha, sigmaRef));
    return rcpp_result_gen;
END_RCPP
}
// meanVecFun
arma::vec meanVecFun(List meanListS, const arma::vec currVol);
RcppExport SEXP _divergenceModelR_meanVecFun(SEXP meanListSSEXP, SEXP currVolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type meanListS(meanListSSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type currVol(currVolSEXP);
    rcpp_result_gen = Rcpp::wrap(meanVecFun(meanListS, currVol));
    return rcpp_result_gen;
END_RCPP
}
// testAffineFilter
List testAffineFilter(const arma::mat testDataMat, const arma::vec testInitState, const arma::mat testInitProcCov, const List testModelParams);
RcppExport SEXP _divergenceModelR_testAffineFilter(SEXP testDataMatSEXP, SEXP testInitStateSEXP, SEXP testInitProcCovSEXP, SEXP testModelParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type testDataMat(testDataMatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type testInitState(testInitStateSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type testInitProcCov(testInitProcCovSEXP);
    Rcpp::traits::input_parameter< const List >::type testModelParams(testModelParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(testAffineFilter(testDataMat, testInitState, testInitProcCov, testModelParams));
    return rcpp_result_gen;
END_RCPP
}
// testSqrtAffineFilter
List testSqrtAffineFilter(const arma::mat testDataMat, const arma::vec testInitState, const arma::mat testInitProcCov, const List testModelParams);
RcppExport SEXP _divergenceModelR_testSqrtAffineFilter(SEXP testDataMatSEXP, SEXP testInitStateSEXP, SEXP testInitProcCovSEXP, SEXP testModelParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type testDataMat(testDataMatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type testInitState(testInitStateSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type testInitProcCov(testInitProcCovSEXP);
    Rcpp::traits::input_parameter< const List >::type testModelParams(testModelParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(testSqrtAffineFilter(testDataMat, testInitState, testInitProcCov, testModelParams));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP cirIntExp(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP luComplex(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_divergenceModelR_affineContract_sqrtFilter", (DL_FUNC) &_divergenceModelR_affineContract_sqrtFilter, 4},
    {"_divergenceModelR_portfolio_sqrtFilter", (DL_FUNC) &_divergenceModelR_portfolio_sqrtFilter, 4},
    {"_divergenceModelR_portfolio_noStock_sqrtFilter", (DL_FUNC) &_divergenceModelR_portfolio_noStock_sqrtFilter, 4},
    {"_divergenceModelR_portfolio_noStock_simpleFilter", (DL_FUNC) &_divergenceModelR_portfolio_noStock_simpleFilter, 4},
    {"_divergenceModelR_DSQ_sqrtFilter", (DL_FUNC) &_divergenceModelR_DSQ_sqrtFilter, 4},
    {"_divergenceModelR_D_sqrtFilter", (DL_FUNC) &_divergenceModelR_D_sqrtFilter, 4},
    {"_divergenceModelR_DS_sqrtFilter", (DL_FUNC) &_divergenceModelR_DS_sqrtFilter, 4},
    {"_divergenceModelR_affineObservationStateHandler", (DL_FUNC) &_divergenceModelR_affineObservationStateHandler, 3},
    {"_divergenceModelR_affineObservationStateHandler_D", (DL_FUNC) &_divergenceModelR_affineObservationStateHandler_D, 3},
    {"_divergenceModelR_affineObservationStateHandler_DS", (DL_FUNC) &_divergenceModelR_affineObservationStateHandler_DS, 3},
    {"_divergenceModelR_affineObservationStateHandler_affineContracts", (DL_FUNC) &_divergenceModelR_affineObservationStateHandler_affineContracts, 3},
    {"_divergenceModelR_affineObservationStateHandler_cumulant", (DL_FUNC) &_divergenceModelR_affineObservationStateHandler_cumulant, 3},
    {"_divergenceModelR_affineObservationStateHandler_cumulant_D", (DL_FUNC) &_divergenceModelR_affineObservationStateHandler_cumulant_D, 3},
    {"_divergenceModelR_affineObservationStateHandler_cumulant_DS", (DL_FUNC) &_divergenceModelR_affineObservationStateHandler_cumulant_DS, 3},
    {"_divergenceModelR_affineObservationStateHandler_optionPortfolios", (DL_FUNC) &_divergenceModelR_affineObservationStateHandler_optionPortfolios, 3},
    {"_divergenceModelR_affineObservationStateHandler_optionPortfolios_noStock", (DL_FUNC) &_divergenceModelR_affineObservationStateHandler_optionPortfolios_noStock, 3},
    {"_divergenceModelR_affineTransitionStateHandler", (DL_FUNC) &_divergenceModelR_affineTransitionStateHandler, 3},
    {"_divergenceModelR_covMatFun", (DL_FUNC) &_divergenceModelR_covMatFun, 3},
    {"_divergenceModelR_cumulant_DS_sqrtFilter", (DL_FUNC) &_divergenceModelR_cumulant_DS_sqrtFilter, 4},
    {"_divergenceModelR_divModelObsNoiseMat", (DL_FUNC) &_divergenceModelR_divModelObsNoiseMat, 5},
    {"_divergenceModelR_divergenceSwapRateCpp", (DL_FUNC) &_divergenceModelR_divergenceSwapRateCpp, 3},
    {"_divergenceModelR_skewnessSwapRateCpp", (DL_FUNC) &_divergenceModelR_skewnessSwapRateCpp, 3},
    {"_divergenceModelR_quarticitySwapRateCpp", (DL_FUNC) &_divergenceModelR_quarticitySwapRateCpp, 3},
    {"_divergenceModelR_glPricer_cpp", (DL_FUNC) &_divergenceModelR_glPricer_cpp, 8},
    {"_divergenceModelR_meanVecFun", (DL_FUNC) &_divergenceModelR_meanVecFun, 2},
    {"_divergenceModelR_testAffineFilter", (DL_FUNC) &_divergenceModelR_testAffineFilter, 4},
    {"_divergenceModelR_testSqrtAffineFilter", (DL_FUNC) &_divergenceModelR_testSqrtAffineFilter, 4},
    {"cirIntExp", (DL_FUNC) &cirIntExp, 4},
    {"luComplex", (DL_FUNC) &luComplex, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_divergenceModelR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
