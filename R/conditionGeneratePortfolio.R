#' Prepare elements for hedged GLT return moment calculation
#' @description Generate list of data frames that describes a sum of hedged GLT portfolio returns or their squares.
#' @param u.vec The frequency vector for which the portfolio is built. Each frequency represents a cross-product of contemporeneous hedged returns.
#' @param tau.js Vector of horizons for which the return is calculated
#' @param deltas Vector of hedging frequencies; if \code{length(detlas) == 1}, \code{rep} is applied to have length equal to \code{length(tau.js)}
#' @param ttm Vector of initial maturity of the portfolio.
#' @param weights A vector of weights to pre-multiply the returns (useful for setting up quadrature approximations)
#' @export
#' @return return.df.list List with two fields: cond.list describing the moment conditions and cond.coeff with coefficients for linear combination.
#'
condGeneratePortfolio <- function(u.vec, tau.js = 1/2 * 5/252, deltas = 1/252/6.5, ttm = 1/12, weights = rep(1,length(tau.js)), N.factors = 3){
  
  # check that tau.js and deltas aren't too big
  stopifnot(max(tau.js) + max(deltas) < ttm)
  
  # make sure we have the right number of u-s
  stopifnot(length(u.vec) == length(ttm))
  
  stopifnot(length(weights) == length(tau.js))
  
  # check if deltas is a single value or more
  if(length(deltas) > 1){
    stopifnot(length(deltas) == length(tau.js))
  }
  else{
    deltas <- rep(deltas, length(tau.js))
  }
  
  return.df.list <- list(cond.list = list(), cond.coeff = numeric(), linTrans = NULL)
  
  # start generating the stuff you need, and keep multiplying each cond list by itself as many times as you have powers.
  for(nn in 1:length(tau.js)){
    loc.list <- hedged.return(u.vec[1],tau.j = tau.js[nn], delta = deltas[nn], ttm = ttm[1], N.factors = N.factors)
    n.terms <- 1
    while(n.terms < length(u.vec)){
      loc.mult <- hedged.return(u.vec[n.terms+1],tau.j = tau.js[nn], delta = deltas[nn], ttm = ttm[n.terms+1], N.factors = N.factors)
      loc.list <- condMultiply(loc.list,loc.mult)
      n.terms <- n.terms+1
    }
    # add term to the condition list
    return.df.list <- condAdd(return.df.list,loc.list,weights=c(1,weights[nn]))
  }
  
  # now a sanity check to see that we are indeed calculating the right combinations
  stopifnot(return.df.list$linTrans[1,] == return.df.list$cond.coeff)
  return(return.df.list)
}