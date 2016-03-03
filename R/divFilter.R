#' @title Divergence filtering
#' @rdname divFilter
#' @description Based on specification, this function calculates (1) matrices for model dynamics calculations, (2) matrices for CF and derivatives calculations. These are used as inputs in a c++ filtering function.
#' @param params.P P measure parameters
#' @param params.Q Q measure parameters
#' @param data
#' @param spec \code{data.frame} with fields governing maturity (t), divergence power (p), divergence type (character)
#' @param N.factors \code{numeric}, number of SV factors
#' @param ... further arguments to ODE solvers
#' @return S x T matrix of divergence swap prices
#' @export divergenceFilter

divergenceFilter <- function(params.P, params.Q, data, spec, N.factors, tfPointers = getPointerToJumpTransform("kouExpJumpTransform"), ...){
  
  # spec has to be ordered by t, p, type
  spec.check <- spec[order(spec$t, spec$p, spec$type),]
  stopifnot(identical(spec,spec.check))
  
  # Translate spec into mkt
  mkt.spec <- data.frame(p = 1, q = 0, r = 0, t = unique(spec$t))
  p.spec <- matrix(0, nrow = length(unique(spec$p)), ncol = N.factors+1)
  p.spec[,1] <- unique(spec$p)
  
  # Solve ODEs for pricing (parameters for observation equations)
  pricing.sol <- odeExtSolveWrap(u = p.spec, params.Q = params.Q, params.P = params.P, mkt = mkt.spec, rtol = 1e-6, atol = 1e-10, N.factors = N.factors, jumpTransform = tfPointers, mod.type = 'standard', ...)
  
  # Solve for stock return and vol factors conditional moments (parameters for state dynamics)
  dynamics.sol <- modelDynamics(u.t.mat = data.frame(u=1,t=1), params.P = params.P, params.Q = params.Q, jumpTransform = tfPointers$TF, dT = 5/252, transform.matrix = TRUE, N.factos = N.factos, N.points = 2, ...)
  
  # Call CPP
  
}