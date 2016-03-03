#' @title VIX index and VIX FUTURES pricing
#' @name vixfut
#' @description Calculate the VIX index and the price of VIX futures contracts with arbitrary maturities in affine jump diffusion models.
#' @param params.Q Q measure parameters.
#' @param t.vec maturity of VIX index, length \code{T}.
#' @param t maturity of futures contract, length 1.
#' @param tau maturity of VIX index, on which the futures are written, length 1.
#' @param v.0 initial variance factor value, size \code{S x (N.factors+1)} for \code{priceVIX}, size \code{1 x (N.factors+1)} for \code{priceVIXFUT}.
#' @param ... arguments to ODE solution routines
#' @details \code{priceVIX} calculates the VIX index for an arbitrary set of starting factor values, dates and maturities. \code{priceVIXFUT} only works for a single index time to maturity and a single futures time to maturity. Time to maturity of the index, on which the futures are written, is given net of the futures maturity, i.e. for calculating the price of a futures contract, that matures in 6 months, on the index value between 6 months and 7 months from now, give \code{t = 6/12} and \code{tau = 1/12}.
#' @export priceVIX
#' @return VIX futures value in 1-factor model with v.0 starting values, VIX value

priceVIX <- function(v.0, t.vec, params.Q, ...){
  
  N.factors <- ncol(v.0)
  u.0 <- matrix(rep(0,N.factors+1),nrow=1)
  cf.vals <- affineCFderivs(u = u.0, params.Q = params.Q, params.P = NULL, t.vec = t.vec, v.0 = v.0, N.factors = N.factors, ...)
  cf.d1 <- drop(cf.vals$cf.d1)
  t.mat <- matrix(rep(t.vec,nrow(v.0)), nrow = length(t.vec), ncol = nrow(v.0))
  
  VIX <- sqrt(-2/t.mat * cf.d1)
  return(VIX)
}

#' @export priceVIXFUT
#' @describeIn vixfut

priceVIXFUT <- function(v.0, t, tau, params.Q, ...){
  
  N.factors <- ncol(v.0)
  mkt <- data.frame(p= 1, r= 0, q= 0, t= tau)
  odestrucs <- ODEstructs(params = params.Q, jumpTransform = expNormJumpTransform, mkt = mkt, N.factors = N.factors)
  
  # Affine coeff derivs: these lads are for the `vix' in `vix futures'
  ode.sol <- odeExtSolveWrap(u = matrix(0,nrow = 1, ncol = N.factors+1), params.Q = params.Q, params.P = NULL, mkt = mkt, N.factors = N.factors,...)
  
  ode.sol <- -2/tau*ode.sol[,,c(paste0("bp1",1:N.factors),"ap1")]
  
  mkt.futures <- data.frame(p=1,q=0,r=0,t=t)
  
  # Funtion to be integrated: the characteristic function at a specific point, multiplied by a magic number
  integrandFooVec <- function(zz){
    # need derivatives of affine coefficients with respect to V args
    # start with wrapper
    cfDerivFoo <- function(zx){
      u.loc <- apply(matrix(zx),1,function(z)return(ode.sol[paste0("bp1",1:N.factors)] * z))
      u.loc <- cbind(0,u.loc)
      loc.cf <- affineCF(u = u.loc, params.Q = params.Q, params.P = NULL, t.vec = mkt.futures$t, v.0 = v.0, N.factors = N.factors, ...)
      return(drop(Re(loc.cf)))
    }
    hh <- 1e-4
    uu <- seq(-5*hh,5*hh,by=hh)
    zz.for.diff <- matrix(zz,nrow=length(zz),ncol=length(uu))
    zz.for.diff <- zz.for.diff + matrix(uu, ncol = length(uu), nrow = length(zz), byrow = TRUE)
    cf.for.diff <- matrix(cfDerivFoo(as.numeric(zz.for.diff)), ncol = ncol(zz.for.diff), nrow = nrow(zz.for.diff))
    cf.d1 <- apply(X = cf.for.diff, MARGIN = 1, FUN = function(cff){
      f0 <- cff[median(order(cff))]
      fVal <- cff[-median(order(cff))]
      fVal <- fVal[-c(1,10)]
      deriv <- fdDerivative(f0 = f0, fVal = fVal, p = 1, k = 7, h = hh)
      return(deriv)
    })
    res.A <- exp(ode.sol["ap1"] * zz)
    res.A.deriv <- ode.sol["ap1"]*res.A
    res <- res.A.deriv * cf.for.diff[,which(uu==0)]
    res <- res + res.A * cf.d1
    res <- res * (-zz)^(-0.5)
    res <- Re(res)
    return(res)
  }
  
  vf <- 1/gamma(0.5) * integrate(f = integrandFooVec, lower = -Inf, upper = 0)$value
  return(vf)
}
