#' @title Divergence pricing
#' @rdname divSwapRates
#' @description Calculation of prices of divergence swaps, divergence-based skewness and kurtosis swaps. All swaps are scale-free, that is they are swaps on the return. The \code{_base} functions calculate prices of components: unscaled divergence, unscaled \eqn{x^p \log{x}} divergence, unscaled \eqn{x^p \log^2{x}} divergence.
#' @param params.Q Q measure parameters
#' @param t.vec T-length vector of divergence swap maturities
#' @param vol.mat S x N.factors matrix of variance states
#' @param p divergence power. For skewness and kurtosis contracts, expansions will be taken around the Dp divergence.
#' @param ... further arguments to ODE solvers
#' @return S x T matrix of divergence swap prices
#' @export divergenceSwapRate

divergenceSwapRate <- function(p, params.Q, t.vec, vol.mat, jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), mod.type, ...){
  
  p0 <- p
  p <- matrix(p,ncol=1)
  p <- cbind(p,matrix(0,nrow=nrow(p), ncol = ncol(vol.mat)))
  if(any(p==0 | p == 1)){
    cfval.list <- affineCFderivs(u = p, params.Q = params.Q, params.P = NULL, t.vec = t.vec, v.0 = vol.mat, N.factors = dim(vol.mat)[2], jumpTransform = jumpTransform, mod.type = mod.type, ...)
    cfval <- cfval.list$cf
  } else {
    cfval <- affineCF(u = p, params.Q = params.Q, params.P = NULL, t.vec = t.vec, v.0 = vol.mat, jumpTransform = jumpTransform$TF, N.factors = dim(vol.mat)[2], CGF = FALSE, mod.type = mod.type, ...)
  }
  
  cfval <- cfval - 1  
  
  # p=0 and p=1 require special treatment: expectations of -Log(X) and X*Log(X)
  ind.0 <- which(p0==0)
  ind.1 <- which(p0==1)
  if(length(ind.0) > 0){
    # Derivative of CGF at 0
    cf.d0 <- - cfval.list$cf.d1[ind.0,,]
    cfval[ind.0,,] <- cf.d0
  }
  if(length(ind.1) > 0){
    # Derivative of CGF at 1
    cf.d0 <- - cfval.list$cf.d1[ind.1,,]
    cfval[ind.1,,] <- cf.d0
  }
  
  dimnames(cfval) <- list(paste0("p=",p0), paste0("t=",t.vec), paste0("V",1:nrow(vol.mat)))
  
  divSR <- cfval
  p0 <- 1/(p0*(p0-1))
  if(length(ind.0)>0){
    p0[ind.0] <- 1
  }
  if(length(ind.1)>0){
    p0[ind.1] <- -1
  }
  p0 <- array(p0,dim=dim(cfval))
  divSR <- p0*cfval
  
  return(divSR)
}

#' @export skewSwapRate
#' @describeIn divSwapRates

skewSwapRate <- function(p, params.Q, t.vec, vol.mat, jumpTransform, mod.type, ...){
  
  p0 <- p
  p <- matrix(p,ncol=1)
  p <- cbind(p,matrix(0,nrow=nrow(p), ncol = ncol(vol.mat)))
  
  ind.0 <- which(p0==0)
  ind.1 <- which(p0==1)
  
  # Two price components E(exp(py)) and E(y*exp(py)). Former first, simply CF
  cfval.list <- affineCFderivs(u = p, params.Q = params.Q, params.P = NULL, t.vec = t.vec, v.0 = vol.mat, N.factors = dim(vol.mat)[2], jumpTransform = jumpTransform, mod.type = mod.type, ...)
  cfval <- cfval.list$cf - 1
  
  # Latter, derivative at p
  cfval.d <- cfval.list$cf.d1
  
  # second derivative at p for 0 and 1
  cfval.d2 <- cfval.list$cf.d2
  # Handle 0 and 1
  
  # First and second derivative at 0
  if(length(ind.0) > 0){
    # Derivative of CGF at 0
    cfval[ind.0,,] <- -cfval.d[ind.0,,]
    # second derivative at 0
    cfval.d[ind.0,,] <- cfval.d2[ind.0,,]
  }
  # First and second derivative at 1
  if(length(ind.1) > 0){
    # Derivative of CGF at 1
    cfval[ind.1,,] <- cfval.d[ind.1,,]
    # Second derivative at 1
    cfval.d[ind.1,,] <- cfval.d2[ind.1,,]
  }
  
  pa <- 1/(p0*(p0-1))
  pa <- array(pa,dim=dim(cfval))
  pb <- -(2*p0-1)/(p0^2*(p0-1)^2)
  pb <- array(pb,dim = dim(cfval))
  if(length(ind.0)>0){
    pa[ind.0,,] <- -0.5
    pb[ind.0,,] <- 1
  }
  if(length(ind.1)>0){
    pa[ind.1,,] <- 0.5
    pb[ind.1,,] <- -1
  }
  divSR <- pa*cfval.d + pb*cfval
  
  dimnames(divSR) <- list(paste0("p=",p0), paste0("t=",t.vec), paste0("V",1:nrow(vol.mat)))
  return(divSR)
}

#' @describeIn divSwapRates
#' @export quartSwapRate

quartSwapRate <- function(p, params.Q, t.vec, vol.mat, jumpTransform, mod.type,  ...){
  
  p0 <- p
  p <- matrix(p,ncol=1)
  p <- cbind(p,matrix(0,nrow=nrow(p), ncol = ncol(vol.mat)))
  
  ind.0 <- which(p0==0)
  ind.1 <- which(p0==1)
  
  
  # Three price components E(exp(py)), E(y*exp(py)), E(y^2*exp(py)). Former first, simply CF
  cfval.list <- affineCFderivs(u = p, params.Q = params.Q, params.P = NULL, t.vec = t.vec, v.0 = vol.mat, N.factors = dim(vol.mat)[2], jumpTransform = jumpTransform, mod.type = mod.type, ...)
  cfval <- cfval.list$cf - 1
  
  # Latter, derivative at p
  cfval.d <- cfval.list$cf.d1
  
  # second derivative at p
  cfval.d2 <- cfval.list$cf.d2
  
  # third derivative at p
  cfval.d3 <- cfval.list$cf.d3
  
  # Handle 0 and 1
  # First and second derivative at 0
  if(length(ind.0) > 0){
    # Derivative of CGF at 0
    cfval[ind.0,,] <- cfval.d[ind.0,,]
    # second derivative at 0
    cfval.d[ind.0,,] <- cfval.d2[ind.0,,]
    # third derivative at 0
    cfval.d2[ind.0,,] <- cfval.d3[ind.0,,]
  }
  # First and second derivative at 1
  if(length(ind.1) > 0){
    # Derivative of CGF at 1
    cfval[ind.1,,] <- cfval.d[ind.1,,]
    # Second derivative at 1
    cfval.d[ind.1,,] <- cfval.d2[ind.1,,]
    # third derivative at 1
    cfval.d2[ind.1,,] <- cfval.d3[ind.1,,]
  }
  
  pa <- (p0^2-2*p0^3+p0^4)/(p0^3*(p0-1)^3)
  pa <- array(pa,dim=dim(cfval))
  pb <- (-2*p0+6*p0^2-4*p0^3)/(p0^3*(p0-1)^3)
  pb <- array(pb,dim = dim(cfval))
  pc <- (2-6*p0+6*p0^2)/(p0^3*(p0-1)^3)
  pc <- array(pc,dim = dim(cfval))
  if(length(ind.0)>0){
    pa[ind.0,,] <- -1/3
    pb[ind.0,,] <- -1
    pc[ind.0,,] <- -2
  }
  if(length(ind.1)>0){
    pa[ind.1,,] <- 1/3
    pb[ind.1,,] <- -1
    pc[ind.1,,] <- 2
  }
  divSR <- pa*cfval.d2 + pb*cfval.d + pc*cfval
  
  dimnames(divSR) <- list(paste0("p=",p0), paste0("t=",t.vec), paste0("V",1:nrow(vol.mat)))
  return(divSR)
}