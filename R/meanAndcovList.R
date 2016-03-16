#' @title Jump-diffusion model dynamics
#' @description  This function creates mean and covariance structures that can be passed to linCombMean to combine into an estimate of the conditional first and second (cross-)moments of the stock return and volatility factors.
#' @param params.P
#' @param params.Q
#' @param dT The Time horizon
#' @param rtol The relative tolerance required in the integration
#' @param N.factors The number of stochastic vol factors assumed.
#' @param N.points Number of points in sparse grid quadratures. 2 usually works pretty good, more than 3 is overkill.
#' @param jumpTransform pointer to jump trasform c++ function
#' @param mod.type character, \code{'standard'} or \code{'casc.vol'}
#' @export
#' @return Returns 2 lists and a matrix. 
modelDynamics <- function(params.P, params.Q, dT=5/252, rtol=1e-14, N.factors = 3, N.points= 2, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, mod.type = 'standard', ...) {
  
  hedge.freq <- 1/78/252
  
  # first calculate mean propagation lists
  stock.struct <- condStockChange(ret.period=dT)
  mean.lists <- list(momentCondition(params.P,params.Q,jumpTransform = jumpTransform,condition.struct=stock.struct,conditional=TRUE,rtol=rtol,...))

  # create function list. First is stock, then U portfolios and then the vol factors
  list.fun <- list(function(x) condStockChange(ret.period=hedge.freq, ret.start = x,N.factors=N.factors))
  
  # now add vol factors
  for (vv in 1:N.factors) {
    list.fun <- append(list.fun,values=function(x) condVolChange(factorNum=vv,ret.start=x,ret.period=hedge.freq,N.factors=N.factors))
    environment(list.fun[[length(list.fun)]]) <- new.env()
    assign("vv",vv,environment(list.fun[[length(list.fun)]]))
  }
  
  # calculate quadrature nodes, using sparse nested grids
  quad.1d <- createSparseGrid(type="KPU",dimension=1,k=N.points)
  quad.2d <- createSparseGrid(type="KPU",dimension=2,k=N.points)
  
  # initialize lists that will hold the structures
  mean.vec <- array(list(),c(1+N.factors))
  
  # we need to seperately do the 2d (autocovariance terms) and 1d quad (variance terms)
  cov.array <- array(list(),c(rep(1+N.factors,2)))
  
  # calculate mean list
  for (nn in 1:length(list.fun)) {
    m.list <- list(a=NULL,b=NULL,linCoeff=NULL)
    for (qq in 1:length(quad.1d$nodes)) {
      m1 <- list.fun[[nn]](dT * quad.1d$nodes[qq])
      ab <- momentCondition(params.P,params.Q,jumpTransform = jumpTransform,condition.struct=m1,conditional=TRUE,rtol=rtol,N.factors=N.factors,...)
      m.list$a <- c(m.list$a,ab$a)
      m.list$b <- rbind(m.list$b,ab$b)
      m.list$linCoeff <- c(m.list$linCoeff,ab$linCoeff * quad.1d$weights[qq] * dT / hedge.freq)
    }
    mean.vec[[nn]] <- m.list
  }
  
  # calculate covariance list
  for (nn1 in 1:length(list.fun)) {
    for (nn2 in 1:nn1) {
      m.list <- list(a=NULL,b=NULL,linCoeff=NULL)
      for (qq in 1:length(quad.2d$weights)) {
        m1 <- list.fun[[nn1]](dT * quad.2d$nodes[qq,1]-hedge.freq)
        m2 <- list.fun[[nn2]](dT * quad.2d$nodes[qq,2])
        m.12 <- condMultiply(m1,m2,N.factors=N.factors)
        ab <- momentCondition(params.P,params.Q,jumpTransform = jumpTransform,condition.struct=m.12,conditional=TRUE,rtol=rtol,N.factors=N.factors,...)
        
        m.list$a <- c(m.list$a, ab$a)
        m.list$b <- rbind(m.list$b, ab$b)
        m.list$linCoeff <- c(m.list$linCoeff, ab$linCoeff * quad.2d$weights[qq] * (dT / hedge.freq)^2)
      }
      for (qq in 1:length(quad.1d$weights)) {
        m1 <- list.fun[[nn1]](dT * quad.1d$nodes[qq])
        m2 <- list.fun[[nn2]](dT * quad.1d$nodes[qq])
        m.12 <- condMultiply(m1,m2,N.factors=N.factors)
        
        ab <- momentCondition(params.P,params.Q,jumpTransform = jumpTransform,condition.struct=m.12,conditional=TRUE,rtol=rtol,N.factors=N.factors,...)
        # we could in principle substract means here to, but it's a smaller order, so not relevant
        m.list$a <- c(m.list$a,ab$a)
        m.list$b <- rbind(m.list$b,ab$b)
        m.list$linCoeff <- c(m.list$linCoeff,ab$linCoeff * quad.1d$weights[qq] * dT / hedge.freq)
      }
      cov.array[[nn1,nn2]] <- m.list
    }
  }
  
  return(list(mean.vec = mean.vec, cov.array = cov.array, cov.list = cov.array[lower.tri(diag(N.factors+1),diag=TRUE)]))
}