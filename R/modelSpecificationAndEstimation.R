#'
#' @title Affine model likelihood
#' @name modLik
#' @description Based on a specification object, this function evaluates the log-likelihood of the filtering result
#' @param data.structure \code{list} with fields \code{spec.mat}, \code{obs.data}
#' @param model.spec \code{list} with fields \code{params.P}, \code{params.Q}, \code{jump.type}, \code{dt}, \code{N.factors}, \code{error}, \code{mkt}
#' @param for.estimation \code{logical}; determines return type (log-lik) or filtering result
#' @param filterFoo \code{function} that handles the filtering, must correspond to model specification and to provided observables.
#' @param par.df \code{data.frame} with columns \code{par.name} with format \code{P$1$kpp} and \code{par.value} with values.
#' @return If \code{for.estimation==TRUE}: log-likelihood value (NOT negative of...), else: list with filtering results
#' @details Not much for now
#' @export

modelLikelihood <- function(data.structure, model.spec, for.estimation = FALSE, filterFoo = DSQ_sqrtFilter){
  
  # Extract some variables from model specification
  N.factors <- model.spec$N.factors
  params.P <- model.spec$params.P
  params.Q <- model.spec$params.Q
  jump.type <- model.spec$jump.type
  mkt.spec <- model.spec$mkt
  error.spec <- model.spec$error
  time.dt <- model.spec$dt
  spec.mat <- data.structure$spec.mat
  data.obs <- data.structure$obs.data
  
  # solve ODEs for pricing
  ode.solutions <- Re(odeExtSolveWrap(
        u = cbind(unique(spec.mat[,"p"]),matrix(0,length(unique(spec.mat[,"p"])),N.factors)), 
        params.P = parListHestonLev$P, 
        params.Q = parListHestonLev$Q, 
        mkt = mkt.spec, 
        jumpTransform = getPointerToJumpTransform(jump.type), 
        N.factors = 1, 
        mod.type = 'standard', 
        rtol = 1e-8, atol = 1e-12
      ))
  
  # Calculate model dynamics coefficients
  model.dynamics <- modelDynamics(
        params.P = params.P, 
        params.Q = params.Q, 
        dT = time.dt, 
        N.factors = N.factors, 
        jumpTransform = getPointerToJumpTransform(jump.type)$TF, 
        N.points = 2,
        mod.type = 'standard', 
        rtol = 1e-8, 
        atol = 1e-12
      )
  
  # Transition equation setup
  vol.cov.array <- model.dynamics$cov.array[-1,-1]
  transition.parameters <- list(
        mean.vec = model.dynamics$mean.vec[-1], 
        cov.array = vol.cov.array[lower.tri(matrix(0,N.factors,N.factors),diag = T)]
      )
  
  # Observation equation setup
  observation.parameters <- list(
        stockParams = list(mean.vec = model.dynamics$mean.vec, cov.array = model.dynamics$cov.array[lower.tri(diag(1+N.factors),diag = T)]), 
        cfCoeffs = ode.solutions, 
        tVec = mkt.spec$t, 
        pVec = unique(spec.mat[,"p"]), 
        cVec = error.spec$cVec, 
        bVec = error.spec$bVec
      )
  
  # Initial state value
  init.state <- affineCF(
        u = rbind(rep(0,N.factors+1),matrix(c(0,rep(1e-5,N.factors)),nrow = 1, ncol = N.factors+1)), 
        params.Q = params.Q, 
        params.P = params.P, 
        t.vec = 10, 
        v.0 = matrix(1,nrow=1,ncol=N.factors), 
        jumpTransform = getPointerToJumpTransform(jump.type)$TF, 
        N.factors = N.factors, 
        CGF = FALSE, 
        mod.type = 'standard', 
        atol = 1e-12, 
        rtol = 1e-8
      )
  init.state <- drop(init.state)
  if(is.null(dim(init.state))){
    init.state <- matrix(init.state,ncol=1)
  }
  init.state <- 1e5 * (init.state[2,]-init.state[1,])
  init.state <- matrix(rep(init.state,2), ncol = 1)
  init.vol <- matrix(0,2*N.factors,2*N.factors)
  init.vol[1:N.factors,1:N.factors] <- covMatFun(
          covListS = vol.cov.array[lower.tri(x = diag(N.factors),diag = T)], 
          covListDim = c(N.factors,N.factors),  
          currVol = initState[1:N.factors,1,drop=F]
        )
  
  filtering.result <- filterFoo(
        dataMat = data.obs, 
        initState = init.state, 
        initProcCov = init.vol, 
        modelParams = list(transition = transition.parameters, observation = observation.parameters)
      )
  
  if(for.estimation){
    return(sum(filtering.result$logL))
  } else {
    return(filtering.result)
  }
}

#' @describeIn modLik
model_translateParameters <- function(par.df){
  
}