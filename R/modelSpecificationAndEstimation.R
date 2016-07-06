#' @title Affine model likelihood based on divergence prices
#' @name modLik
#' @description Functions for preparing model specification, evaluating the likelihood for the DSQ filtering.
#' @param data.structure \code{list} with fields \code{spec.mat}, \code{obs.data}, \code{noise.cov.cube}
#' @param model.spec \code{list} with fields \code{params.P}, \code{params.Q}, \code{jump.type}, \code{dt}, \code{N.factors}, \code{error}, \code{mkt}
#' @param for.estimation \code{logical}; determines return type (log-lik) or filtering result
#' @param filterFoo \code{function} that handles the filtering, must correspond to model specification and to provided observables.
#' @param N.points \code{integer}, number of integration points for double quadrature in moments of the state and stock price.
#' @param penalized \code{FALSE} by default, if \code{TRUE}, the Feller constraint is imposed as a penalty on the likelihood. Otherwise, infinity is returned (discontinuous at boundary).
#' @return \code{model_Likelihood} if \code{for.estimation==TRUE}: log-likelihood value (NOT negative of...), else: list with filtering results
#' @details Not much for now
#' @export

model_Likelihood <- function(data.structure, model.spec, for.estimation = FALSE, filterFoo = DSQ_sqrtFilter, N.points = 5, penalized = FALSE, penalty = 1e12){
  
  # Extract some variables from model specification
  N.factors <- model.spec$N.factors
  params.P <- model.spec$params.P
  params.Q <- model.spec$params.Q
  jump.type <- model.spec$jump.type
  time.dt <- data.structure$dt
  mkt.spec <- data.structure$mkt
  spec.mat <- data.structure$spec.mat
  data.obs <- data.structure$obs.data
  
  # Test parameters for vol factor Feller conditions
  feller.check <- model_fellerConditionCheck(params.P = params.P, params.Q = params.Q, N.factors = N.factors)
  logl.penalty <- 0
  if(!penalized){
    if(!any(feller.check$p)){
      print("P-measure Feller condition not met")
      return(-1e15)
    }
    if(!any(feller.check$q)){
      print("Q-measure Feller condition not met")
      return(-1e15)
    }
  } else {
    logl.penalty <- -penalty * sum( as.numeric(!feller.check$p) * pmin(0,feller.check$pval)^2 + as.numeric(!feller.check$q) * pmin(0, feller.check$qval)^2)
  }
  
  # solve ODEs for pricing
  ode.solutions <- tryCatch(
        Re(odeExtSolveWrap(
          u = cbind(unique(spec.mat[,"p"]),matrix(0,length(unique(spec.mat[,"p"])),N.factors)), 
          params.Q = params.Q, 
          mkt = mkt.spec, 
          jumpTransform = getPointerToJumpTransform(jump.type), 
          N.factors = N.factors, 
          mod.type = 'standard', 
          rtol = 1e-12, atol = 1e-28
        ))
      ,error = function(e){
        print(e)
        return(-1e15)
        })
  
  # Calculate model dynamics coefficients
  model.dynamics <- tryCatch(
      modelDynamics(
        params.P = params.P, 
        params.Q = params.Q, 
        dT = time.dt, 
        N.factors = N.factors, 
        jumpTransform = getPointerToJumpTransform(jump.type)$TF, 
        N.points = N.points,
        mod.type = 'standard', 
        rtol = 1e-12, 
        atol = 1e-30
      )
    , error = function(e){
      print(e)
      return(-1e15)
    })
  
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
        divNoiseCube = data.structure$noise.cov.cube,
        tVec = mkt.spec$t, 
        pVec = unique(spec.mat[,"p"])
      )
  
  # Long-term means under Q
  q.init.state <- matrix(0,nrow = 2, ncol = N.factors)
  for(kk in 1:N.factors){
    deriv.u <- matrix(0,nrow=2,ncol = N.factors+1)
    deriv.u[2,kk+1] <- 1e-5
    init.state.loc <- tryCatch(affineCF(
      u = deriv.u, 
      params.Q = params.Q, 
      params.P = NULL, 
      t.vec = 2, 
      v.0 = matrix(1,nrow=1,ncol=N.factors), 
      jumpTransform = getPointerToJumpTransform(jump.type)$TF, 
      N.factors = N.factors, 
      CGF = FALSE, 
      mod.type = 'standard', 
      atol = 1e-12, 
      rtol = 1e-10
    ),
    error = function(e){
      print(e)
      return(-1e15)
    })
    q.init.state[,kk] <- drop(init.state.loc)
  }
  
  q.init.state <- drop(q.init.state)
  if(is.null(dim(q.init.state))){
    q.init.state <- matrix(q.init.state,ncol=1)
  }
  q.init.state <- 1e5 * (q.init.state[2,]-q.init.state[1,])
  q.init.state <- Re(matrix(rep(q.init.state,2), ncol = 1))
  q.init.state[which(is.infinite(q.init.state))] <- 100.0
  
  # Initial state values -- long-term means under P
  init.state <- matrix(0,nrow = 2, ncol = N.factors)
  for(kk in 1:N.factors){
    deriv.u <- matrix(0,nrow=2,ncol = N.factors+1)
    deriv.u[2,kk+1] <- 1e-5
    init.state.loc <- tryCatch(affineCF(
      u = deriv.u, 
      params.Q = params.Q, 
      params.P = params.P, 
      t.vec = 2, 
      v.0 = matrix(1,nrow=1,ncol=N.factors), 
      jumpTransform = getPointerToJumpTransform(jump.type)$TF, 
      N.factors = N.factors, 
      CGF = FALSE, 
      mod.type = 'standard', 
      atol = 1e-12, 
      rtol = 1e-10
    ),
    error = function(e){
      print(e)
      return(-1e15)
    })
    init.state[,kk] <- drop(init.state.loc)
  }
  
  init.state <- drop(init.state)
  if(is.null(dim(init.state))){
    init.state <- matrix(init.state,ncol=1)
  }
  init.state <- 1e5 * (init.state[2,]-init.state[1,])
  init.state <- Re(matrix(rep(init.state,2), ncol = 1))
  
  # check stationarity under P and Q (approximate)
  eta.P <- unlist(sapply(params.P,function(x) x$eta))
  eta.Q <- unlist(sapply(params.Q,function(x) x$eta))
  stat.penalty <- 0
  stat.penalty <- as.numeric(init.state[1:N.factors]/eta.P >= 10) * pmax(0, init.state[1:N.factors]/eta.P - 10)^2
  stat.penalty <- c(stat.penalty, as.numeric(q.init.state[1:N.factors]/eta.Q >= 10) * pmax(0, q.init.state[1:N.factors]/eta.Q - 10)^2)
  stat.penalty <- -penalty * sum(stat.penalty)
  logl.penalty <- logl.penalty + stat.penalty

  init.vol <- matrix(0,2*N.factors,2*N.factors)
  init.vol[1:N.factors,1:N.factors] <- covMatFun(
          covListS = vol.cov.array[lower.tri(x = diag(N.factors),diag = T)], 
          covListDim = c(N.factors,N.factors),  
          currVol = init.state[1:N.factors,1,drop=F]
        )
  
  filtering.result <- tryCatch(filterFoo(
        dataMat = data.obs, 
        initState = init.state, 
        initProcCov = init.vol, 
        modelParams = list(transition = transition.parameters, observation = observation.parameters)
      ),
        error = function(e) {
          print(e)
          list(logL = -1e15)
        }
      )
  
  if(for.estimation){
    return(sum(filtering.result$logL + logl.penalty))
  } else {
    return(filtering.result)
  }
}

#' @describeIn modLik
#' @param par.vec vector with model parameter values
#' @param par.names parameter names, character vector equal in length to par
#' @param par.restr parameter equality restrictions, data.frame; par.vec and par.restr have to exhaust the model parameter set together.
#' @param N.factors integer, number of SV factors
#' @return \code{model_translateParameters} \code{list} with fields \code{P} and \code{Q}, input for all ODE calling functions.
#' @export
model_translateParameters <- function(par.vec, par.names = names(par.vec), par.restr, N.factors){
  
  # Assign parameter names
  names(par.vec) <- par.names
  
  # Check whether parameters are real-valued.
  stopifnot(is.double(par.vec))
  
  # Check for missing parameter names
  par.names <- names(par.vec)
  if(any(is.na(par.names))){
    stop("Unnamed parameters in par.vec vector")
  }
  
  # Initialise parameter holders; fill the easy ones.
  params.Q <- params.P <- vector(mode = "list", length = N.factors + 1)
  names(params.Q) <- names(params.P) <- c(as.character(1:N.factors),"jmp")
  
  # Fill the params vector with P parameters that are missing, they should be in the par.restr
  for(kk in 1:N.factors){
    which.constr <- which(grepl(pattern = paste0("P.",kk), x = par.restr$par.name))
    if(length(which.constr) > 0){
      loc.vec <- as.vector(par.restr[which.constr,"par.value"])
      names(loc.vec) <- par.restr[which.constr, "par.name"]
      # overwrite parameters which are in the par vector but are constrained
      for(nn in 1:length(loc.vec)){
        if(names(loc.vec)[nn] %in% names(par.vec)){
          par.vec[names(loc.vec)[nn]] <- loc.vec[nn]
        }else{
          par.vec <- c(par.vec, loc.vec[nn]) 
        } 
      }
    }
  }
  
  # Start filling params.P list with parameters
  for(kk in 1:N.factors){
    params.P[[as.character(kk)]] <- as.list(par.vec[which(grepl(pattern=paste0("P\\$",kk,"\\$"),names(par.vec)))])
    names(params.P[[as.character(kk)]]) <- substr(names(par.vec)[which(grepl(pattern=paste0("P\\$",kk,"\\$"),names(par.vec)))],5,8)
  }
  
  # Fill par.vec with jump P parameters that might be constrained
  which.jmp.constr <- which(grepl(pattern = "P.jmp", x = par.restr$par.name))
  if(length(which.jmp.constr) > 0){
    loc.vec <- as.vector(par.restr[which.jmp.constr, "par.value"])
    names(loc.vec) <- par.restr[which.jmp.constr,"par.name"]
    for(nn in 1:length(loc.vec)){
      if(names(loc.vec)[nn] %in% names(par.vec)){
        par.vec[names(loc.vec)[nn]] <- loc.vec[nn]
      }else{
        par.vec <- c(par.vec, loc.vec[nn]) 
      } 
    }
  }
  
  # fill params.P list with jmp parameters
  params.P$jmp <- as.list(par.vec[which(grepl(pattern = "P\\$jmp\\$", names(par.vec)) & !grepl(pattern = "lprop", names(par.vec)))])
  names(params.P$jmp) <- substring(names(par.vec)[which(grepl(pattern = "P\\$jmp\\$", names(par.vec)) & !grepl(pattern = "lprop", names(par.vec)))],7)
  params.P$jmp$lprop <- par.vec[which(grepl(pattern = "P\\$jmp\\$lprop", names(par.vec)))]
  params.P$jmp$lprop <- params.P$jmp$lprop[order(names(params.P$jmp$lprop))]
  
  # Fill par.vec with Q parameters that might be missing, they should be in par.restr or constrained to be the same as P
  for(kk in 1:N.factors){
    which.constr <- which(grepl(pattern = paste0("Q.",kk), x = par.restr$par.name))
    if(length(which.constr) > 0){
      loc.vec <- as.vector(par.restr[which.constr,"par.value"])
      for(nn in 1:length(loc.vec)){
        if(loc.vec[nn] == -999){
          loc.vec[nn] <- par.vec[sub("Q\\$","P\\$",par.restr[which.constr[nn],"par.name"])]
        } 
      }
      names(loc.vec) <- par.restr[which.constr, "par.name"]
      # overwrite parameters which are in the par vector but are constrained
      for(nn in 1:length(loc.vec)){
        if(names(loc.vec)[nn] %in% names(par.vec)){
          par.vec[names(loc.vec)[nn]] <- loc.vec[nn]
        }else{
          par.vec <- c(par.vec, loc.vec[nn])   
        }
      }
    }
  }
  
  # Fill params.Q list with sv parameters
  for(kk in 1:N.factors){
    params.Q[[as.character(kk)]] <- as.list(par.vec[which(grepl(pattern=paste0("Q\\$",kk,"\\$"),names(par.vec)))])
    names(params.Q[[as.character(kk)]]) <- substr(names(par.vec)[which(grepl(pattern=paste0("Q\\$",kk,"\\$"),names(par.vec)))],5,8)
  }
  
  # Fill Q jump parameters into par.vec, from constraints, including P and Q constraints
  which.jmp.constr <- which(grepl(pattern = "Q.jmp", x = par.restr$par.name))
  if(length(which.jmp.constr) > 0){
    loc.vec <- as.vector(par.restr[which.jmp.constr, "par.value"])
    names(loc.vec) <- par.restr[which.jmp.constr,"par.name"]
    for(pp in 1:length(loc.vec)){
      if(loc.vec[pp] == -999){
        loc.name <- names(loc.vec)[pp]
        loc.name <- sub("Q","P",loc.name)
        loc.vec[pp] <- par.vec[loc.name]
      }
    }
    for(nn in 1:length(loc.vec)){
      if(names(loc.vec)[nn] %in% names(par.vec)){
        par.vec[names(loc.vec)[nn]] <- loc.vec[nn]
      }else{
        par.vec <- c(par.vec,loc.vec[nn]) 
      }
    }
  }
  
  # Fill params.Q$jmp with parameters
  params.Q$jmp <- as.list(par.vec[which(grepl(pattern = "Q\\$jmp\\$", names(par.vec)) & !grepl(pattern = "lprop", names(par.vec)))])
  names(params.Q$jmp) <- substring(names(par.vec)[which(grepl(pattern = "Q\\$jmp\\$", names(par.vec)) & !grepl(pattern = "lprop", names(par.vec)))],7)
  params.Q$jmp$lprop <- par.vec[which(grepl(pattern = "Q\\$jmp\\$lprop", names(par.vec)))]
  params.Q$jmp$lprop <- params.Q$jmp$lprop[order(names(params.Q$jmp$lprop))]
  
  # Check if all P factor parameters are defined. They all have to be numeric.
  svNames <- c("kpp","eta","lmb","rho","phi")
  jmpNames <- c("lvec",paste0("lprop[",1:N.factors,"]"),"rhoc","muYc","sigmaYc","muSc")
  
  for(M in c("P","Q")){
    for(pp in svNames){
      for(fs in as.character(1:N.factors)){
        if(!is.numeric(eval(parse(text = paste0("params.",M,"[[as.character(",fs,")]]$",pp))))){
          stop(paste0("Error at parameter translation: parameter ", paste0(M,"$",fs,"$",pp) ," has non-numeric value"))
        } 
      }
    }
  }
  
  for(M in c("P","Q")){
    for(pp in jmpNames){   
      if(!is.numeric(eval(parse(text = paste0("params.",M,"$","jmp","$",pp))))){
        stop(paste0("Error at parameter translation: parameter ", paste0(M,"$","jmp","$",pp) ," has non-numeric value"))
      } 
    }
  }
  
  par.list <- list(P = params.P, Q = params.Q)
  return(par.list)
}

#' @describeIn modLik
#' @return \code{model_makeDefaultParameterStructures} returns \code{data.frame} \code{par.restr} and \code{character} vector \code{par.names}
#' @export
model_makeDefaultParameterStructures <- function(N.factors, pq.equality = c("Q$jmp$lvec",paste0("Q$jmp$lprop.",1:N.factors),paste0("Q$",1:N.factors,"$eta"))){
  
  # Fill in the basic constraints we agreed upon
  loc.df <- data.frame(par.name = character(), par.value = numeric())
  for(kk in 1:N.factors){
    loc.df <- rbind(loc.df ,data.frame(par.name = paste0("P$",kk,"$eta"), par.value = 1)) # normalisation
  }
  for(nm in pq.equality){
    loc.df <- rbind(loc.df, data.frame(par.name = nm, par.value = -999))
  }
  
  # Add required P-Q constraints
  for(kk in 1:N.factors){
    loc.df <- rbind(loc.df, data.frame(par.name = paste0("Q$",kk,"$rho"), par.value = -999))
    loc.df <- rbind(loc.df, data.frame(par.name = paste0("Q$",kk,"$lmb"), par.value = -999))
    loc.df <- rbind(loc.df, data.frame(par.name = paste0("Q$",kk,"$phi"), par.value = -999))
    loc.df <- rbind(loc.df, data.frame(par.name = paste0("P$",kk,"$erp"), par.value = 0))
  }
  par.restr <- loc.df
  
  #   names <- c('P$svFast$kpp','P$svFast$lmb','P$svFast$rho','P$svFast$eta','P$svSlow$kpp','P$svSlow$lmb','P$svSlow$rho','P$svSlow$eta','Q$svFast$kpp','Q$svFast$eta','Q$svSlow$kpp','Q$svSlow$eta','P$jmp$lvec','P$jmp$lprop.1','P$jmp$lprop.2','P$jmp$muY','P$jmp$sigmaY','Q$jmp$muY','Q$jmp$sigmaY')
    
  var.names <- as.vector(vapply(1:N.factors,FUN.VALUE=character(5),function(x) return(paste0(paste0(x,"$"),c("kpp","eta","rho","phi","lmb")))))
  var.names <- c("1$erp0", var.names)
    
  P.names <- paste0("P$",var.names[which(!grepl(pattern = "eta", var.names))])
  Q.names <- paste0("Q$", var.names[which(grepl(pattern = "kpp|eta", var.names))])
    
  int.names.P <- paste0("P$jmp$",c("lvec",paste0("lprop.",1:N.factors)))
  int.names.Q <- paste0("Q$jmp$",c("lvec",paste0("lprop.",1:N.factors)))
    
  jmp.names <- c(paste0("P$jmp$", c("muSc","muYc","sigmaYc","rhoc")), paste0("Q$jmp$", c("muSc","muYc","sigmaYc","rhoc")))
    
  names <- c(P.names, int.names.P, jmp.names, Q.names, int.names.Q)
  
  which.restr <- which(names %in% par.restr$par.name)
  if(length(which.restr) > 0){
    names <- names[-which.restr]
  }
  return(list(par.names=names, par.restr = par.restr))
}

#' @describeIn modLik
#' @return \code{model_wrapLikelihood} wraps the likelihood function so that it only accepts a parameter vector argument -- use this for optimizers that do not allow passing extra arguments to the optimised function.
#' @export

model_wrapLikelihood <- function(data.structure, model.spec, for.estimation = FALSE, filterFoo = DSQ_sqrtFilter, N.points = 5, penalized = FALSE, penalty){
  
  retFoo <- function(par.vec){
    par.list <- model_translateParameters(par.vec = par.vec, par.names = model.spec$par.names, par.restr = model.spec$par.restr, N.factors = model.spec$N.factors)
    model.spec$params.P <- par.list$P
    model.spec$params.Q <- par.list$Q
    
    logLik <- model_Likelihood(data.structure = data.structure, model.spec = model.spec, for.estimation = for.estimation, filterFoo = filterFoo, N.points = N.points, penalized = penalized, penalty = penalty)
    
    return(logLik)
  }
 
  return(retFoo)
}

#' @describeIn modLik
#' @return \code{model_fellerConditionCheck} list with two logical vectors reporting whether the Feller conditions are satisfied
#' @export
model_fellerConditionCheck <- function(params.P, params.Q, N.factors){
    
  p <- logical(N.factors)
  pval <- numeric(N.factors)
  for(ff in 1:N.factors){
    pval[ff] <- 2 * params.P[[as.character(ff)]]$kpp * params.P[[as.character(ff)]]$eta - params.P[[as.character(ff)]]$lmb^2
    p[ff] <- pval[ff] > 0
  }
  
  q <- logical(N.factors)
  qval <- numeric(N.factors)
  for(ff in 1:N.factors){
    qval[ff] <- 2 * params.Q[[as.character(ff)]]$kpp * params.Q[[as.character(ff)]]$eta - params.Q[[as.character(ff)]]$lmb^2
    q[ff] <- qval[ff] > 0
  }
  return(list(p=p,q=q,pval= pval, qval= qval))
}

#' @describeIn modLik
#' @param noisePar vector of noise variance magnitudes, equal to number of observed pfolios
#' @return \code{model_Likelihood_extraNoise} \code{list} with fields \code{P} and \code{Q}, input for all ODE calling functions.
#' @export

model_Likelihood_extraNoise <- function(data.structure, model.spec, for.estimation = FALSE, filterFoo = DSQ_sqrtFilter, N.points = 5, penalized = FALSE, penalty = 1e12){
  
  # Extract some variables from model specification
  N.factors <- model.spec$N.factors
  params.P <- model.spec$params.P
  params.Q <- model.spec$params.Q
  jump.type <- model.spec$jump.type
  time.dt <- data.structure$dt
  mkt.spec <- data.structure$mkt
  spec.mat <- data.structure$spec.mat
  data.obs <- data.structure$obs.data
  noise.par <- model.spec$noise.par
  
  # Test parameters for vol factor Feller conditions
  feller.check <- model_fellerConditionCheck(params.P = params.P, params.Q = params.Q, N.factors = N.factors)
  logl.penalty <- 0
  if(!penalized){
    if(!any(feller.check$p)){
      print("P-measure Feller condition not met")
      return(-1e15)
    }
    if(!any(feller.check$q)){
      print("Q-measure Feller condition not met")
      return(-1e15)
    }
  } else {
    logl.penalty <- -penalty * sum( as.numeric(!feller.check$p) * pmin(0,feller.check$pval)^2 + as.numeric(!feller.check$q) * pmin(0, feller.check$qval)^2)
  }
  
  # solve ODEs for pricing
  ode.solutions <- tryCatch(
    Re(odeExtSolveWrap(
      u = cbind(unique(spec.mat[,"p"]),matrix(0,length(unique(spec.mat[,"p"])),N.factors)), 
      params.Q = params.Q, 
      mkt = mkt.spec, 
      jumpTransform = getPointerToJumpTransform(jump.type), 
      N.factors = N.factors, 
      mod.type = 'standard', 
      rtol = 1e-12, atol = 1e-28
    ))
    ,error = function(e){
      print(e)
      return(-1e15)
    })
  
  # Calculate model dynamics coefficients
  model.dynamics <- tryCatch(
    modelDynamics(
      params.P = params.P, 
      params.Q = params.Q, 
      dT = time.dt, 
      N.factors = N.factors, 
      jumpTransform = getPointerToJumpTransform(jump.type)$TF, 
      N.points = N.points,
      mod.type = 'standard', 
      rtol = 1e-12, 
      atol = 1e-30
    )
    , error = function(e){
      print(e)
      return(-1e15)
    })
  
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
    divNoiseCube = data.structure$noise.cov.cube,
    tVec = mkt.spec$t, 
    pVec = unique(spec.mat[,"p"]),
    errSdParVec = noise.par
  )
  
  # Long-term means under Q
  q.init.state <- matrix(0,nrow = 2, ncol = N.factors)
  for(kk in 1:N.factors){
    deriv.u <- matrix(0,nrow=2,ncol = N.factors+1)
    deriv.u[2,kk+1] <- 1e-5
    init.state.loc <- tryCatch(affineCF(
      u = deriv.u, 
      params.Q = params.Q, 
      params.P = NULL, 
      t.vec = 2, 
      v.0 = matrix(1,nrow=1,ncol=N.factors), 
      jumpTransform = getPointerToJumpTransform(jump.type)$TF, 
      N.factors = N.factors, 
      CGF = FALSE, 
      mod.type = 'standard', 
      atol = 1e-12, 
      rtol = 1e-10
    ),
    error = function(e){
      print(e)
      return(-1e15)
    })
    q.init.state[,kk] <- drop(init.state.loc)
  }
  
  q.init.state <- drop(q.init.state)
  if(is.null(dim(q.init.state))){
    q.init.state <- matrix(q.init.state,ncol=1)
  }
  q.init.state <- 1e5 * (q.init.state[2,]-q.init.state[1,])
  q.init.state <- Re(matrix(rep(q.init.state,2), ncol = 1))
  q.init.state[which(is.infinite(q.init.state))] <- 100.0
  
  # Initial state values -- long-term means under P
  init.state <- matrix(0,nrow = 2, ncol = N.factors)
  for(kk in 1:N.factors){
    deriv.u <- matrix(0,nrow=2,ncol = N.factors+1)
    deriv.u[2,kk+1] <- 1e-5
    init.state.loc <- tryCatch(affineCF(
      u = deriv.u, 
      params.Q = params.Q, 
      params.P = params.P, 
      t.vec = 2, 
      v.0 = matrix(1,nrow=1,ncol=N.factors), 
      jumpTransform = getPointerToJumpTransform(jump.type)$TF, 
      N.factors = N.factors, 
      CGF = FALSE, 
      mod.type = 'standard', 
      atol = 1e-12, 
      rtol = 1e-10
    ),
    error = function(e){
      print(e)
      return(-1e15)
    })
    init.state[,kk] <- drop(init.state.loc)
  }
  
  init.state <- drop(init.state)
  if(is.null(dim(init.state))){
    init.state <- matrix(init.state,ncol=1)
  }
  init.state <- 1e5 * (init.state[2,]-init.state[1,])
  init.state <- Re(matrix(rep(init.state,2), ncol = 1))
  
  # check stationarity under P and Q (approximate)
  eta.P <- unlist(sapply(params.P,function(x) x$eta))
  eta.Q <- unlist(sapply(params.Q,function(x) x$eta))
  stat.penalty <- 0
  stat.penalty <- as.numeric(init.state[1:N.factors]/eta.P >= 10) * pmax(0, init.state[1:N.factors]/eta.P - 10)^2
  stat.penalty <- c(stat.penalty, as.numeric(q.init.state[1:N.factors]/eta.Q >= 10) * pmax(0, q.init.state[1:N.factors]/eta.Q - 10)^2)
  stat.penalty <- -penalty * sum(stat.penalty)
  logl.penalty <- logl.penalty + stat.penalty
  #   if(any(init.state[1:N.factors]/eta.P >= 5 )){
  #     print("Non-stationary P-measure volatility")
  #     return(-1e15)
  #   }
  
  init.vol <- matrix(0,2*N.factors,2*N.factors)
  init.vol[1:N.factors,1:N.factors] <- covMatFun(
    covListS = vol.cov.array[lower.tri(x = diag(N.factors),diag = T)], 
    covListDim = c(N.factors,N.factors),  
    currVol = init.state[1:N.factors,1,drop=F]
  )
  
  filtering.result <- tryCatch(filterFoo(
    dataMat = data.obs, 
    initState = init.state, 
    initProcCov = init.vol, 
    modelParams = list(transition = transition.parameters, observation = observation.parameters)
  ),
  error = function(e) {
    print(e)
    list(logL = -1e15)
  }
  )
  
  if(for.estimation){
    return(sum(filtering.result$logL + logl.penalty))
  } else {
    return(filtering.result)
  }
}

#' @describeIn modLik
#' @return \code{model_wrapLikelihood_extraNoise} wraps the likelihood function with extra noise so that it only accepts a parameter vector argument -- use this for optimizers that do not allow passing extra arguments to the optimised function.
#' @export

model_wrapLikelihood_extraNoise <- function(data.structure, model.spec, for.estimation = FALSE, filterFoo = DSQ_sqrtFilter, N.points = 5, penalized = FALSE, penalty){
  
  retFoo <- function(par.vec){
    noise.par <- tail(par.vec,nrow(data.structure$spec.mat))
    par.vec <- head(par.vec,-nrow(data.structure$spec.mat))
    model.spec$noise.par <- noise.par
    par.list <- model_translateParameters(par.vec = par.vec, par.names = model.spec$par.names, par.restr = model.spec$par.restr, N.factors = model.spec$N.factors)
    model.spec$params.P <- par.list$P
    model.spec$params.Q <- par.list$Q
    
    logLik <- model_Likelihood_extraNoise(data.structure = data.structure, model.spec = model.spec, for.estimation = for.estimation, filterFoo = filterFoo, N.points = N.points, penalized = penalized, penalty = penalty)
    
    return(logLik)
  }
  
  return(retFoo)
}

#' @describeIn modLik
#' @param noisePar vector of noise variance magnitudes, equal to number of observed pfolios
#' @return \code{model_Likelihood_portfolio_extraNoise} \code{list} with fields \code{P} and \code{Q}, input for all ODE calling functions.
#' @export

model_Likelihood_portfolio_extraNoise <- function(data.structure, model.spec, for.estimation = FALSE, filterFoo = divergenceModelR:::portfolio_sqrtFilter, N.points = 5, penalized = FALSE, penalty = 1e12, N.GL.points = 64){
  
  # quadrature setup
  gl <- statmod::gauss.quad(n = N.GL.points, kind = "laguerre", alpha = 0)
  
  # Extract some variables from model specification
  N.factors <- model.spec$N.factors
  params.P <- model.spec$params.P
  params.Q <- model.spec$params.Q
  jump.type <- model.spec$jump.type
  time.dt <- data.structure$dt
  mkt.list <- data.structure$mkt.list
  wts.list <- data.structure$wts.list
  strikeMat.list <- data.structure$strikeMat.list
  data.obs <- data.structure$obs.data
  noise.par <- model.spec$noise.par
  
  # Test parameters for vol factor Feller conditions
  feller.check <- model_fellerConditionCheck(params.P = params.P, params.Q = params.Q, N.factors = N.factors)
  logl.penalty <- 0
  if(!penalized){
    if(!any(feller.check$p)){
      print("P-measure Feller condition not met")
      return(-1e15)
    }
    if(!any(feller.check$q)){
      print("Q-measure Feller condition not met")
      return(-1e15)
    }
  } else {
    logl.penalty <- -penalty * sum( as.numeric(!feller.check$p) * pmin(0,feller.check$pval)^2 + as.numeric(!feller.check$q) * pmin(0, feller.check$qval)^2)
  }
  
  # collate maturities into mkt.spec
  mkt.spec <- do.call(what = rbind, args = mkt.list)
  mkt.spec <- cbind(t = unique(mkt.spec[,"t"]), r = 0, q = 0, p = 0)
  mkt.spec <- mkt.spec[order(mkt.spec[,"t"]),]
  
  # solve ODEs for pricing
  ode.solutions <- tryCatch(
    jumpDiffusionODEs(
      u = cbind(1i * gl$nodes,matrix(0,length(gl$nodes),N.factors)), 
      params = params.Q, 
      mkt = as.data.frame(mkt.spec), 
      jumpTransform = getPointerToJumpTransform(jump.type)$TF, 
      N.factors = N.factors, 
      mod.type = 'standard', 
      rtol = 1e-12, atol = 1e-28
    )
    ,error = function(e){
      print(e)
      return(-1e15)
    })
  
  ode.solutions.list <- vector(mode = "list", length = length(mkt.list))
  for(ind in 1:length(mkt.list)){
    loc.t <- mkt.list[[ind]][,"t"]
    t.ind <- which(mkt.spec[,"t"] %in% loc.t)
    ode.solutions.list[[ind]] <- ode.solutions[,t.ind,,drop=FALSE]
  }
  
  # Calculate model dynamics coefficients
  model.dynamics <- tryCatch(
    modelDynamics(
      params.P = params.P, 
      params.Q = params.Q, 
      dT = time.dt, 
      N.factors = N.factors, 
      jumpTransform = getPointerToJumpTransform(jump.type)$TF, 
      N.points = N.points,
      mod.type = 'standard', 
      rtol = 1e-12, 
      atol = 1e-30
    )
    , error = function(e){
      print(e)
      return(-1e15)
    })
  
  # Transition equation setup
  vol.cov.array <- model.dynamics$cov.array[-1,-1]
  transition.parameters <- list(
    mean.vec = model.dynamics$mean.vec[-1], 
    cov.array = vol.cov.array[lower.tri(matrix(0,N.factors,N.factors),diag = T)]
  )
  
  # Observation equation setup
  observation.parameters <- list(
    stockParams = list(mean.vec = model.dynamics$mean.vec, cov.array = model.dynamics$cov.array[lower.tri(diag(1+N.factors),diag = T)]), 
    cfCoeffs = ode.solutions.list,
    divNoiseCube = data.structure$noise.cov.cube,
    strikeMats = strikeMat.list,
    mkts = mkt.list,
    wts = wts.list,
    errSdParVec = noise.par,
    quadWeights = gl$weights,
    quadNodes = gl$nodes
  )
  
  # Long-term means under Q
  q.init.state <- matrix(0,nrow = 2, ncol = N.factors)
  for(kk in 1:N.factors){
    deriv.u <- matrix(0,nrow=2,ncol = N.factors+1)
    deriv.u[2,kk+1] <- 1e-5
    init.state.loc <- tryCatch(affineCF(
      u = deriv.u, 
      params.Q = params.Q, 
      params.P = NULL, 
      t.vec = 2, 
      v.0 = matrix(1,nrow=1,ncol=N.factors), 
      jumpTransform = getPointerToJumpTransform(jump.type)$TF, 
      N.factors = N.factors, 
      CGF = FALSE, 
      mod.type = 'standard', 
      atol = 1e-12, 
      rtol = 1e-10
    ),
    error = function(e){
      print(e)
      return(-1e15)
    })
    q.init.state[,kk] <- drop(init.state.loc)
  }
  
  q.init.state <- drop(q.init.state)
  if(is.null(dim(q.init.state))){
    q.init.state <- matrix(q.init.state,ncol=1)
  }
  q.init.state <- 1e5 * (q.init.state[2,]-q.init.state[1,])
  q.init.state <- Re(matrix(rep(q.init.state,2), ncol = 1))
  q.init.state[which(is.infinite(q.init.state))] <- 100.0
  
  # Initial state values -- long-term means under P
  init.state <- matrix(0,nrow = 2, ncol = N.factors)
  for(kk in 1:N.factors){
    deriv.u <- matrix(0,nrow=2,ncol = N.factors+1)
    deriv.u[2,kk+1] <- 1e-5
    init.state.loc <- tryCatch(affineCF(
      u = deriv.u, 
      params.Q = params.Q, 
      params.P = params.P, 
      t.vec = 2, 
      v.0 = matrix(1,nrow=1,ncol=N.factors), 
      jumpTransform = getPointerToJumpTransform(jump.type)$TF, 
      N.factors = N.factors, 
      CGF = FALSE, 
      mod.type = 'standard', 
      atol = 1e-12, 
      rtol = 1e-10
    ),
    error = function(e){
      print(e)
      return(-1e15)
    })
    init.state[,kk] <- drop(init.state.loc)
  }
  
  init.state <- drop(init.state)
  if(is.null(dim(init.state))){
    init.state <- matrix(init.state,ncol=1)
  }
  init.state <- 1e5 * (init.state[2,]-init.state[1,])
  init.state <- Re(matrix(rep(init.state,2), ncol = 1))
  
  # check stationarity under P and Q (approximate)
  eta.P <- unlist(sapply(params.P,function(x) x$eta))
  eta.Q <- unlist(sapply(params.Q,function(x) x$eta))
  stat.penalty <- 0
  stat.penalty <- as.numeric(init.state[1:N.factors]/eta.P >= 10) * pmax(0, init.state[1:N.factors]/eta.P - 10)^2
  stat.penalty <- c(stat.penalty, as.numeric(q.init.state[1:N.factors]/eta.Q >= 10) * pmax(0, q.init.state[1:N.factors]/eta.Q - 10)^2)
  stat.penalty <- -penalty * sum(stat.penalty)
  logl.penalty <- logl.penalty + stat.penalty
  
  init.vol <- matrix(0,2*N.factors,2*N.factors)
  init.vol[1:N.factors,1:N.factors] <- covMatFun(
    covListS = vol.cov.array[lower.tri(x = diag(N.factors),diag = T)], 
    covListDim = c(N.factors,N.factors),  
    currVol = init.state[1:N.factors,1,drop=F]
  )
  
  filtering.result <- tryCatch(filterFoo(
    dataMat = data.obs, 
    initState = init.state, 
    initProcCov = init.vol, 
    modelParams = list(transition = transition.parameters, observation = observation.parameters)
  ),
  error = function(e) {
    print(e)
    list(logL = -1e15)
  }
  )
  
  if(for.estimation){
    return(sum(filtering.result$logL + logl.penalty))
  } else {
    return(filtering.result)
  }
}

#' @describeIn modLik
#' @return \code{model_wrapLikelihood_portfolio_extraNoise} wraps the likelihood function with extra noise so that it only accepts a parameter vector argument -- use this for optimizers that do not allow passing extra arguments to the optimised function.
#' @export

model_wrapLikelihood_portfolio_extraNoise <- function(data.structure, model.spec, for.estimation = FALSE, filterFoo = DSQ_sqrtFilter, N.points = 5, penalized = FALSE, penalty, N.GL.points = 64){
  
  retFoo <- function(par.vec){
    noise.par <- tail(par.vec,ncol(data.structure$obs.data)-1)
    par.vec <- head(par.vec,-(ncol(data.structure$obs.data)-1))
    model.spec$noise.par <- noise.par
    par.list <- model_translateParameters(par.vec = par.vec, par.names = model.spec$par.names, par.restr = model.spec$par.restr, N.factors = model.spec$N.factors)
    model.spec$params.P <- par.list$P
    model.spec$params.Q <- par.list$Q
    
    logLik <- model_Likelihood_portfolio_extraNoise(data.structure = data.structure, model.spec = model.spec, for.estimation = for.estimation, filterFoo = filterFoo, N.points = N.points, penalized = penalized, penalty = penalty, N.GL.points = N.GL.points)
    
    return(logLik)
  }
  
  return(retFoo)
}