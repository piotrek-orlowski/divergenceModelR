#' @title Affine model evaluation
#' @name modelEvaluation
#' @description Functions for evaluating fit of models estimated on divergence data
#' @param data.structure \code{list} with fields \code{spec.mat}, \code{obs.data}, \code{noise.cov.cube}
#' @param model.spec \code{list} with fields \code{params.P}, \code{params.Q}, \code{jump.type}, \code{dt}, \code{N.factors}, \code{error}, \code{mkt}
#' @return \code{evaluation_fittedPrices} returns a list with fields \code{obs} (the observed returns and swap prices) and \code{fit} (the predicted returns and fitted swap prices.)
#' @export

evaluation_fittedPrices <- function(par.vec, model_Likelihood = model_Likelihood, par.list = NULL, data.structure, model.spec, filterFoo = DSQ_sqrtFilter, handlerFoo = affineObservationStateHandler, N.points = 5){
  if(is.null(par.vec) & is.null(par.list)){
    stop("No parameters supplied: both par.vec and par.list are NULL")
  }
  
  # Get parameters from vector into lists
  par.list <- model_translateParameters(par.vec = par.vec[1:length(model.spec$par.names)], par.names = model.spec$par.names, par.restr = model.spec$par.restr, N.factors = model.spec$N.factors)
  model.spec$params.P <- par.list$P
  model.spec$params.Q <- par.list$Q
  
  # Evaluate likelihood
  logLik <- model_Likelihood(data.structure = data.structure, model.spec = model.spec, for.estimation = FALSE, filterFoo = filterFoo, N.points = N.points, penalized = FALSE)
  
  # Prepare data for calling the affine observation state handler
  model.dynamics <- modelDynamics(params.P = model.spec$params.P, params.Q = model.spec$params.Q, dT = data.structure$dt, N.factors = model.spec$N.factors, jumpTransform = getPointerToJumpTransform(model.spec$jump.type)$TF, N.points = N.points, mod.type = 'standard')
  
  sol.derivs <- Re(odeExtSolveWrap(u = cbind(unique(data.structure$spec.mat$p),matrix(0,nrow = length(unique(data.structure$spec.mat$p)), ncol = model.spec$N.factors)), params.P = NULL, params.Q = model.spec$params.Q, mkt = data.structure$mkt, jumpTransform = getPointerToJumpTransform(model.spec$jump.type), N.factors = model.spec$N.factors, mod.type = 'standard', rtol = 1e-12, atol = 1e-28))
  
  obsParameters <- list(stockParams = list(mean.vec = model.dynamics$mean.vec, cov.array = model.dynamics$cov.array[lower.tri(diag(1+model.spec$N.factors) ,diag = T)]), cfCoeffs = sol.derivs, tVec = data.structure$mkt$t, pVec = unique(data.structure$spec.mat$p), divNoiseCube = data.structure$noise.cov.cube, errSdParVec = model.spec$noise.par)

  swap.fit <- t(handlerFoo(stateMat = t(logLik$estimState[-1,,drop=F]), modelParameters = obsParameters, iterCount = 1)$yhat)
  
  pred.ret <- apply(head(logLik$estimState[,1:model.spec$N.factors,drop=FALSE],-1), 1, meanVecFun, meanListS = model.dynamics$mean.vec[1])
  
  return(list(obs = data.structure$obs.data, fit = swap.fit[,-1], exp.ret = pred.ret, lik.object = logLik))
}

#' @rdname modelEvaluation
#' @details \code{evaluation_plotting} plots the following: (1) time series of filtered factor values, (2) time series of fitted and observed swap rates, (3) Q-Q plots of filtered factors and bootstrap sample from the estimated model factors, (4) autocorrelation functions of pricing errors, (5)
#' @return \code{evaluation_plotting} does not return.
#' @export
evaluation_plotting <- function(par.vec, model_Likelihood= model_Likelihood, par.list = NULL, data.structure, model.spec, filterFoo = DSQ_sqrtFilter, handlerFoo = affineObservationStateHandler, N.points = 5, plot.path){
  
  # par lists
  par.list <- model_translateParameters(par.vec = par.vec[1:length(model.spec$par.names)], par.names = model.spec$par.names, par.restr = model.spec$par.restr, model.spec$N.factors)
  model.spec$params.P <- par.list$P
  model.spec$params.Q <- par.list$Q
  
  # fitted values
  model.results <- evaluation_fittedPrices(par.vec = par.vec, model_Likelihood = model_Likelihood, par.list = par.list, data.structure = data.structure, model.spec = model.spec, filterFoo = filterFoo, N.points = N.points, handlerFoo = handlerFoo)
  
  fit <- model.results$fit
  
  dates <- data.structure$dates$date
  
  # plot fitted values
  N.plots <- ncol(fit)
  for(npl in 1:N.plots){
    tikzDevice::tikz(file = paste0(plot.path,"swap-price-",npl,".tex"), width = 5, height = 3.7, standAlone = TRUE, sanitize = FALSE)
      loc.data <- cbind(data.structure$obs.data[,npl+1], fit[,npl])
      title <- paste0("Price of ",as.character(data.structure$spec.mat$type[npl])," swap rate \n $\\tau=$",sprintf("%1.2f",data.structure$spec.mat$t[npl]),", $p=$", sprintf("%1.1f",data.structure$spec.mat$p[npl]))
      plot(dates, loc.data[,1], ylim = range(loc.data), type = 'l', lwd = 1.5, col = 'black', main = title, xlab = "Date", ylab = "Swap rate")
      lines(dates, loc.data[,2], lwd = 1.5, lty = 2, col = "darkorange")
    dev.off()
    tools::texi2pdf(file = paste0(plot.path,"swap-price-",npl,".tex"), clean = TRUE, quiet = TRUE)
  }
  
  # plot filtered states
  for(npl in 1:model.spec$N.factors){
    tikzDevice::tikz(file = paste0(plot.path,"vol-factor-",npl,".tex"), width = 5, height = 3.7, standAlone = TRUE, sanitize = FALSE)
      if(model.spec$params.P[[as.character(npl)]]$phi == 0){
        title <- paste0("Filtered factor $V_",npl,"$")
        plot(dates, model.results$lik.object$estimState[-1,npl], type = 'l', col = 'black', lwd = 1.5, xlab = 'Date', ylab =paste0("$V_",npl,"$"), main = title)  
      } else {
        title <- paste0("Filtered factor $\\phi^2_",npl,"V_",npl,"$")
        plot(dates, model.spec$params.P[[as.character(npl)]]$phi^2*model.results$lik.object$estimState[-1,npl], type = 'l', col = 'black', lwd = 1.5, xlab = 'Date', ylab =paste0("$\\phi^2_",npl,"V_",npl,"$"), main = title)  
      }
    dev.off()
    tools::texi2pdf(file = paste0(plot.path,"vol-factor-",npl,".tex"), clean = TRUE, quiet = TRUE)
  }
  
  # plot pricing error ACFs
  for(npl in 1:N.plots){
    tikzDevice::tikz(file = paste0(plot.path,"fit-err-acf-",npl,".tex"), width = 6, height = 3.7, standAlone = TRUE, sanitize = FALSE)
    loc.data <- data.structure$obs.data[,npl+1] - fit[,npl]
      layout(t(c(1:2)))
      par(mar = c(4.1,4.1,5.1,2.1))
      title <- paste0("Pricing error \n of ",as.character(data.structure$spec.mat$type[npl])," swap rate \n $\\tau=$",sprintf("%1.2f",data.structure$spec.mat$t[npl]),", $p=$", sprintf("%1.1f",data.structure$spec.mat$p[npl]))
      acf(loc.data, demean = T, main = title, lwd = 2)
      title <- paste0("Absolute pricing error \n of ",as.character(data.structure$spec.mat$type[npl])," swap rate \n $\\tau=$",sprintf("%1.2f",data.structure$spec.mat$t[npl]),", $p=$", sprintf("%1.1f",data.structure$spec.mat$p[npl]))
      acf(abs(loc.data), demean = T, main = title, lwd = 2)
    dev.off()
    tools::texi2pdf(file = paste0(plot.path,"fit-err-acf-",npl,".tex"), clean = TRUE, quiet = TRUE)
  }
  
#   # stationary distribution of factors
  vol.paths <- affineSimulate(paramsList = par.list, N.factors = model.spec$N.factors, t.days = 504, t.freq = 1, freq.subdiv = 20, rng.seed = 10984, jumpGeneratorPtr = getPointerToGenerator(model.spec$jump.type), jumpTransformPtr = getPointerToJumpTransform(fstr = model.spec$jump.type)$TF, mod.type = 'standard', nrepl = 2e3)$V.array
#   
#   sim.v <- vol.paths[nrow(vol.paths), seq(2,ncol(vol.paths),by=(1+model.spec$N.factors))]
#   
}

#' @rdname modelEvaluation
#' @export
evaluation_VCOVsandwich <- function(){
  
}