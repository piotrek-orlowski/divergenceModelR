#' @title Integrated variance risk premia
#' @description Given the parameters and initial values, this function produces the term structure of the integrated variance risk premium. Definition of IVRP as in reference below. This function treats the case of no volatility jumps if \code{control$use.CIR = TRUE}.
#' @param params.structure a named list of named lists that describe the \code{P} and \code{Q} behaviour of the stochastic volatility model.
#' @param v.0 initial condition for \code{IVRP} calculation, if \code{NULL}, long-run factor mean values are used instead.
#' @param tT vector of length at least 2, start- and end-point of calculation. If the vector is lonver, the VRP will be returned at the intermediate values as well. Otherwise, the output values are determined via \code{rp.control$grid.d}.
#' @param grid.d determines the number of intervals into which the \code{tT} period is divided
#' @export
#' @return Returns a list with fields \code{tot}, \code{cont} and \code{jmp}, referring to, respectively, the total risk premium, its continuous and discontinuous parts.
#' @references Ait-Sahalia, Y. and Karaman, M., and Mancini, L., ``The Term Structure of Variance Risk Premia and the Expectation Hypothesis''

integrVarRP <- function(params.structure, v.0 = NULL, tT = c(0,1), grid.d = 1e3, N.factors = 3){
  
  params.vector <- NULL
  
  params.P <- params.structure$P
  params.Q <- params.structure$Q
  
  if(is.null(params.structure) & is.null(params.vector)){
    stop("No parameters were provided.")
  }
  
  if(any(tT < 0)){
    stop("Negative time provided in tT.")
  }
  
  if(is.null(v.0) & !is.null(params.structure)){
    v.0 <- numeric()
    for(kk in 1:N.factors){
      v.0[kk] <- params.P[[as.character(kk)]]$eta
    }
  }
  ## To calculate the risk premia, you will have to calculate the second moment of the jumps, as well as the expected values of the variance factors under P and Q for the range of times from 0 to t.
  
  J.2.P <- (params.P$jmp$muYc + params.P$jmp$muSc * params.P$jmp$rhoc)^2 + params.P$jmp$sigmaYc^2
  J.2.Q <- (params.Q$jmp$muYc + params.Q$jmp$muSc * params.Q$jmp$rhoc)^2 + params.Q$jmp$sigmaYc^2
  
  
  locWrapper <- function(ttm,v.nought,sv){
    return(.Call("cirIntExp",c(sv$kpp,sv$eta),ttm[1],ttm[2],v.nought))
  }
  
  mkt <- list(p = 1, q = 0, r = 0, t = 0)
  mkt$t <- seq(tT[1], tT[2], length.out = grid.d)
  mkt$t <- c(mkt$t, tT[-c(1,length(tT))])
  mkt$t <- mkt$t[order(mkt$t)]
  
  ttm.mat <- t(apply(X = matrix(2:length(mkt$t),length(mkt$t)-1,1), MARGIN = 1, FUN = function(x){
    return(c(mkt$t[x-1],mkt$t[x]))
  }))
  
  ## calculate integrated moments with the use of the CPP-based function for those factors that don't jump. If P$jmp$muSc == Q$jmp$muSc == 1e-6 or less, assume that the first factor doesn't jump either.
  
  if(params.P$jmp$muSc == params.Q$jmp$muSc & params.P$jmp$muSc < 1e-6){
    start.at <- 1
  }else{
    start.at <- 2
  }
  
  v.mom.int.Q <- v.mom.int.P <- matrix(0,nrow = nrow(ttm.mat), ncol = N.factors)
  for(nn in start.at:N.factors){
    v.mom.int.P[,nn]  <- apply(X = ttm.mat, MARGIN = 1, FUN = locWrapper, v.nought = v.0[nn], sv = params.P[[as.character(nn)]])
    v.mom.int.P[,nn] <- cumsum(v.mom.int.P[,nn])
    v.mom.int.Q[,nn]  <- apply(X = ttm.mat, MARGIN = 1, FUN = locWrapper, v.nought = v.0[nn], sv = params.Q[[as.character(nn)]])
    v.mom.int.Q[,nn] <- cumsum(v.mom.int.Q[,nn])
  }
  
  # 
  if(start.at > 1){
    ## Prepare CF evaluator
    cfDerivAtZero <- function(measure){
      h.step <- 1e-3
      loc.u <- rep(0,N.factors+1)
      loc.u[2] <- h.step
      mkt.loc <- mkt
      mkt.loc$t <- mkt.loc$t[-1]
      if(measure == "P"){
        loc.sol <- twoFactorJumpODEsSolveP(matrix(loc.u,nrow=1,byrow=TRUE), params.P, params.Q, mkt.loc, N.factors = N.factors, rtol = 1e-15) 
      }
      else if(measure == "M"){
        params.M <- params.Q
        params.M$jmp <- params.P$jmp
        loc.sol <- twoFactorJumpODEsSolveP(matrix(loc.u,nrow=1,byrow=TRUE), params.M, params.Q, mkt.loc, N.factors = N.factors, rtol = 1e-15) 
      }
      else{
        loc.sol <- twoFactorJumpODEsSolve(matrix(loc.u,nrow=1,byrow=TRUE), params.Q, mkt.loc, N.factors = N.factors, rtol = 1e-15)
      }
      res <- 1/h.step * expm1(Re(loc.sol[,,"a"] + loc.sol[,,"b1"]*v.0[1]))
#       res <- loc.sol[,,"a"]
#       for(kk in 1:N.factors){
#         res <- res + v.0[kk] * loc.sol[,,paste0("b",kk)]
#       }
#       res <- 1/h.step * expm1(Re(res))
      return(res)
    }
    v.mom.int.P[,1] <- cfDerivAtZero("P")
    v.mom.int.P[,1] <- cumsum(v.mom.int.P[,1] * (tT[2] - tT[1])/grid.d) # integration
    v.mom.int.Q[,1] <- cfDerivAtZero("Q")
    v.mom.int.Q[,1] <- cumsum(v.mom.int.Q[,1] * (tT[2] - tT[1])/grid.d) # integration
    v.mom.int.M <- cfDerivAtZero("M") # for jump factor rp decomposition
    v.mom.int.M <- cumsum(v.mom.int.M * (tT[2] - tT[1])/grid.d) # integration
  }
  
  phi.vec <- vapply(params.P[as.character(1:N.factors)],FUN.VALUE = numeric(1), function(x) return(x$phi))

  ## Initialise RP return list
  ivrp <- list()
  ivrp$ttm <- mkt$t[-1]
  ivrp$jmpP <- (params.P$jmp$lvec[1] * ivrp$ttm)
  ivrp$jmpP <- ivrp$jmpP + rowSums(t(apply(v.mom.int.P, 1, function(x) return(x * params.P$jmp$lprop))))
  ivrp$jmpP <- ivrp$jmpP * J.2.P
  ivrp$jmpQ <- (params.Q$jmp$lvec[1] * ivrp$ttm)
  ivrp$jmpQ <- ivrp$jmpQ + rowSums(t(apply(v.mom.int.Q, 1, function(x) return(x * params.Q$jmp$lprop))))
  ivrp$jmpQ <- ivrp$jmpQ * J.2.Q
  
  ivrp$P <- ivrp$jmpP + rowSums(v.mom.int.P * matrix(phi.vec^2,nrow = nrow(v.mom.int.P), ncol = N.factors, byrow=T))
  ivrp$Q <- ivrp$jmpQ + rowSums(v.mom.int.Q * matrix(phi.vec^2,nrow = nrow(v.mom.int.Q), ncol = N.factors, byrow=T))

  ivrp$pureInt <- params.P$jmp$lvec[1] * ivrp$ttm - params.Q$jmp$lvec[1] * ivrp$ttm
  ivrp$pureInt <- ivrp$pureInt + rowSums(t(apply(v.mom.int.P, 1, function(x) return(x * params.P$jmp$lprop)))) - rowSums(t(apply(v.mom.int.Q, 1, function(x) return(x * params.Q$jmp$lprop))))
  ivrp$pureInt <- ivrp$pureInt * J.2.P / ivrp$ttm

  ivrp$pureJmp <- (params.P$jmp$lvec[1] * ivrp$ttm + rowSums(t(apply(v.mom.int.Q, 1, function(x) return(x * params.Q$jmp$lprop))))) * (J.2.P - J.2.Q) / ivrp$ttm
  
  ivrp$pureJmpByFactor <- ((t(apply(v.mom.int.Q, 1, function(x) return(x * params.Q$jmp$lprop))))) * (J.2.P - J.2.Q) / ivrp$ttm
  
  ivrp$pureIntByFactor <- 0
  ivrp$pureIntByFactor <- ivrp$pureIntByFactor + (t(apply(v.mom.int.P, 1, function(x) return(x * params.P$jmp$lprop)))) - (t(apply(v.mom.int.Q, 1, function(x) return(x * params.Q$jmp$lprop))))
  ivrp$pureIntByFactor <- ivrp$pureIntByFactor * J.2.P / ivrp$ttm

  ivrp$factors <- (v.mom.int.P - v.mom.int.Q)/ivrp$ttm
  ivrp$varJmpDecomp <- cbind(v.mom.int.P[,1] - v.mom.int.M, v.mom.int.M - v.mom.int.Q[,1])/ivrp$ttm
  
  if(length(tT) > 2){
    loc.tT <- tT[-c(1,length(tT))]
    ivrp$cont <- ivrp$cont[which(loc.tT %in% mkt$t[-1])]
    ivrp$jmp <- ivrp$jmp[which(loc.tT %in% mkt$t[-1])]
    ivrp$ttm <- mkt$t[which(loc.tT %in% mkt$t[-1])]
  }
  
  # the return on the variance-swap
  ivrp$ret <- (ivrp$P-ivrp$Q)/(ivrp$Q)
  
  # the "dollar-value" of the VRP
  ivrp$dollar <- (ivrp$P-ivrp$Q)/ivrp$ttm
  
  # how much we would earn on a "jump-swap"
  ivrp$jmp <- (ivrp$jmpP-ivrp$jmpQ)/ivrp$jmpQ
  
  # how much we would earn on a "diffusive-swap"
  ivrp$cont <- (rowSums(v.mom.int.P * matrix(phi.vec^2,nrow = nrow(v.mom.int.P), ncol = N.factors, byrow=T)) - rowSums(v.mom.int.Q * matrix(phi.vec^2,nrow = nrow(v.mom.int.Q), ncol = N.factors, byrow=T)))#/(rowSums(v.mom.int.Q * matrix(phi.vec^2,nrow = nrow(v.mom.int.Q), ncol = N.factors, byrow=T)))
  
  return(ivrp)
}