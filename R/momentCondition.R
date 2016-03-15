#' This function calculates the general iterative moment condition that can be used for moment based model calibration. Futures prices (with some long expiry) are assumed for the stock!
#' @param params.P a list of structures that describe the P behaviour of the volatility factors + underlying jumps
#' @param params.Q a list of structures that describe the Q behaviour of the volatility factors + underlying jumps
#' @param condition.struct a list containing the data frames of relative timings, the frequency + maturity pairs and the required powers  / combination coefficients / linear transformation matrix 
#' @param atol absolute tolerance required of the solution (this should be set very low in general)
#' @param rtol relative tolerance required of the solution (the number of correct digits we require)
#' @param ... further arguments passed to \code{\link{twoFactorJumpODEsSolveP}}
#' @export
#' @return Returns the value (unconditional P expectation) of the specified product.
momentCondition <- function(params.P, params.Q, condition.struct, conditional=FALSE, atol = 1e-30, rtol=1e-12, rtol.Q=rtol, N.factors = (length(params.P)-1), jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, mod.type = 'standard',...) {
  # make sure constants can be left out
  stopifnot(abs(sum(condition.struct$linTrans[1,])) < 1e-12)
  
  # first make sure the maturities are aligned for all condition lists
  N <- length(condition.struct$cond.list)
  
  a <- matrix(0,nrow=N,ncol = 1)
  b <- matrix(0,nrow=N,ncol = N.factors + 1)
  
  # now create market data frame, we assume everything is "futurized", so dividend yields & interest rates will be set to 0
  mkt <- data.frame(p=1, q=0, r=0, t = 1)
  
  # now fill up all condition lists with the same maturities
  for (nc in 1:N) {    
    c.list <- condition.struct$cond.list[[nc]]
    
    # make sure conditional list is ordered
    c.list <- c.list[order(c.list$t),]
    
    # now get times to maturity
    t.all <- unique(c.list$t)
    
    # now add initial condition, which is either just before the current time point (referenced at t=0), or way back (unconditional case)
    if (!conditional) {
      # now add the finalizing condition
      t.all <- c(t.all[1]-15,t.all)
    } else {
      # conditionality is always assumed to be at time 0!!! Some structures give back t.all that doesn't have 0, so we have to set an absolute reference point.
      t.all <- sort(unique(c(0,t.all)))
    }
    
    # now iterate backwards
    for (tt in (t.all[order(t.all,decreasing=TRUE)])) {
      cond.t <- subset(c.list,t==tt)
      
      a.loc <- 0
      b.loc <- matrix(0,nrow=1,ncol=1+N.factors)
      
      for (nn in 1:nrow(cond.t)) {
        # if required frequency is not 0 or 1, then we need to calculate the implied c.f.
        if (!(cond.t[nn,"u"] %in% c(0,1))) { 
          mkt$t <- cond.t[nn,"tau"]
          # calculate implied l.t.
          q <- Re(as.vector(jumpDiffusionODEs(matrix(c(cond.t[nn,"u"],rep(0,N.factors)),1,N.factors+1),params=params.Q,mkt=mkt,atol=atol,rtol=rtol.Q, N.factors = N.factors,jumpTransform = jumpTransform, mod.type = mod.type,...)))
          # update constants
          a.loc <- a.loc + q[N.factors+1] * cond.t[nn,"p"]
          
          # update state dependent part
          b.loc[1+1:N.factors] <-  b.loc[1+1:N.factors] + q[1:N.factors] * cond.t[nn,"p"]
        }
        b.loc[1] <- b.loc[1] + cond.t[nn,"u"]*cond.t[nn,"p"]
        for (vv in 1:N.factors) {
          b.loc[1+vv] <- b.loc[1+vv] + cond.t[nn,paste0("v",vv)]
        }
      }
      
      # now update the a and b values
      a[nc,1] <- a[nc,1] + a.loc
      b[nc,] <- b[nc,] + b.loc
      
      if (tt > min(t.all)) {
        # now calculate next maturity (remember, we are going backwards) and update using the P ODE solver
        t.next <- t.all[which(t.all==tt)-1]
        # look at how much the time step is
        mkt$t <- tt-t.next
        
        # now calculate the P ODE solution
        q <- jumpDiffusionODEsP(b[nc,,drop=FALSE],params.P,params.Q,mkt,atol=atol,rtol=rtol,N.factors = N.factors,jumpTransform = jumpTransform, mod.type = mod.type,...)
        
        # update a and b
        a[nc,1] <- a[nc,1] + q[1,1,"a"]
        for (nn in 1:N.factors) {
          b[nc,1+nn] <- q[1,1,paste0("b",nn)]
        }
      }
    }
  }
  
  a <- Re(a)
  b <- Re(b)
  
  if (!conditional) {
    return(linCombMean(a,b[,-1,drop=FALSE],condition.struct$linTrans[1,],rep(1,N.factors)))
  } else {
    return(list(a=a,b=b[,-1,drop=FALSE],linCoeff = condition.struct$linTrans[1,]))
  }
}