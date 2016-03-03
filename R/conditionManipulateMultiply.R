#' This function takes two data frame lists which define elements of moment conditions, and two vectors of corresponding linear combination coefficients, and combines them into a single list of data frames which would result from multiplying the two expressions rational in GLT-portfolios.
#' @param cond.list.1 A list of data frames defining elements of the moment condition.
#' @param cond.list.2 A list of data frames defining elements of the moment condition.
#' @param coeffs.1 A list of coefficients for taking a linear combination of the individual moment conditions defined by \code{cond.list.1}; must be of equal length as \code{cond.list.1}
#' @param coeffs.2 A list of coefficients for taking a linear combination of the individual moment conditions defined by \code{cond.list.2}; must be of equal length as \code{cond.list.2}
#' @param lag Positive real, the time distance between elements of \code{cond.list.1} and \code{cond.list.2}, the former happening ``later''.
#' @export
#' @return List with named elements, \code{cond.list} is a list akin to the input lists, \code{cond.coeff} is a list of corresponding coefficients.
#'

condMultiply <- function(condStruct.1,condStruct.2, lag = 0, compactify = FALSE, N.factors = 3) {
  cond.list.1 <- condStruct.1$cond.list
  cond.list.2 <- condStruct.2$cond.list
  coeffs.1 <- condStruct.1$cond.coeff
  coeffs.2 <- condStruct.2$cond.coeff
  linTrans.1 <- condStruct.1$linTrans
  linTrans.2 <- condStruct.2$linTrans
   
  coeff.mult <- NULL
  cond.mult <- list()
  linTrans.mult <- NULL
  
  i <- 1
  for (l1 in 1:length(cond.list.1)) {
    for (l2 in 1:length(cond.list.2)) {
      coeff.mult <- c(coeff.mult, coeffs.1[l1]*coeffs.2[l2])
      c.1 <- condRegularize(cond.list.1[[l1]], N.factors)
      c.1$t <- c.1$t + lag
      c.2 <- condRegularize(cond.list.2[[l2]], N.factors)
      list.added <- rbind(c.2,c.1)
      
      # now make sure all frames have the same dimnames (will be required for comparison)    
      cond.mult[[i]] <- condRegularize(list.added, N.factors)
      
      i <- i+1
    }
  }
  linTrans.mult <- kronecker(linTrans.1,linTrans.2)
  if (compactify & nrow(as.matrix(linTrans.mult))>1) {
    return(condCompactify(cond.mult,coeff.mult,linTrans.mult))
  } else {
    return(list(cond.list = cond.mult, cond.coeff = coeff.mult, linTrans = linTrans.mult))
  }
}