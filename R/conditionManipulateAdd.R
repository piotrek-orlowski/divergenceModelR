#' This function takes two data frame lists which define elements of moment conditions, and two vectors of corresponding linear combination coefficients, and combines them into a single list of data frames which would result from multiplying the two expressions rational in GLT-portfolios.
#' @param cond.list.1 A list of data frames defining elements of the moment condition.
#' @param cond.list.2 A list of data frames defining elements of the moment condition.
#' @param coeffs.1 A list of coefficients for taking a linear combination of the individual moment conditions defined by \code{cond.list.1}; must be of equal length as \code{cond.list.1}
#' @param coeffs.2 A list of coefficients for taking a linear combination of the individual moment conditions defined by \code{cond.list.2}; must be of equal length as \code{cond.list.2}
#' @param lag Positive real, the time distance between elements of \code{cond.list.1} and \code{cond.list.2}, the former happening ``later''.
#' @export
#' @return List with named elements, \code{cond.list} is a list akin to the input lists, \code{cond.coeff} is a list of corresponding coefficients.
#'

condAdd <- function(condStruct.1, condStruct.2, weights = c(1,1)) {
  coeff.add <- NULL
  cond.add <- list()
  linTrans.add <- NULL
  
  # check the trivial cases first
  if (length(condStruct.1$cond.list)==0 & length(condStruct.2$cond.list) == 0) {
    return(list(cond.list = cond.add, cond.coeff = coeff.add, linTrans = linTrans.add))
  }
  if (length(condStruct.1$cond.list)==0 & length(condStruct.2$cond.list) >0) {
    return(list(cond.list = condStruct.2$cond.list, cond.coeff = condStruct.2$cond.coeff*weights[2], linTrans = condStruct.2$linTrans*weights[2]))
  }
  if (length(condStruct.1$cond.list) > 0 & length(condStruct.2$cond.list) ==0) {
    return(list(cond.list = condStruct.1$cond.list, cond.coeff = condStruct.1$cond.coeff*weights[1], linTrans = condStruct.1$linTrans*weights[1]))
  }
  
  # if we haven't returned yet, then both lists are not null. In this case combine.
  coeff.add <- c(weights[1] * condStruct.1$cond.coeff, weights[2] * condStruct.2$cond.coeff)
  cond.add <- c(condStruct.1$cond.list, condStruct.2$cond.list)
  
  # now combine the linear transformation matrices
  linTrans.add <- condStruct.1$linTrans * weights[1]
  
  for (ii in 1:length(condStruct.1$cond.list)) {
    linTrans.add <- rbind(linTrans.add,0)
  }
  
  for (ii in 1:length(condStruct.2$cond.list)) {
    linTrans.add <- cbind(linTrans.add,weights[2]*c(rep(condStruct.2$linTrans[1,ii],length(condStruct.1$cond.list)),condStruct.2$linTrans[,ii]))
  }
  
  return(list(cond.list = cond.add, cond.coeff = coeff.add, linTrans = linTrans.add))
}