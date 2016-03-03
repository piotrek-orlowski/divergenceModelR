#' This function takes a list of moment-condition defining data frames and a vector of corresponding coefficients, equal in length to the list. The data frame list is searched for identical entries, which are made unique, and the corresponding coefficients are summed.
#' @param cond.list A list of data frames defining elements of the moment condition.
#' @param cond.coeff A list of coefficients for taking a linear combination of the individual moment conditions defined by \code{cond.list}; must be of equal length as \code{cond.list}
#' @export
#' @return List with named elements, \code{cond.list} is a list akin to the input list, \code{cond.coeff} is a list of corresponding coefficients.
#'

condCompactify <- function(cond.list, cond.coeff, cond.linTrans) {
  cond.compact <- list()
  coeff.compact <- NULL
  linTrans.compact <- cond.linTrans
  n.elements <- 1
  
  stopifnot(length(cond.list) == length(cond.coeff) & length(cond.list) == nrow(cond.linTrans))
  
  # we will use this matrix to reduce the linear transformation matrix
  reductorMat <- diag(rep(1,nrow(cond.linTrans)))
  
  for (nn in seq(length(cond.list),1,-1)) {
    c.list  <- cond.list[[nn]]
    # if the current panel is not yet in the compactified list then add & count the number of occurences
    if (!any(sapply(cond.compact, function(x) return(identical(x,c.list))))) {
      cond.compact[[n.elements]] <- c.list
      n.elements <- n.elements + 1
      #coeff.compact <- c(coeff.compact,sum(cond.coeff[sapply(cond.list,function(x) return(identical(x,c.list)))]))
      
      redundant.variables <- sapply(cond.list,function(x) return(identical(x,c.list)))
      redundant.variables[nn] <- FALSE
      diag(reductorMat) <- diag(reductorMat) * (!redundant.variables)
      reductorMat[redundant.variables,nn]<- 1
    }
  }
  # we need to reverse the elements
  cond.out <- cond.compact
  for (nn in 1:length(cond.compact)) {cond.out[[length(cond.compact)+1-nn]] <- cond.compact[[nn]]}
  linTrans.reduced <- as.matrix(cond.linTrans %*% reductorMat)
  linTrans.reduced <- linTrans.reduced[,colSums(abs(linTrans.reduced))>0]

  linTrans.lu <- .Call("luComplex",list(rMat = Re(linTrans.reduced), iMat = Im(linTrans.reduced)))
  
  if (max(abs(linTrans.lu$iU))>1e-6) {warning("complex values not supported!")}
  #linTrans.lu <- (linTrans.lu$rU + 1i*linTrans.lu$iU)
  linTrans.lu <- linTrans.lu$rU
  
  linTrans.out <- as.matrix(linTrans.lu)
  return(list(cond.list = cond.out, cond.coeff = linTrans.out[1,], linTrans = linTrans.out))
}