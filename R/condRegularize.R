#' Ensure that condition lists have uniform columns
#' 
#' This function takes a condition list and fills it up the missing columns with 0s.
#' @param cond.struct The condition structure to be regularized
#' @param N.factors The number of volatility factors assumed in the calculations.
#' @export
#' @return Returns a regularized condition list.
condRegularize <- function(cond.struct,N.factors = 3) {
  col.names <- c("t","u","tau","p",paste("v",seq(1,N.factors),sep=""))
  
  # times ALWAYS have to be specified
  stopifnot("t" %in% names(cond.struct))
  
  # now missing fields should be filled with 0-s
  cond.struct[,setdiff(col.names, names(cond.struct))] <- 0
  
  # now let's have the proper column order, although this is not strictly relevant (R also takes care of this)
  cond.struct <- cond.struct[,col.names]
  
  # sort the rows according to time and add reference 0 timepoint
  cond.reg <- rbind(0,cond.struct[order(cond.struct$t),])
  
  # remove row names to ensure we can use the identical function
  rownames(cond.reg) <- NULL
  
  return(cond.reg)
}