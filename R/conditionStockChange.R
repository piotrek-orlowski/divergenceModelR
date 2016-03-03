#' Return condition list for stock return
#' 
#' @param ret.period The period over which we are interested in the change
#' @export
#' @return A standard cond struct list
condStockChange <- function(ret.period = 5/252, ret.start=0, N.factors = 3) {
  cond.list <- list()
  cond.list[[1]] <- data.frame(t=ret.start+ret.period, u = 1, tau = 1,  p=1)
  cond.list[[2]] <- data.frame(t=ret.start, u = 1, tau = 1, p=1)
  cond.coeff <- c(1,-1)
  return(list(cond.list = lapply(cond.list, function(x) condRegularize(x,N.factors)), cond.coeff = cond.coeff, linTrans = rbind(cond.coeff,c(0,1))))
}