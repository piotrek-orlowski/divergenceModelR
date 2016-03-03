#' Return condition list for approximate volatility change
#' 
#' @param factorNum Integer, specifying which factor we are talking about
#' @param ret.period The period over which we are interested in the change
#' @param lin.u A small number such that expm(lin.u * V) is nearly equal to V
#' @export
#' @return A standard cond struct list
condVolChange <- function(factorNum=1, ret.period = 5/252, ret.start = 0, lin.u = 1e-3, N.factors=3) {
  cond.list <- list()

  cond.list[[1]] <- data.frame(t=c(0,ret.period) + ret.start)
  cond.list[[1]][,paste0("v",factorNum)] <- c(-1,1) * lin.u

  cond.list[[2]] <- data.frame(t=0, u = 0, p=1)
  
  cond.coeff <- c(1,-1) / lin.u
  return(list(cond.list = lapply(cond.list, function(x) condRegularize(x,N.factors)), cond.coeff = cond.coeff, linTrans = rbind(cond.coeff,c(0,1/lin.u))))
}