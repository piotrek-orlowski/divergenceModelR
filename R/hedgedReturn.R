#' Create condition structure for hedged option return
#' @param u The frequency corresponding to the option portfolio
#' @param tau.j Start of the return period
#' @param delta Length of the return period
#' @param ttm The time-to-maturity of the option at initiation
#' @export
#' @return A usual structure list
hedged.return <- function(u,tau.j,delta,ttm, N.factors=3){
  cond.list <- list()
  cond.list[[1]] <- data.frame(t=c(0,tau.j+delta), u = u, tau = c(ttm,ttm-tau.j-delta), p = c(-1,1))
  cond.list[[2]] <- data.frame(t=c(0,tau.j), u = u, tau = c(ttm,ttm-tau.j), p=c(-1,1))
  cond.list[[3]] <- data.frame(t=c(0,tau.j,tau.j,tau.j+delta), u = c(u,1,u,1), tau = c(ttm,1,ttm-tau.j,1), p = c(-1,-1,1,1))
  cond.coeff <- c(1,-1+u,-u)
  return(list(cond.list = lapply(cond.list, function(x) condRegularize(x,N.factors)), cond.coeff = cond.coeff, linTrans = rbind(c(1,-1+u,-u),c(0,1,-1+u),c(0,0,1))))
}