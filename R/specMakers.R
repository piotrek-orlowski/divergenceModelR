#' @title Model setup -- specification generators
#' @name modSetup
#' @description Functions which fix model and data specifications for estimation. These functions are not intended for use for the public. They might get separated from the package in the future and placed in a separate archive.
#' @return \code{model.spec} and \code{data.structure} objects

#' @details \code{spec_1FtoyModel} Specifies a single-factor model with leverage effect, exponential jumps in volatility and correlated double exponential jumps in the underlying, with same tail parameter on both sides. 
#' @export
spec_1FtoyModel <- function(){
  
  N.factors <- 1
  model.spec <- model_makeDefaultParameterStructures(N.factors = N.factors, pq.equality = c("Q$jmp$lvec",paste0("Q$jmp$lprop.", 1:N.factors)))
  model.spec$par.names <- c(model.spec$par.names[1], as.character(model.spec$par.restr$par.name[1]), model.spec$par.names[-1])
  model.spec$par.restr <- model.spec$par.restr[-1,]
  model.spec$N.factors <- N.factors
  model.spec$jump.type <- 'kouExpJumpTransform'
  
  return(model.spec)
}

#' @describeIn modSetup
#' @details \code{specData_DSQ_1M_6M_0115} returns loads the divergence sample from 2001 to 2015, p= 0.5, maturities 1/12 and 1/2 years. Includes SP500 returns from OMTR and the bootstrapped error varcov matrices.
#' @export
specData_DSQ_1M_6M_0115 <- function(path.to.data){
  
  load(paste0(path.to.data,"divergence-prices-omtr-asymptotic-iv-for-estimation.RData"))
  load(paste0(path.to.data,"SP500_daily.RData"))
  
  pVec <- 0.5
  tVec <- c(1/12,1/2)
  
  mkt <- data.frame(p=1,q=0,r=0,t=tVec)
  
  require(dplyr)
  
  date.subset <- which(unique(db$day) >= as.Date("2001-01-01"))
  ret.start.date <- unique(db$day)[date.subset[1]-1]
  
  data.spec <- expand.grid(p= unique(db$u), t= unique(db$t), type = c('div','skew','quart'))
  obs.spec <- expand.grid(p= pVec, t= tVec, type = c('div','skew','quart'))
  
  db <- db %>% filter(day >= as.Date("2001-01-01")) %>% filter(u %in% pVec) %>% filter(t %in% tVec)
  db <- db %>% mutate(div = div/t)
  
  obs.data <- db %>% filter(t == 1/12) %>% select(day,div,skew_scaled,quart_scaled)
  obs.data <- cbind(obs.data, db %>% filter(t == 1/2) %>% select(div,skew_scaled,quart_scaled))
  colnames(obs.data) <- c("day","div_0.083_0.5","skew_0.083_0.5","quart_0.083_0.5","div_0.5_0.5","skew_0.5_0.5","quart_0.5_0.5")
  obs.data <- obs.data %>% select(day, div_0.083_0.5, div_0.5_0.5, skew_0.083_0.5, skew_0.5_0.5, quart_0.083_0.5, quart_0.5_0.5)
  
  SP500_daily <- SP500_daily %>% filter(date >= ret.start.date, date <= max(obs.data$day))
  date.ranges <- cbind.data.frame(start= c(ret.start.date,head(obs.data$day,-1)), end = obs.data$day)
  SP500_weekly <- apply(X = date.ranges, MARGIN = 1, FUN = function(x){
    loc.sp <- SP500_daily %>% filter(date > x["start"], date <= x["end"])
    ret.df <- data.frame(date = x["end"], ret = prod(1 + loc.sp$ret)-1)
    return(ret.df)
  })
  SP500_weekly <- rbind_all(SP500_weekly)
  SP500_weekly <- SP500_weekly %>% mutate(date = as.Date(date))  
  obs.data <- SP500_weekly %>% inner_join(obs.data, by = c("date"="day"))
  
  data.spec <- cbind(data.spec,index = 1:nrow(data.spec))
  obs.spec <- data.spec %>% inner_join(obs.spec)
  
  noise.cov.cube <- db.vcov[obs.spec$index,obs.spec$index,date.subset]
  noise.cov.cube[which(obs.spec$type=='div'),which(obs.spec$type=='div'),] <- noise.cov.cube[which(obs.spec$type=='div'),which(obs.spec$type=='div'),] * array(outer(1/tVec, 1/tVec), dim = c(length(tVec), length(tVec), dim(noise.cov.cube)[3]))
  dt <- 7/365
  
  spec.mat <- obs.spec %>% select(-index)
  return(list(dt = dt, noise.cov.cube = noise.cov.cube, spec.mat = spec.mat, obs.data = as.matrix(obs.data[,-1]), dates = obs.data[,1], mkt = mkt))
}

#' @describeIn modSetup
#' @details \code{spec_3FsepIntModel} Specifies a three-factor model with leverage effect, exponential jumps in volatility and correlated double exponential jumps in the underlying, with same tail parameter on both sides. One of the factors only drives jump intensity, not the continuous volatility
#' @export
spec_3FsepIntModel <- function(){
  
  N.factors <- 3
  model.spec <- model_makeDefaultParameterStructures(N.factors = N.factors, pq.equality = c("Q$jmp$lvec",paste0("Q$jmp$lprop.", 1:(N.factors-1)),"Q$jmp$muYc",paste0("Q$",2,"$eta")))
  model.spec$par.names <- c(model.spec$par.names[1], as.character(model.spec$par.restr$par.name[1]), model.spec$par.names[-1])
  model.spec$par.restr <- model.spec$par.restr[-1,]
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.jmp.muYc",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$jmp$muYc", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.3.phi",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$3$phi", par.value = 0),model.spec$par.restr)
  model.spec$N.factors <- N.factors
  model.spec$jump.type <- 'kouExpJumpTransform'
  
  return(model.spec)  
}