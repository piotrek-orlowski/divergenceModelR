#' @title Model setup -- specification generators
#' @name modSetup
#' @description Functions which fix model and data specifications for estimation. These functions are not intended for use for the public. They might get separated from the package in the future and placed in a separate archive.
#' @return \code{model.spec} and \code{data.structure} objects

#' @details \code{specData_DSQ_1M_6M_0115_extraNoise} returns loads the divergence sample from 2001 to 2015, p= 0.5, maturities 1/12 and 1/2 years. Includes SP500 returns from OMTR and the bootstrapped error correlation matrices.
#' @export
specData_DSQ_1M_6M_0115_extraNoise <- function(path.to.data){
  
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
  SP500_weekly <- bind_rows(SP500_weekly)
  SP500_weekly <- SP500_weekly %>% mutate(date = as.Date(date))  
  obs.data <- SP500_weekly %>% inner_join(obs.data, by = c("date"="day"))
  
  data.spec <- cbind(data.spec,index = 1:nrow(data.spec))
  obs.spec <- data.spec %>% inner_join(obs.spec)
  
  noise.cov.cube <- db.vcov[obs.spec$index,obs.spec$index,date.subset]
  noise.cov.cube[which(obs.spec$type=='div'),which(obs.spec$type=='div'),] <- noise.cov.cube[which(obs.spec$type=='div'),which(obs.spec$type=='div'),] * array(outer(1/tVec, 1/tVec), dim = c(length(tVec), length(tVec), dim(noise.cov.cube)[3]))
  dt <- 7/365
  
  noise.cov.cube <- apply(noise.cov.cube,3,cov2cor)
  noise.cov.cube <- array(as.numeric(noise.cov.cube), dim = c(ncol(obs.data)-2,ncol(obs.data)-2,nrow(obs.data)))
  
  spec.mat <- obs.spec %>% select(-index)
  return(list(dt = dt, noise.cov.cube = noise.cov.cube, spec.mat = spec.mat, obs.data = as.matrix(obs.data[,-1]), dates = obs.data[,1], mkt = mkt))
}

#' @rdname modSetup
#' @details \code{specData_DS_1M_6M_0115_extraNoise} returns loads the divergence sample from 2001 to 2015, p= 0.5, maturities 1/12 and 1/2 years. Includes SP500 returns from OMTR and the bootstrapped error correlation matrices. Excludes quarticity
#' @export
specData_DS_1M_6M_0115_extraNoise <- function(path.to.data){
  
  load(paste0(path.to.data,"divergence-prices-omtr-asymptotic-iv-for-estimation.RData"))
  load(paste0(path.to.data,"SP500_daily.RData"))
  
  pVec <- 0.5
  tVec <- c(1/12,1/2)
  
  mkt <- data.frame(p=1,q=0,r=0,t=tVec)
  
  require(dplyr)
  
  date.subset <- which(unique(db$day) >= as.Date("2001-01-01"))
  ret.start.date <- unique(db$day)[date.subset[1]-1]
  
  data.spec <- expand.grid(p= unique(db$u), t= unique(db$t), type = c('div','skew','quart'))
  obs.spec <- expand.grid(p= pVec, t= tVec, type = c('div','skew'))
  
  db <- db %>% filter(day >= as.Date("2001-01-01")) %>% filter(u %in% pVec) %>% filter(t %in% tVec)
  db <- db %>% mutate(div = div/t)
  
  obs.data <- db %>% filter(t == 1/12) %>% select(day,div,skew_scaled)
  obs.data <- cbind(obs.data, db %>% filter(t == 1/2) %>% select(div,skew_scaled))
  colnames(obs.data) <- c("day","div_0.083_0.5","skew_0.083_0.5","div_0.5_0.5","skew_0.5_0.5")
  obs.data <- obs.data %>% select(day, div_0.083_0.5, div_0.5_0.5, skew_0.083_0.5, skew_0.5_0.5)
  
  SP500_daily <- SP500_daily %>% filter(date >= ret.start.date, date <= max(obs.data$day))
  date.ranges <- cbind.data.frame(start= c(ret.start.date,head(obs.data$day,-1)), end = obs.data$day)
  SP500_weekly <- apply(X = date.ranges, MARGIN = 1, FUN = function(x){
    loc.sp <- SP500_daily %>% filter(date > x["start"], date <= x["end"])
    ret.df <- data.frame(date = x["end"], ret = prod(1 + loc.sp$ret)-1)
    return(ret.df)
  })
  SP500_weekly <- bind_rows(SP500_weekly)
  SP500_weekly <- SP500_weekly %>% mutate(date = as.Date(date))  
  obs.data <- SP500_weekly %>% inner_join(obs.data, by = c("date"="day"))
  
  data.spec <- cbind(data.spec,index = 1:nrow(data.spec))
  obs.spec <- data.spec %>% inner_join(obs.spec)
  
  noise.cov.cube <- db.vcov[obs.spec$index,obs.spec$index,date.subset]
  noise.cov.cube[which(obs.spec$type=='div'),which(obs.spec$type=='div'),] <- noise.cov.cube[which(obs.spec$type=='div'),which(obs.spec$type=='div'),] * array(outer(1/tVec, 1/tVec), dim = c(length(tVec), length(tVec), dim(noise.cov.cube)[3]))
  dt <- 7/365
  
  noise.cov.cube <- apply(noise.cov.cube,3,cov2cor)
  noise.cov.cube <- array(as.numeric(noise.cov.cube), dim = c(ncol(obs.data)-2,ncol(obs.data)-2,nrow(obs.data)))
  
  spec.mat <- obs.spec %>% select(-index)
  return(list(dt = dt, noise.cov.cube = noise.cov.cube, spec.mat = spec.mat, obs.data = as.matrix(obs.data[,-1]), dates = obs.data[,1], mkt = mkt))
}

#' @rdname modSetup
#' @details \code{specData_D_1M_6M_0115_extraNoise} returns loads the divergence sample from 2001 to 2015, p= 0.5, maturities 1/12 and 1/2 years. Includes SP500 returns from OMTR and the bootstrapped error correlation matrices. Excludes skewness and quarticity
#' @export

specData_D_1M_6M_0115_extraNoise <- function(path.to.data){
  
  load(paste0(path.to.data,"divergence-prices-omtr-asymptotic-iv-for-estimation.RData"))
  load(paste0(path.to.data,"SP500_daily.RData"))
  
  pVec <- 0.5
  tVec <- c(1/12,1/2)
  
  mkt <- data.frame(p=1,q=0,r=0,t=tVec)
  
  require(dplyr)
  
  date.subset <- which(unique(db$day) >= as.Date("2001-01-01"))
  ret.start.date <- unique(db$day)[date.subset[1]-1]
  
  data.spec <- expand.grid(p= unique(db$u), t= unique(db$t), type = c('div','skew','quart'))
  obs.spec <- expand.grid(p= pVec, t= tVec, type = c('div'))
  
  db <- db %>% filter(day >= as.Date("2001-01-01")) %>% filter(u %in% pVec) %>% filter(t %in% tVec)
  db <- db %>% mutate(div = div/t)
  
  obs.data <- db %>% filter(t == 1/12) %>% select(day,div)
  obs.data <- cbind(obs.data, db %>% filter(t == 1/2) %>% select(div))
  colnames(obs.data) <- c("day","div_0.083_0.5","div_0.5_0.5")
  obs.data <- obs.data %>% select(day, div_0.083_0.5, div_0.5_0.5)
  
  SP500_daily <- SP500_daily %>% filter(date >= ret.start.date, date <= max(obs.data$day))
  date.ranges <- cbind.data.frame(start= c(ret.start.date,head(obs.data$day,-1)), end = obs.data$day)
  SP500_weekly <- apply(X = date.ranges, MARGIN = 1, FUN = function(x){
    loc.sp <- SP500_daily %>% filter(date > x["start"], date <= x["end"])
    ret.df <- data.frame(date = x["end"], ret = prod(1 + loc.sp$ret)-1)
    return(ret.df)
  })
  SP500_weekly <- bind_rows(SP500_weekly)
  SP500_weekly <- SP500_weekly %>% mutate(date = as.Date(date))  
  obs.data <- SP500_weekly %>% inner_join(obs.data, by = c("date"="day"))
  
  data.spec <- cbind(data.spec,index = 1:nrow(data.spec))
  obs.spec <- data.spec %>% inner_join(obs.spec)
  
  noise.cov.cube <- db.vcov[obs.spec$index,obs.spec$index,date.subset]
  noise.cov.cube[which(obs.spec$type=='div'),which(obs.spec$type=='div'),] <- noise.cov.cube[which(obs.spec$type=='div'),which(obs.spec$type=='div'),] * array(outer(1/tVec, 1/tVec), dim = c(length(tVec), length(tVec), dim(noise.cov.cube)[3]))
  dt <- 7/365
  
  noise.cov.cube <- apply(noise.cov.cube,3,cov2cor)
  noise.cov.cube <- array(as.numeric(noise.cov.cube), dim = c(ncol(obs.data)-2,ncol(obs.data)-2,nrow(obs.data)))
  
  spec.mat <- obs.spec %>% select(-index)
  return(list(dt = dt, noise.cov.cube = noise.cov.cube, spec.mat = spec.mat, obs.data = as.matrix(obs.data[,-1]), dates = obs.data[,1], mkt = mkt))
}

#' @rdname modSetup
#' @details \code{spec_1FtoyModel_extraNoise}  Specifies a single-factor model with leverage effect, exponential jumps in volatility and correlated double exponential jumps in the underlying, with same tail parameter on both sides. There are variance parameters for observation noise.
#' @export

spec_1FtoyModel_extraNoise <- function(U){
  
  N.factors <- 1
  model.spec <- model_makeDefaultParameterStructures(N.factors = N.factors, pq.equality = c("Q$jmp$lvec",paste0("Q$jmp$lprop.", 1:N.factors)))
  erp.pos <- which(grepl("erp",as.character(model.spec$par.restr$par.name)))
  model.spec$par.names <- c(model.spec$par.names[1], as.character(model.spec$par.restr$par.name[erp.pos]), model.spec$par.names[-1])
  model.spec$par.restr <- model.spec$par.restr[-erp.pos,]
  model.spec$N.factors <- N.factors
  model.spec$jump.type <- 'kouExpJumpTransform'
  
  model.spec$par.names <- model.spec$par.names[order(model.spec$par.names)]
  model.spec$par.restr <- model.spec$par.restr[order(as.character(model.spec$par.restr$par.name)),]
  
  model.spec$scaleFoo <- function(par.vec){
    par.vec[1] <- -0.2 + 0.4*par.vec[1]
    par.vec[2] <- -0.2 + 0.3 * par.vec[2]
    par.vec[3] <- 1e-2 + 20 * par.vec[3]
    par.vec[4] <- 1e-4 + 4 * par.vec[4]
    par.vec[5] <- 1e-4 + 0.5 * par.vec[5]
    par.vec[6] <- -1 + 2*par.vec[6]
    par.vec[7] <- 1e-2 + 20 * par.vec[7]
    par.vec[8] <- 1e-2 + 20 * par.vec[8]
    par.vec[9] <- 1/(1e-4 + par.vec[9])
    par.vec[10] <- -0.2 + 0.3 * par.vec[10]
    par.vec[11] <- -0.2 + 0.4 * par.vec[11]
    par.vec[12] <- 1e-4 + 0.3 * par.vec[12]
    par.vec[13] <- 1e-4 + 4 * par.vec[13]
    par.vec[14] <- 1e-2 + 16 * par.vec[14]
    par.vec[15] <- 1/(1e-4 + par.vec[15])
    par.vec[14] <- -0.15 + 0.3 * par.vec[14]
    par.vec[15] <- -0.2 + 0.4 * par.vec[15]
    par.vec[16] <- 1e-4 + 0.3 * par.vec[16]
    
    for(kk in 19:(19+U-1)){
      par.vec[kk] <- 1e-4 + 1e-1 * par.vec[kk]
    }
    
    return(par.vec)
  }
  
  return(model.spec)
}

#' @rdname modSetup
#' @param U number of observed pfolios
#' @details \code{spec_3FsepIntModel_erp_extraNoise} Specifies a three-factor model with leverage effect, exponential jumps in volatility and correlated double exponential jumps in the underlying, with same tail parameter on both sides. One of the factors only drives jump intensity, not the continuous volatility. There are variance parameters for observation noise. Plus, there is a more flexible parametrisation of the equity risk premium parameters, with some Q parameters restricted. Finally, the first vol factor has pure-jump dynamics
#' @export

spec_3FsepIntModel_erp_extraNoise <- function(U){
  
  N.factors <- 3
  model.spec <- model_makeDefaultParameterStructures(N.factors = N.factors, pq.equality = c("Q$jmp$lvec",paste0("Q$jmp$lprop.",c(2,3)),"Q$jmp$muYc",paste0("Q$",c(1,3),"$eta")))
  model.spec$par.names <- c(model.spec$par.names, "P$2$erp")
  model.spec$par.restr <- model.spec$par.restr[-which(grepl("P.2.erp",model.spec$par.restr$par.name)),]
  model.spec$par.names <- c(model.spec$par.names, "P$1$erp")
  model.spec$par.restr <- model.spec$par.restr[-which(grepl("P.1.erp",model.spec$par.restr$par.name)),]
  model.spec$par.names <- c(model.spec$par.names, "P$3$erp")
  model.spec$par.restr <- model.spec$par.restr[-which(grepl("P.3.erp",model.spec$par.restr$par.name)),]
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.jmp.muYc",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$jmp$muYc", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.3.phi",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$3$phi", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.3.rho",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$3$rho", par.value = 0),model.spec$par.restr)
  # model.spec$par.names <- model.spec$par.names[-which(grepl("P.2.rho",model.spec$par.names))]
  # model.spec$par.restr <- rbind(data.frame(par.name = "P$2$rho", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.jmp.lprop.2",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$jmp$lprop.2", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.1.lmb",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$1$lmb", par.value = 1e-2),model.spec$par.restr)
  model.spec$N.factors <- N.factors
  model.spec$jump.type <- 'kouExpJumpTransform'
  
  model.spec$par.names <- model.spec$par.names[order(model.spec$par.names)]
  model.spec$par.restr <- model.spec$par.restr[order(as.character(model.spec$par.restr$par.name)),]
  
  model.spec$scaleFoo <- function(par.vec){
    
    par.vec[1] <- -0.3 + 0.6*par.vec[1]           # P$1$erp
    par.vec[2] <- -0.3 + 0.6*par.vec[2]           # P$1$erp0
    par.vec[3] <- 1e-2 + 16 * par.vec[3]          # P$1$kpp
    # par.vec[4] <- 1e-4 + 4 * par.vec[4]           # P$1$lmb
    par.vec[4] <- 1e-4 + 0.3 * par.vec[4]         # P$1$phi
    par.vec[5] <- -1 + 2*par.vec[5]               # P$1$rho
    par.vec[6] <- -0.3 + 0.6 * par.vec[6]           # P$2$erp
    par.vec[7] <- 1e-2 + 16 * par.vec[7]          # P$2$kpp
    par.vec[8] <- 1e-4 + 4 * par.vec[8]         # P$2$lmb
    par.vec[9] <- 1e-4 + 0.3 * par.vec[9]         # P$2$phi
    par.vec[10] <- -1 + 2*par.vec[10]               # P$2$rho
    par.vec[11] <- -0.3 + 0.6*par.vec[11]           # P$3$erp
    par.vec[12] <- 1e-2 + 16 * par.vec[12]        # P$3$kpp
    par.vec[13] <- 1e-4 + 4 * par.vec[13]         # P$3$lmb
    par.vec[14] <- 1e-4 + 20 * par.vec[14]        # P$jmp$lprop.1
    # par.vec[15] <- 1e-4 + 20 * par.vec[15]        # P$jmp$lprop.2
    par.vec[15] <- 1e-4 + 20 * par.vec[15]        # P$jmp$lprop.3
    par.vec[16] <- 1e-4 + 20 * par.vec[16]        # P$jmp$lvec
    par.vec[17] <- 1+1/(1e-2 + par.vec[17])       # P$jmp$muSc
    par.vec[18] <- -0.2 + 0.25 * par.vec[18]       # P$jmp$rhoc
    par.vec[19] <- 1e-4 + 0.12 * par.vec[19]       # P$jmp$sigmaYc
    #par.vec[20] <- 1e-2 + 2 * par.vec[20]         # Q$1$eta
    par.vec[20] <- 1e-2 + 16 * par.vec[20]        # Q$1$kpp
    par.vec[21] <- 1e-2 + 2 * par.vec[21]         # Q$2$eta
    par.vec[22] <- 1e-2 + 16 * par.vec[22]        # Q$2$kpp
    par.vec[23] <- 1e-2 + 16 * par.vec[23]         # Q$3$kpp
    par.vec[24] <- 1e-4 + 20 * par.vec[24]        # Q$jmp$lprop.1
    par.vec[25] <- 1+1/(1e-2 + par.vec[25])       # Q$jmp$muSc
    par.vec[26] <- -0.2 + 0.25 * par.vec[26]       # Q$jmp$rhoc
    par.vec[27] <- 1e-4 + 0.12 * par.vec[27]       # Q$jmp$sigmaYc
    
    for(kk in 28:(28+U-1)){
      par.vec[kk] <- 1e-4 + 1e-1 * par.vec[kk]
    }
    
    return(par.vec)
  }
  
  return(model.spec)  
}

#' @rdname modSetup
#' @details \code{spec_3FsepIntModel_portfolio_DJ} Specifies a three-factor model with leverage effect, exponential jumps in volatility and correlated double exponential jumps in the underlying, with same tail parameter on both sides. One of the factors only drives jump intensity, not the continuous volatility. There is a more flexible parametrisation of the equity risk premium parameters, with some Q parameters restricted. Finally, the first vol factor has pure-jump dynamics
#' @export

spec_3FsepIntModel_portfolio_DJ <- function(){
  
  N.factors <- 3
  model.spec <- model_makeDefaultParameterStructures(N.factors = N.factors, pq.equality = c("Q$jmp$lvec",paste0("Q$jmp$lprop.",c(2,3)),"Q$jmp$muYc",paste0("Q$",c(1,3),"$eta"),"Q$jmp$muStock"))
  
  model.spec$par.names <- model.spec$par.names[-which(grepl("sigma|mu", model.spec$par.names))]
  
  model.spec$par.names <- c(model.spec$par.names, c("P$jmp$gammaProp","Q$jmp$gammaProp","P$jmp$muVol","P$jmp$muStock2","Q$jmp$muVol","Q$jmp$muStock2"))
  model.spec$par.restr <- rbind(data.frame(par.name = "P$jmp$muStock", par.value = 0),model.spec$par.restr)
  
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.jmp.lvec",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$jmp$lvec", par.value = 0),model.spec$par.restr)
  
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.3.phi",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$3$phi", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.3.rho",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$3$rho", par.value = 0),model.spec$par.restr)
  
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.jmp.lprop.2",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$jmp$lprop.2", par.value = 0),model.spec$par.restr)
  
  model.spec$N.factors <- N.factors
  model.spec$jump.type <- 'oneSidedExponential_2'
  
  model.spec$par.names <- model.spec$par.names[order(model.spec$par.names)]
  model.spec$par.restr <- model.spec$par.restr[order(as.character(model.spec$par.restr$par.name)),]
  
  model.spec$scaleFoo <- function(par.vec){
    
    par.vec[1] <- -0.1 + 0.3*par.vec[1]           # P$1$erp0
    par.vec[2] <- 1e-2 + 16 * par.vec[2]          # P$1$kpp
    par.vec[3] <- 1e-4 + 4 * par.vec[3]           # P$1$lmb
    par.vec[4] <- 1e-4 + 0.3 * par.vec[4]         # P$1$phi
    par.vec[5] <- -1 + 2*par.vec[5]               # P$1$rho
    par.vec[6] <- 1e-2 + 16 * par.vec[6]          # P$2$kpp
    par.vec[7] <- 1e-4 + 4 * par.vec[7]         # P$2$lmb
    par.vec[8] <- 1e-4 + 0.3 * par.vec[8]         # P$2$phi
    par.vec[9] <- -1 + 2*par.vec[9]               # P$2$rho
    par.vec[10] <- 1e-2 + 16 * par.vec[10]        # P$3$kpp
    par.vec[11] <- 1e-4 + 4 * par.vec[11]         # P$3$lmb
    par.vec[12] <- 1e-4 + 7 * par.vec[12]        # P$jmp$gammaProp
    par.vec[13] <- 1e-4 + 20 * par.vec[13]        # P$jmp$lprop.1
    par.vec[14] <- 1e-4 + 20 * par.vec[14]        # P$jmp$lprop.3
    par.vec[15] <- 1e-4 + 0.3 * par.vec[15]       # P$jmp$muStock2
    par.vec[16] <- 1e-4 + 0.5 * par.vec[16]       # P$jmp$muVol
    par.vec[17] <- -0.2 + 0.25 * par.vec[17]       # P$jmp$rhoc
    par.vec[18] <- 1e-2 + 16 * par.vec[18]        # Q$1$kpp
    par.vec[19] <- 1e-2 + 2 * par.vec[19]         # Q$2$eta
    par.vec[20] <- 1e-2 + 16 * par.vec[20]        # Q$2$kpp
    par.vec[21] <- 1e-2 + 16 * par.vec[21]         # Q$3$kpp
    par.vec[22] <- 1e-4 + 7 * par.vec[22]        # Q$jmp$gammaProp
    par.vec[23] <- 1e-4 + 20 * par.vec[23]        # Q$jmp$lprop.1
    par.vec[24] <- 1e-4 + 0.3 * par.vec[24]       # Q$jmp$muStock2
    par.vec[25] <- 1e-4 + 0.5 * par.vec[25]       # Q$jmp$muVol
    par.vec[26] <- -0.2 + 0.25 * par.vec[26]       # Q$jmp$rhoc
    
    par.vec[27] <- 1      # noise covariance matrix scaling
    
    return(par.vec)
  }
  
  return(model.spec)  
}

#' @rdname modSetup
#' @details \code{spec_3FsepIntModel_erp} Specifies a three-factor model with leverage effect, exponential jumps in volatility and correlated double exponential jumps in the underlying, with same tail parameter on both sides. One of the factors only drives jump intensity, not the continuous volatility. There is a more flexible parametrisation of the equity risk premium parameters, with some Q parameters restricted. Finally, the first vol factor has pure-jump dynamics
#' @export

spec_3FsepIntModel_erp_portfolio <- function(){
  
  N.factors <- 3
  model.spec <- model_makeDefaultParameterStructures(N.factors = N.factors, pq.equality = c("Q$jmp$lvec",paste0("Q$jmp$lprop.",c(2,3)),"Q$jmp$muYc",paste0("Q$",c(1,3),"$eta")))
  model.spec$par.names <- c(model.spec$par.names, "P$2$erp")
  model.spec$par.restr <- model.spec$par.restr[-which(grepl("P.2.erp",model.spec$par.restr$par.name)),]
  model.spec$par.names <- c(model.spec$par.names, "P$1$erp")
  model.spec$par.restr <- model.spec$par.restr[-which(grepl("P.1.erp",model.spec$par.restr$par.name)),]
  model.spec$par.names <- c(model.spec$par.names, "P$3$erp")
  model.spec$par.restr <- model.spec$par.restr[-which(grepl("P.3.erp",model.spec$par.restr$par.name)),]
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.jmp.muYc",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$jmp$muYc", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.3.phi",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$3$phi", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.3.rho",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$3$rho", par.value = 0),model.spec$par.restr)
  # model.spec$par.names <- model.spec$par.names[-which(grepl("P.2.rho",model.spec$par.names))]
  # model.spec$par.restr <- rbind(data.frame(par.name = "P$2$rho", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.jmp.lprop.2",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$jmp$lprop.2", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.1.lmb",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$1$lmb", par.value = 1e-2),model.spec$par.restr)
  model.spec$N.factors <- N.factors
  model.spec$jump.type <- 'kouExpJumpTransform'
  
  model.spec$par.names <- model.spec$par.names[order(model.spec$par.names)]
  model.spec$par.restr <- model.spec$par.restr[order(as.character(model.spec$par.restr$par.name)),]
  
  model.spec$scaleFoo <- function(par.vec){
    
    par.vec[1] <- -0.3 + 0.6*par.vec[1]           # P$1$erp
    par.vec[2] <- -0.3 + 0.6*par.vec[2]           # P$1$erp0
    par.vec[3] <- 1e-2 + 16 * par.vec[3]          # P$1$kpp
    # par.vec[4] <- 1e-4 + 4 * par.vec[4]           # P$1$lmb
    par.vec[4] <- 1e-4 + 0.3 * par.vec[4]         # P$1$phi
    par.vec[5] <- -1 + 2*par.vec[5]               # P$1$rho
    par.vec[6] <- -0.3 + 0.6 * par.vec[6]           # P$2$erp
    par.vec[7] <- 1e-2 + 16 * par.vec[7]          # P$2$kpp
    par.vec[8] <- 1e-4 + 4 * par.vec[8]         # P$2$lmb
    par.vec[9] <- 1e-4 + 0.3 * par.vec[9]         # P$2$phi
    par.vec[10] <- -1 + 2*par.vec[10]               # P$2$rho
    par.vec[11] <- -0.3 + 0.6*par.vec[11]           # P$3$erp
    par.vec[12] <- 1e-2 + 16 * par.vec[12]        # P$3$kpp
    par.vec[13] <- 1e-4 + 4 * par.vec[13]         # P$3$lmb
    par.vec[14] <- 1e-4 + 20 * par.vec[14]        # P$jmp$lprop.1
    # par.vec[15] <- 1e-4 + 20 * par.vec[15]        # P$jmp$lprop.2
    par.vec[15] <- 1e-4 + 20 * par.vec[15]        # P$jmp$lprop.3
    par.vec[16] <- 1e-4 + 20 * par.vec[16]        # P$jmp$lvec
    par.vec[17] <- 1+1/(1e-2 + par.vec[17])       # P$jmp$muSc
    par.vec[18] <- -0.2 + 0.25 * par.vec[18]       # P$jmp$rhoc
    par.vec[19] <- 1e-4 + 0.12 * par.vec[19]       # P$jmp$sigmaYc
    #par.vec[20] <- 1e-2 + 2 * par.vec[20]         # Q$1$eta
    par.vec[20] <- 1e-2 + 16 * par.vec[20]        # Q$1$kpp
    par.vec[21] <- 1e-2 + 2 * par.vec[21]         # Q$2$eta
    par.vec[22] <- 1e-2 + 16 * par.vec[22]        # Q$2$kpp
    par.vec[23] <- 1e-2 + 16 * par.vec[23]         # Q$3$kpp
    par.vec[24] <- 1e-4 + 20 * par.vec[24]        # Q$jmp$lprop.1
    par.vec[25] <- 1+1/(1e-2 + par.vec[25])       # Q$jmp$muSc
    par.vec[26] <- -0.2 + 0.25 * par.vec[26]       # Q$jmp$rhoc
    par.vec[27] <- 1e-4 + 0.12 * par.vec[27]       # Q$jmp$sigmaYc
    
    par.vec[28] <- 5e2 + 5e3 * par.vec[28]      # noise covariance matrix scaling
    
    return(par.vec)
  }
  
  return(model.spec)  
}


#' @rdname modSetup
#' @details \code{spec_3FsepIntModel_smpl} specifies a three-factor model with minimalist erp form and two jump transforms: one for asset/vol, second for vol/intensity
#' @export

spec_3FsepIntModel_smpl <- function(){
  
  N.factors <- 3
  model.spec <- model_makeDefaultParameterStructures(N.factors = N.factors, pq.equality = c("Q$jmp$lvec",paste0("Q$jmp$lprop.",c(1:3)),paste0("Q$",c(1,2,3),"$eta"),"Q$jmp$muVol","Q$jmp$muInt","Q$jmp$muStock"))
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.2.phi",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$2$phi", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.2.rho",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$2$rho", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.3.rho",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$3$rho", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.3.lmb",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$3$lmb", par.value = 1e-4),model.spec$par.restr)
  
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.jmp.lprop.1",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$jmp$lprop.1", par.value = 0),model.spec$par.restr)
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.jmp.lprop.3",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$jmp$lprop.3", par.value = 0),model.spec$par.restr)
  
  model.spec$par.names <- model.spec$par.names[-which(grepl("P.jmp.lvec",model.spec$par.names))]
  model.spec$par.restr <- rbind(data.frame(par.name = "P$jmp$lvec", par.value = 0),model.spec$par.restr)
  
  model.spec$par.names <- c(model.spec$par.names, c("P$jmp$gammaProp","P$jmp$muStock","P$jmp$muVol","P$jmp$muInt","P$jmp$muVol2","Q$jmp$muVol2","Q$jmp$gammaProp"))
  
  model.spec$par.names <- model.spec$par.names[-which(grepl("*muSc|*muYc|*sigmaYc",model.spec$par.names))]
  
  model.spec$N.factors <- N.factors
  model.spec$jump.type <- 'oneSidedExponential'
  
  model.spec$par.names <- model.spec$par.names[order(model.spec$par.names)]
  model.spec$par.restr <- model.spec$par.restr[order(as.character(model.spec$par.restr$par.name)),]
  
  model.spec$scaleFoo <- function(par.vec){
    
    par.vec[1] <- 0.0 + 0.1*par.vec[1]           # P$1$erp0
    par.vec[2] <- 1e-2 + 16 * par.vec[2]          # P$1$kpp
    par.vec[3] <- 1e-4 + 4 * par.vec[3]           # P$1$lmb
    par.vec[4] <- 1e-4 + 0.3 * par.vec[4]         # P$1$phi
    par.vec[5] <- -1 + 2*par.vec[5]               # P$1$rho
    par.vec[6] <- 1e-2 + 16 * par.vec[6]          # P$2$kpp
    par.vec[7] <- 1e-4 + 4 * par.vec[7]         # P$2$lmb
    # par.vec[9] <- 1e-4 + 0.3 * par.vec[9]         # P$2$phi
    # par.vec[10] <- -1 + 2*par.vec[10]               # P$2$rho
    # par.vec[11] <- -0.3 + 0.6*par.vec[11]           # P$3$erp
    par.vec[8] <- 1e-2 + 16 * par.vec[8]        # P$3$kpp
    par.vec[9] <- 1e-4 + 0.3 * par.vec[9]         # P$3$phi
    # par.vec[13] <- 1e-4 + 4 * par.vec[13]         # P$3$lmb
    par.vec[10] <- 1e-4 + 5 * par.vec[10]        # P$jmp$gammaProp
    par.vec[11] <- 1e-4 + 20 * par.vec[11]        # P$jmp$lprop.2
    # par.vec[15] <- 1e-4 + 20 * par.vec[15]        # P$jmp$lprop.3
    par.vec[12] <- 1e-4 + 0.5 * par.vec[12]        # P$jmp$muInt
    par.vec[13] <- 1e-4 + 0.2 * par.vec[13]        # P$jmp$muStock
    par.vec[14] <- 1e-4 + 0.5 * par.vec[14]        # P$jmp$muVol
    par.vec[15] <- 1e-4 + 0.5 * par.vec[15]        # P$jmp$muVol2
    # par.vec[17] <- 1+1/(1e-2 + par.vec[17])       # P$jmp$muSc
    par.vec[16] <- -0.2 + 0.25 * par.vec[16]       # P$jmp$rhoc
    #par.vec[20] <- 1e-2 + 2 * par.vec[20]         # Q$1$eta
    par.vec[17] <- 1e-2 + 16 * par.vec[17]        # Q$1$kpp
    par.vec[18] <- 1e-2 + 16 * par.vec[18]        # Q$2$kpp
    par.vec[19] <- 1e-2 + 16 * par.vec[19]         # Q$3$kpp
    par.vec[20] <- 1e-4 + 5 * par.vec[20]        # Q$jmp$gammaProp
    par.vec[21] <- 1e-4 + 0.5 * par.vec[21]        # Q$jmp$muVol2
    par.vec[22] <- -0.2 + 0.25 * par.vec[22]       # Q$jmp$rhoc
    
    par.vec <- c(par.vec,1) # add that for handling the noisematrix in stateHandler
    return(par.vec)
  }
  
  return(model.spec)  
}

# @rdname modSetup
# @details \code{specData_DS_1M_6M_0115_cor_p0.0_p0.5} returns loads the divergence sample from 2001 to 2015, p= 0.5, maturities 1/12 and 1/2 years. Includes SP500 returns from OMTR and the bootstrapped error correlation matrices. Excludes quarticity
# @export

# specData_DS_1M_6M_0115_cor_p0.0_p0.5 <- function(path.to.data){
#   
#   load(paste0(path.to.data,"divergence-prices-omtr-asymptotic-iv-for-estimation.RData"))
#   load(paste0(path.to.data,"SP500_daily.RData"))
#   
#   pVec <- c(0,0.5)
#   tVec <- c(1/12,1/2)
#   
#   mkt <- data.frame(p=1,q=0,r=0,t=tVec)
#   
#   require(dplyr)
#   
#   date.subset <- which(unique(db$day) >= as.Date("2001-01-01"))
#   ret.start.date <- unique(db$day)[date.subset[1]-1]
#   
#   data.spec <- expand.grid(p= unique(db$u), t= unique(db$t), type = c('div','skew','quart'))
#   obs.spec <- expand.grid(p= pVec, t= tVec, type = c('div','skew'))
#   
#   db <- db %>% filter(day >= as.Date("2001-01-01")) %>% filter(u %in% pVec) %>% filter(t %in% tVec)
#   db <- db %>% mutate(div = div/t)
#   
#   obs.data <- db %>% select(day) %>% distinct
#   
#   for(tp in c("div","skew_scaled")){
#     for(tt in tVec){
#       for(pp in pVec){
#         tmp <- db %>% filter(t == tt, u == pp) %>% select_("day",tp)
#         colnames(tmp) <- c("day",paste0(tp,"_",sprintf("%1.3f",tt),"_",sprintf("%1.2f",pp)))
#         obs.data <- inner_join(obs.data,tmp,by="day")
#       }
#     }
#   }
#   # obs.data <- db %>% filter(t == 1/12) %>% select(day,div,skew_scaled)
#   # obs.data <- cbind(obs.data, db %>% filter(t == 1/2) %>% select(div,skew_scaled))
#   # colnames(obs.data) <- c("day","div_0.083_0.5","skew_0.083_0.5","div_0.5_0.5","skew_0.5_0.5")
#   # obs.data <- obs.data %>% select(day, div_0.083_0.5, div_0.5_0.5, skew_0.083_0.5, skew_0.5_0.5)
#   
#   SP500_daily <- SP500_daily %>% filter(date >= ret.start.date, date <= max(obs.data$day))
#   date.ranges <- cbind.data.frame(start= c(ret.start.date,head(obs.data$day,-1)), end = obs.data$day)
#   SP500_weekly <- apply(X = date.ranges, MARGIN = 1, FUN = function(x){
#     loc.sp <- SP500_daily %>% filter(date > x["start"], date <= x["end"])
#     ret.df <- data.frame(date = x["end"], ret = prod(1 + loc.sp$ret)-1)
#     return(ret.df)
#   })
#   SP500_weekly <- bind_rows(SP500_weekly)
#   SP500_weekly <- SP500_weekly %>% mutate(date = as.Date(date))  
#   obs.data <- SP500_weekly %>% inner_join(obs.data, by = c("date"="day"))
#   
#   data.spec <- cbind(data.spec,index = 1:nrow(data.spec))
#   obs.spec <- data.spec %>% inner_join(obs.spec)
#   
#   noise.cov.cube <- db.vcov[obs.spec$index,obs.spec$index,date.subset]
#   noise.cov.cube[which(obs.spec$type=='div'),which(obs.spec$type=='div'),] <- noise.cov.cube[which(obs.spec$type=='div'),which(obs.spec$type=='div'),] * array(outer(1/rep(tVec,each=length(pVec)), 1/rep(tVec,each=length(pVec))), dim = c(length(tVec)*length(pVec), length(tVec)*length(pVec), dim(noise.cov.cube)[3]))
#   dt <- 7/365
#   
#   noise.cov.cube <- apply(noise.cov.cube,3,cov2cor)
#   noise.cov.cube <- array(as.numeric(noise.cov.cube), dim = c(ncol(obs.data)-2,ncol(obs.data)-2,nrow(obs.data)))
#   
#   spec.mat <- obs.spec %>% select(-index)
#   return(list(dt = dt, noise.cov.cube = noise.cov.cube, spec.mat = spec.mat, obs.data = as.matrix(obs.data[,-1]), dates = obs.data[,1], mkt = mkt))
# }

#' @rdname modSetup
#' @details \code{specData_DSQ_3Mat} returns loads the divergence sample from 2001 to 2015, p= 0.5, maturities 1/12 and 1/2 years. Includes SP500 returns from OMTR and the bootstrapped error correlation matrices. Excludes quarticity
#' @export
specData_PF_DSQ_3Mat <- function(path.to.data){
  
  load(paste0(path.to.data,"spx-pfolio-drgs-noise.RData"))
  load(paste0(path.to.data,"SP500_daily.RData"))
  
  date.subset <- which(unique(div.db$date) >= as.Date("2001-01-01"))
  ret.start.date <- unique(div.db$date)[date.subset[1]-1]
  
  mkt.list <- mkt.list[date.subset]
  wts.list <- wts.list[date.subset]
  strike.list <- strike.list[date.subset]
  
  div.db <- div.db %>% dplyr::filter(date >= as.Date("2001-01-01"))
  div.db <- div.db %>% dplyr::select(div_mid, skew_mid, quart_mid)
  
  obs.data <- matrix(as.numeric(t(as.matrix(div.db[,-1]))), nrow = nrow(div.db)/3, ncol = (ncol(div.db)-1)*3, byrow = T)
  obs.data <- data.frame(day = unique(div.db$date), obs.data)
  
  SP500_daily <- SP500_daily %>% dplyr::filter(date >= ret.start.date, date <= max(obs.data$day))
  date.ranges <- cbind.data.frame(start= c(ret.start.date,head(obs.data$day,-1)), end = obs.data$day)
  SP500_weekly <- apply(X = date.ranges, MARGIN = 1, FUN = function(x){
    loc.sp <- SP500_daily %>% dplyr::filter(date > x["start"], date <= x["end"])
    ret.df <- data.frame(date = x["end"], ret = prod(1 + loc.sp$ret)-1)
    return(ret.df)
  })
  SP500_weekly <- bind_rows(SP500_weekly)
  SP500_weekly <- SP500_weekly %>% dplyr::mutate(date = as.Date(date))  
  obs.data <- SP500_weekly %>% dplyr::inner_join(obs.data, by = c("date"="day"))
  
  noise.cov.cube <- noise.cov.cube[,,(dim(noise.cov.cube)[3] - nrow(obs.data) + 1):dim(noise.cov.cube)[3], drop = F]
  noise.cov.cube <- noise.cov.cube/16.0
  
  dt <- 7/365
  
  return(list(dt = dt, noise.cov.cube = noise.cov.cube, obs.data = as.matrix(obs.data[,-1]), dates = obs.data[,1], mkt.list = mkt.list, wts.list = wts.list, strikeMat.list = strike.list))
}


#' @rdname modSetup
#' @details \code{specData_PF_D_3Mat} returns loads the divergence sample from 2001 to 2015, p= 0.5, maturities 1/12 and 1/2 years. Only divergence portfolios are taken. Includes SP500 returns from OMTR and the bootstrapped error correlation matrices. Excludes quarticity
#' @export
specData_PF_D_3Mat <- function(path.to.data){
  
  res <- specData_PF_DSQ_3Mat(path.to.data)
  
  res$wts.list <- lapply(res$wts.list, function(x) return(x[,,1,drop=F]))
  res$obs.data <- res$obs.data[,c(1,2,5,8)]
  res$noise.cov.cube <- res$noise.cov.cube[c(1,4,7),c(1,4,7),]
  
  return(res)
}


