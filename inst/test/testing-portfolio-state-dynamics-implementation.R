#### ---- DYNAMICS ----
library(divergenceModelR)
load("D:/git/hfAffine/RpackageDev/affineOption/data/heston.params.RData")
parListHestonLev$Q$`1`$phi <- parListHestonLev$P$`1`$phi <- 0.3
parListHestonLev$Q$`1`$lmb <- parListHestonLev$P$`1`$lmb <- 1.2
parListHestonLev$P$`1`$erp <- 0
parListHestonLev$P$`1`$erp0 <- 0
parListHestonLev$P$`1`$rho <- -0.8
parListHestonLev$Q$`1`$rho <- -0.8

aa <- modelDynamics(params.P = parListHestonLev$P, params.Q = parListHestonLev$Q, dT = 5/252, N.factors = 1, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.points = 5, mod.type = 'standard')

testMean <- meanVecFun(meanListS = aa$mean.vec, currVol = matrix(1,1,1))
testCov <- covMatFun(covListS = aa$cov.list, covListDim = c(2,2) ,  currVol = matrix(1e-9,1,1))
testCov <- covMatFun(covListS = aa$cov.array[c(1,2,4)], covListDim = c(2,2) ,  currVol = matrix(1,1,1))

calc.par <- list(mean.vec = aa$mean.vec[-1], cov.array = aa$cov.array[2,2])

meancov <- affineTransitionStateHandler(stateMat = rbind(c(0.5,1,2),c(1,1,1)), modelParameters = calc.par,1)

#### ---- OBSERVATION ----
library(abind)
library(statmod)
mkt.spec <- data.frame(t=c(1/12,0.5),r=0.0,q=0,p=1)

mkt.mat <- as.matrix(mkt.spec)

strikeMat <- array(seq(-0.5,0.5,by=0.05), dim = c(length(seq(-0.5,0.5,by=0.05)),2,3))
strikeMat <- aperm(strikeMat, c(2,1,3))

wtsCube <- abind(exp(-2*strikeMat[,,1])*strikeMat[,,1]*(strikeMat[,,1] <=0),exp(-2*strikeMat[,,1]), along = 3)

gl <- gauss.quad(n = 64, kind = "laguerre", alpha = 0)

cfCoeff <- jumpDiffusionODEs(u = cbind(1i*gl$nodes, matrix(0,nrow = length(gl$nodes), ncol = 1)), params = parListHestonLev$Q, mkt = mkt.spec, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, rtol = 1e-14, atol = 1e-28, N.factors = 1, mod.type = 'standard')

obsParameters <- list(stockParams = list(mean.vec = aa$mean.vec, cov.array = aa$cov.array[lower.tri(diag(2) ,diag = T)]), cfCoeffs = list(cfCoeff), strikeMats = list(strikeMat), mkts = list(mkt.mat), wts = list(wtsCube), quadWeights = gl$weights, quadNodes = gl$nodes, divNoiseCube = array(1,dim=c(4,4,1)), errSdParVec = rep(0,4))

cfVals <- affineCF(u = cbind(1i*gl$nodes, matrix(0,nrow = length(gl$nodes), ncol = 1)), params.Q = parListHestonLev$Q, params.P = NULL, t.vec = NULL, v.0 = matrix(c(1,0.5,1.5),nrow=3), jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.factors = 1, mod.type = 'standard', mkt = mkt.spec)

optpr <- glPricer_cpp(strikeMat = strikeMat, mkt = mkt.mat, glWts = gl$weights, glNodes = gl$nodes, cfVals = cfVals, Nfactors = 1)

stobs <- affineObservationStateHandler_optionPortfolios(stateMat = cbind(rbind(1,1),rbind(0.5,1),rbind(1.5,1)), modelParameters = obsParameters, iterCount = 0)

#### ---- TEST WHOLE FILTER ----
load("D:/git/hfAffine/RpackageDev/affineOption/data/heston.params.RData")
parListHestonLev$Q$`1`$phi <- parListHestonLev$P$`1`$phi <- 0.3
parListHestonLev$Q$`1`$lmb <- parListHestonLev$P$`1`$lmb <- 1.2
parListHestonLev$P$`1`$erp <- 0
parListHestonLev$P$`1`$erp0 <- 0.05
parListHestonLev$P$jmp$lvec <- parListHestonLev$Q$jmp$lvec <- 1
parListHestonLev$P$jmp$lprop <- parListHestonLev$Q$jmp$lprop <- 1
parListHestonLev$P$jmp$muYc <- -0.03
parListHestonLev$Q$jmp$muYc <- -0.05
parListHestonLev$P$jmp$sigmaYc <- 0.05
parListHestonLev$Q$jmp$sigmaYc <- 0.05
parListHestonLev$P$jmp$muSc <- 0.15
parListHestonLev$Q$jmp$muSc <- 0.20
parListHestonLev$P$jmp$rhoc <- -0.1
parListHestonLev$Q$`1`$rho <- parListHestonLev$P$`1`$rho <- -0.99

mod.dyn <- modelDynamics(params.P = parListHestonLev$P, params.Q = parListHestonLev$Q, dT = 5/252, N.factors = 1, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.points = 2, mod.type = 'standard')

# generate test data
set.seed(123555)
initState <- matrix(c(1,1), ncol = 1)
initVol <- matrix(0,2,2)
initVol[1,1] <- covMatFun(covListS = mod.dyn$cov.array[c(1,2,4)], covListDim = c(2,2) ,  currVol = initState[1,1,drop=F])[2,2]

Ndays <- 3000

paths <- affineSimulate(paramsList = parListHestonLev, N.factors = 1, t.days = Ndays, t.freq = 1, freq.subdiv = 78)
stock.ret <- paths$sim.arrays[[1]]$S.array[seq(1,Ndays,by=5),1]
stock.ret <- diff(stock.ret)/head(stock.ret,-1)

vols <- tail(paths$sim.arrays[[1]]$V.array[seq(1,Ndays,by=5),1],-1)

strikeMat <- array(seq(-0.5,0.5,by=0.05), dim = c(length(seq(-0.5,0.5,by=0.01)),2,length(vols)))
strikeMat <- aperm(strikeMat, c(2,1,3))

wtsCube <- abind(exp(-2*strikeMat[,,1])*strikeMat[,,1]*(strikeMat[,,1] <=0),exp(-2*strikeMat[,,1]), along = 3)

cfCoeffs <- strikeMat.list <- mkt.list <- wts.list <- vector(mode = "list", length = length(vols))

for(kk in 1:length(vols)){
  strikeMat.list[[kk]] <- strikeMat[,,kk,drop=FALSE]
  wts.list[[kk]] <- wtsCube
  mkt.list[[kk]] <- mkt.mat
}

cfCoeff <- jumpDiffusionODEs(u = cbind(1i*gl$nodes, matrix(0,nrow = length(gl$nodes), ncol = 1)), params = parListHestonLev$Q, mkt = mkt.spec, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, rtol = 1e-14, atol = 1e-28, N.factors = 1, mod.type = 'standard')

cfVals <- affineCF(u = cbind(1i*gl$nodes, matrix(0,nrow = length(gl$nodes), ncol = 1)), params.Q = parListHestonLev$Q, params.P = NULL, t.vec = NULL, v.0 = matrix(vols,nrow=length(vols)), jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.factors = 1, mod.type = 'standard', mkt = mkt.spec)

optpr <- glPricer_cpp(strikeMat = strikeMat, mkt = mkt.mat, glWts = gl$weights, glNodes = gl$nodes, cfVals = cfVals, Nfactors = 1)

optpf <- apply(X = optpr, MARGIN = 3, FUN = function(x){
  loc <- wtsCube * abind(x,x,along = 3)
  res <- as.numeric(apply(loc,c(1,3),sum))
  return(res)
})

optpf <- t(optpf)
optpf <- optpf[,c(1,3,2,4)]

obsData <- cbind(stock.ret, optpf)

noiseCube <- array(0,c(ncol(obsData)-1, ncol(obsData)-1, length(vols)))
for(kk in 1:length(vols)){
  noiseCube[,,kk] <- diag(c(0.0005,0.005,0.008,0.03)^2)
  cfCoeffs[[kk]] <- cfCoeff
}

obsDataTrue <- obsData
for(kk in 1:length(vols)){
  err <- diag(c(0.0005,0.005,0.008,0.03)) %*% rnorm(ncol(obsData)-1)
  obsData[kk,-1] <- obsData[kk,-1] + err
}

N.factors<-1

transition.parameters <- list(
  mean.vec = mod.dyn$mean.vec[-1], 
  cov.array = mod.dyn$cov.array[2,2]
)


# Observation equation setup
observation.parameters <- list(
  stockParams = list(mean.vec = mod.dyn$mean.vec, cov.array = mod.dyn$cov.array[lower.tri(diag(1+N.factors),diag = T)]), 
  cfCoeffs = cfCoeffs,
  divNoiseCube = noiseCube,
  strikeMats = strikeMat.list,
  mkts = mkt.list,
  wts = wts.list,
  errSdParVec = rep(1,4),
  quadWeights = gl$weights,
  quadNodes = gl$nodes
)

modelParams <- list(transition = transition.parameters, observation = observation.parameters)

testMean <- meanVecFun(meanListS = mod.dyn$mean.vec, currVol = matrix(1,1,1))
testCov <- covMatFun(covListS = mod.dyn$cov.list, covListDim = c(2,2) ,  currVol = matrix(1,1,1))

#### RUN ####
icov <- matrix(0,2,2)
icov[1,1] <- testCov[2,2]
system.time(test <- divergenceModelR:::portfolio_sqrtFilter(dataMat = obsData, initState = rbind(1,1), initProcCov = icov, modelParams = modelParams))

#### GET LIKELIHOOD ####

modLik <- model_wrapLikelihood_portfolio_extraNoise(data.structure = specData_PF_DSQ_3Mat("d:/git/LugTeX/swapRateModel/optData/"), model.spec = model.spec, for.estimation = T, filterFoo = divergenceModelR:::portfolio_sqrtFilter, N.points = 3, penalized = T, penalty = 1e8)
