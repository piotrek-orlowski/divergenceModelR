#### ---- DYNAMICS ----
library(divergenceModelR)
load("D:/git/hfAffine/RpackageDev/affineOption/data/heston.params.RData")
parListHestonLev$Q$`1`$phi <- parListHestonLev$P$`1`$phi <- 0.3
parListHestonLev$Q$`1`$lmb <- parListHestonLev$P$`1`$lmb <- 0.8

aa <- modelDynamics(params.P = parListHestonLev$P, params.Q = parListHestonLev$Q, dT = 5/252, N.factors = 1, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.points = 2, mod.type = 'standard')

testMean <- meanVecFun(meanListS = aa$mean.vec, currVol = matrix(1,1,1))
testCov <- covMatFun(covListS = aa$cov.list, covListDim = c(2,2) ,  currVol = matrix(1,1,1))
testCov <- covMatFun(covListS = aa$cov.array[c(1,2,4)], covListDim = c(2,2) ,  currVol = matrix(1,1,1))

calc.par <- list(mean.vec = aa$mean.vec[-1], cov.array = aa$cov.array[2,2])

meancov <- affineTransitionStateHandler(stateMat = rbind(c(0.5,1,2),c(0,0,0)), modelParameters = calc.par)

#### ---- OBSERVATION ----
mkt.spec <- data.frame(p=1,q=0,r=0,t=c(0.5,1))

sol.check.derivs <- Re(odeExtSolveWrap(u = cbind(0.5,0), params.P = parListHestonLev$P, params.Q = parListHestonLev$Q, mkt = mkt.spec, jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, mod.type = 'standard', rtol = 1e-12, atol = 1e-22))

obsParameters <- list(stockParams = list(mean.vec = aa$mean.vec, cov.array = aa$cov.array[lower.tri(diag(2) ,diag = T)]), cfCoeffs = sol.check.derivs, tVec = mkt.spec$t, pVec = 0.5, cVec = rep(0.8,3), bVec = rep(0.05,3))

stobs <- affineObservationStateHandler(stateMat = rbind(c(0.98,1,1.1),c(1,1,1)), modelParameters = obsParameters)

#### ---- TEST WHOLE FILTER ----
load("D:/git/hfAffine/RpackageDev/affineOption/data/heston.params.RData")
parListHestonLev$Q$`1`$phi <- parListHestonLev$P$`1`$phi <- 0.3
parListHestonLev$Q$`1`$lmb <- parListHestonLev$P$`1`$lmb <- 0.8

mod.dyn <- modelDynamics(params.P = parListHestonLev$P, params.Q = parListHestonLev$Q, dT = 5/252, N.factors = 1, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.points = 2, mod.type = 'standard')


# generate test data
set.seed(123555)
initState <- matrix(c(1,1), ncol = 1)
initVol <- matrix(0,2,2)
initVol[1,1] <- covMatFun(covListS = aa$cov.array[c(1,2,4)], covListDim = c(2,2) ,  currVol = initState[1,1,drop=F])[2,2]

Ndays <- 400

paths <- affineSimulate(paramsList = parListHestonLev, N.factors = 1, t.days = Ndays, t.freq = 1, freq.subdiv = 78)
stock.ret <- paths$S.array[seq(1,Ndays,by=5),"F"]
stock.ret <- diff(stock.ret)/head(stock.ret,-1)

vols <- tail(paths$V.array[seq(1,Ndays,by=5),"v1"],-1)
dvrg <- divergenceSwapRate(p = 1/2, params.Q = parListHestonLev$Q, t.vec = mkt.spec$t, vol.mat = matrix(vols), mod.type = 'standard')
skew <- skewSwapRate(p = 1/2, params.Q = parListHestonLev$Q, t.vec = mkt.spec$t, vol.mat = matrix(vols), jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), mod.type = 'standard')
quart <- quartSwapRate(p = 1/2, params.Q = parListHestonLev$Q, t.vec = mkt.spec$t, vol.mat = matrix(vols), jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), mod.type = 'standard')

obsData <- cbind(stock.ret,t(dvrg[1,,]),t(skew[1,,])/t(dvrg[1,,]^1.5),t(quart[1,,])/t(dvrg[1,,]^2))

obsDataTrue <- obsData
for(kk in 1:length(vols)){
  err <- t(chol(divModelObsNoiseMat(corrs = rep(0.0,3), bpars = c(0.005,0.05,0.075), spotVar = vols[kk], matVec = mkt.spec$t, U = 1))) %*% rnorm(6)
  obsData[kk,-1] <- obsData[kk,-1] + err
}

ode.solutions <- Re(odeExtSolveWrap(u = cbind(0.5,0), params.P = parListHestonLev$P, params.Q = parListHestonLev$Q, mkt = mkt.spec, jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, mod.type = 'standard', rtol = 1e-12, atol = 1e-22))

obsParameters <- list(stockParams = list(mean.vec = mod.dyn$mean.vec, cov.array = mod.dyn$cov.array[lower.tri(diag(2) ,diag = T)]), cfCoeffs = ode.solutions, tVec = mkt.spec$t, pVec = 0.5, cVec = rep(0.0,3), bVec = c(0.005,0.05,0.075))

transParameters <- list(mean.vec = mod.dyn$mean.vec[-1], cov.array = mod.dyn$cov.array[2,2])

#### ---- CALL FILTER ----
sink(file = 'ftestlog.txt', append = F, type = 'output')
for(kk in 1:1000){
  print(paste("iteration",kk,"started"))
  ftest <- testAffineFilter(testDataMat = obsData, testInitState = initState, testInitProcCov = initVol, testModelParams = list(transition = transParameters, observation = obsParameters))
  ftestSqrt <- testSqrtAffineFilter(testDataMat = obsData, testInitState = initState, testInitProcCov = initVol, testModelParams = list(transition = transParameters, observation = obsParameters))
}
sink(NULL)
ftestsame <- testAffineFilter(testDataMat = obsData, testInitState = initState, testInitProcCov = initVol, testModelParams = list(transition = transParameters, observation = obsParameters))

