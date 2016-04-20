#### ---- DYNAMICS ----
library(divergenceModelR)
load("D:/git/hfAffine/RpackageDev/affineOption/data/heston.params.RData")
parListHestonLev$Q$`1`$phi <- parListHestonLev$P$`1`$phi <- 0.3
parListHestonLev$Q$`1`$lmb <- parListHestonLev$P$`1`$lmb <- 1.2

aa <- modelDynamics(params.P = parListHestonLev$P, params.Q = parListHestonLev$Q, dT = 5/252, N.factors = 1, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.points = 2, mod.type = 'standard')

testMean <- meanVecFun(meanListS = aa$mean.vec, currVol = matrix(1,1,1))
testCov <- covMatFun(covListS = aa$cov.list, covListDim = c(2,2) ,  currVol = matrix(1e-9,1,1))
testCov <- covMatFun(covListS = aa$cov.array[c(1,2,4)], covListDim = c(2,2) ,  currVol = matrix(1,1,1))

calc.par <- list(mean.vec = aa$mean.vec[-1], cov.array = aa$cov.array[2,2])

meancov <- affineTransitionStateHandler(stateMat = rbind(c(0.5,1,2),c(0,0,0)), modelParameters = calc.par,1)

#### ---- OBSERVATION ----
mkt.spec <- data.frame(p=1,q=0,r=0,t=c(1/12,0.5))
sol.check.derivs <- Re(odeExtSolveWrap(u = rbind(cbind(0.5,0),cbind(0,0)), params.P = parListHestonLev$P, params.Q = parListHestonLev$Q, mkt = mkt.spec, jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, mod.type = 'standard', rtol = 1e-12, atol = 1e-22))

obsParameters <- list(stockParams = list(mean.vec = aa$mean.vec, cov.array = aa$cov.array[lower.tri(diag(2) ,diag = T)]), cfCoeffs = sol.check.derivs, tVec = mkt.spec$t, pVec = c(0.5,0), cVec = rep(0.8,3), bVec = rep(0.05,3), divNoiseCube = array(1,dim=c(12,12,1)))

stobs <- affineObservationStateHandler(stateMat = rbind(c(0.98,1,1.1),c(1,1,1)), modelParameters = obsParameters, iterCount = 0)

#### ---- TEST WHOLE FILTER ----
load("D:/git/hfAffine/RpackageDev/affineOption/data/heston.params.RData")
parListHestonLev$Q$`1`$phi <- parListHestonLev$P$`1`$phi <- 0.3
parListHestonLev$Q$`1`$lmb <- parListHestonLev$P$`1`$lmb <- 1.2

mod.dyn <- modelDynamics(params.P = parListHestonLev$P, params.Q = parListHestonLev$Q, dT = 5/252, N.factors = 1, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.points = 2, mod.type = 'standard')


# generate test data
set.seed(123555)
initState <- matrix(c(1,1), ncol = 1)
initVol <- matrix(0,2,2)
initVol[1,1] <- covMatFun(covListS = aa$cov.array[c(1,2,4)], covListDim = c(2,2) ,  currVol = initState[1,1,drop=F])[2,2]

Ndays <- 4000

paths <- affineSimulate(paramsList = parListHestonLev, N.factors = 1, t.days = Ndays, t.freq = 1, freq.subdiv = 78)
stock.ret <- paths$S.array[seq(1,Ndays,by=5),"F"]
stock.ret <- diff(stock.ret)/head(stock.ret,-1)

vols <- tail(paths$V.array[seq(1,Ndays,by=5),"v1"],-1)
dvrg <- divergenceSwapRate(p = c(0,1/2), params.Q = parListHestonLev$Q, t.vec = mkt.spec$t, vol.mat = matrix(vols), mod.type = 'standard')
skew <- skewSwapRate(p = c(0,1/2), params.Q = parListHestonLev$Q, t.vec = mkt.spec$t, vol.mat = matrix(vols), jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), mod.type = 'standard')
quart <- quartSwapRate(p = c(0,1/2), params.Q = parListHestonLev$Q, t.vec = mkt.spec$t, vol.mat = matrix(vols), jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), mod.type = 'standard')

# obsData <- cbind(stock.ret,t(dvrg[1,,]),t(skew[1,,])/t(dvrg[1,,]^1.5),t(quart[1,,])/t(dvrg[1,,]^2))
if(dim(dvrg)[2]==1 & dim(dvrg)[1]==1){
  obsData <- cbind(stock.ret,drop(dvrg),drop(skew)/drop(dvrg^1.5),drop(quart)/drop(dvrg^2))
} else if(dim(dvrg)[2]==1 & dim(dvrg)[1]>1) {
  obsData <- cbind(stock.ret,t(drop(dvrg)),t(drop(skew)/drop(dvrg^1.5)),t(drop(quart)/drop(dvrg^2)))
} else {
  obsData <- cbind(stock.ret,matrix(dvrg,nrow=dim(dvrg)[3],ncol = sum(dim(dvrg)[1:2])),matrix(drop(skew)/drop(dvrg^1.5),nrow=dim(dvrg)[3],ncol = sum(dim(dvrg)[1:2])),matrix(drop(quart)/drop(dvrg^2),nrow=dim(dvrg)[3],ncol = sum(dim(dvrg)[1:2])))
}

obsDataTrue <- obsData
for(kk in 1:length(vols)){
  err <- t(chol(divModelObsNoiseMat(corrs = rep(0.9,3), bpars = c(0.005,0.025,0.055), spotVar = vols[kk], matVec = mkt.spec$t, U = 2))) %*% rnorm(ncol(obsData)-1)
  obsData[kk,-1] <- obsData[kk,-1] + err
}

ode.solutions <- Re(odeExtSolveWrap(u = rbind(cbind(0,0),cbind(0.5,0)), params.P = parListHestonLev$P, params.Q = parListHestonLev$Q, mkt = mkt.spec, jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, mod.type = 'standard', rtol = 1e-12, atol = 1e-22))

obsParameters <- list(stockParams = list(mean.vec = mod.dyn$mean.vec, cov.array = mod.dyn$cov.array[lower.tri(diag(2) ,diag = T)]), cfCoeffs = ode.solutions, tVec = mkt.spec$t, pVec = c(0,0.5), cVec = rep(0.9,3), bVec = c(0.005,0.025,0.055))

transParameters <- list(mean.vec = mod.dyn$mean.vec[-1], cov.array = mod.dyn$cov.array[2,2])

#### ---- CALL FILTER ----
# sink(file = 'ftestlog.txt', append = F, type = 'output')
# t.loc <- Sys.time()
# for(kk in 1:500){
#   # print(paste("iteration",kk,"started"))
#   ftest <- testAffineFilter(testDataMat = obsData, testInitState = initState, testInitProcCov = initVol, testModelParams = list(transition = transParameters, observation = obsParameters))
#   ftestSqrt <- testSqrtAffineFilter(testDataMat = obsData, testInitState = initState, testInitProcCov = initVol, testModelParams = list(transition = transParameters, observation = obsParameters))
# }
# t.loc <- difftime(Sys.time(),t.loc,units = 'secs')
# sink(NULL)
# ftestsame <- testAffineFilter(testDataMat = obsData, testInitState = initState, testInitProcCov = initVol, testModelParams = list(transition = transParameters, observation = obsParameters))

#### ---- CALL MOD LIK ----
data.structure <- list(obs.data = obsData, spec.mat = expand.grid(t = mkt.spec$t, p = c(0,0.5), type = c("div","skew","quart")))

model.spec <- list(params.P = parListHestonLev$P, params.Q = parListHestonLev$Q, jump.type = 'expNormJumpTransform', dt = 5/252, N.factors  = 1, error = list(cVec=rep(0.9,3), bVec = c(0.005,0.025,0.055)), mkt = mkt.spec)
# = 1, error = list(cVec=rep(0.9,3), bVec = c(0.005,0.015,0.015)), mkt = mkt.spec)

# lik.test <- modelLikelihood(data.structure = data.structure, model.spec = model.spec, for.estimation = TRUE, filterFoo = divergenceModelR:::DSQ_sqrtFilter)

# meancov <- affineTransitionStateHandler(stateMat = rbind(c(0.5,1,2),c(0,0,0)), modelParameters = calc.par)

lik.test.notest <- modelLikelihood(data.structure = data.structure, model.spec = model.spec, for.estimation = F, filterFoo = divergenceModelR:::DSQ_sqrtFilter)


### callers tests ----
load("D:/git/hfAffine/RpackageDev/affineOption/data/heston.params.RData")
parListHestonLev$Q$`1`$phi <- parListHestonLev$P$`1`$phi <- 0.3
parListHestonLev$Q$`1`$lmb <- parListHestonLev$P$`1`$lmb <- 1.2
parListHestonLev$P$jmp$lvec <- parListHestonLev$Q$jmp$lvec <- 1
parListHestonLev$P$jmp$lprop <- parListHestonLev$Q$jmp$lprop <- 1
parListHestonLev$P$jmp$muYc <- -0.03
parListHestonLev$Q$jmp$muYc <- -0.05
parListHestonLev$P$jmp$sigmaYc <- 0.05
parListHestonLev$Q$jmp$sigmaYc <- 0.05
parListHestonLev$P$jmp$muSc <- 1/0.15
parListHestonLev$Q$jmp$muSc <- 1/0.20

nm.list <- spec_1FtoyModel()$par.names
par.unlist <- unlist(parListHestonLev)
names(par.unlist) <- gsub(pattern = "\\.",replacement = "$", names(par.unlist))
names(par.unlist) <- gsub(pattern = "lprop",replacement = "lprop.1", names(par.unlist))
par.vec <- par.unlist[nm.list]
par.vec.2 <- par.vec
par.vec.2["P$1$erp"] <- 0.5
modLik <- model_wrapLikelihood(data.structure = specData_DSQ_1M_6M_0115("D:/git/lugtex/swapRateModel/optData/"), model.spec = spec_1FtoyModel(), for.estimation = T, filterFoo = DSQ_sqrtFilter, N.points = 3)

system.time(
  eval.test <- modLik(par.vec = par.vec)
)

library(parallel)
cl <- makeCluster(2)
clusterEvalQ(cl, library(divergenceModelR))
par.list <- list(par.vec, par.vec.2)

system.time(
  par.eval.test <- parLapply(cl = cl, X = par.list, fun = modLik)
)

library(nloptr)
