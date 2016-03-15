###
#
# TEST FILTER WITH 2F MODEL
#
###

# ---- SETUP ----
library(divergenceModelR)

# ---- MODEL PARAMETERS ----
## Case #3 -- 2F no volatility jumps, Inspired by Andersen, Todorov and Fusari (2015)
par.A1 <- vector(mode = "list", length = 3)
names(par.A1) <- c("1","2","jmp")
par.A1[["1"]] <- list(eta = 1, lmb = 1, kpp = 1, rho = -0.9, phi = 0.2, erp = 0)
par.A1[["2"]] <- list(eta = 1, lmb = 2, kpp = 10, rho = -0.9, phi = 0.2, erp = 0)
par.A1[["jmp"]] <- list(lvec = 2, lprop = c(200,100), muYc = -0.0003, sigmaYc = 3e-3, muSc = 1/0.0003, rhoc = -0.2)

parList <- list(P = par.A1, Q = par.A1)
parList$Q$jmp[c("muYc","sigmaYc","muSc")] <- c(-4e-4,4e-3,1/0.00035)
# ---- SIMULATE DATA ----
set.seed(123555)
# initState <- matrix(c(1,1), ncol = 1)
# initVol <- matrix(0,2,2)
# initVol[1,1] <- covMatFun(covListS = aa$cov.array[c(1,2,4)], covListDim = c(2,2) ,  currVol = initState[1,1,drop=F])[2,2]

Ndays <- 4000
mkt.spec <- data.frame(p=1,r=0,q=0,t=c(1/12,1/2))

paths <- affineSimulate(paramsList = parList, N.factors = 2, t.days = Ndays, t.freq = 1, freq.subdiv = 78, jumpGeneratorPtr = getPointerToGenerator("kouExpJumpTransform"), jumpTransformPtr = getPointerToJumpTransform('kouExpJumpTransform')$TF)

stock.ret <- paths$S.array[seq(1,Ndays,by=5),"F"]
stock.ret <- diff(stock.ret)/head(stock.ret,-1)

vols <- tail(paths$V.array[seq(1,Ndays,by=5),c("v1","v2")],-1)

dvrg <- divergenceSwapRate(p = c(0.1,1/2), params.Q = parList$Q, t.vec = mkt.spec$t, vol.mat = vols, mod.type = 'standard', jumpTransform = getPointerToJumpTransform('kouExpJumpTransform'))

skew <- skewSwapRate(p = c(0.1,1/2), params.Q = parList$Q, t.vec = mkt.spec$t, vol.mat = vols, jumpTransform = getPointerToJumpTransform('kouExpJumpTransform'), mod.type = 'standard')

quart <- quartSwapRate(p = c(0.1,1/2), params.Q = parList$Q, t.vec = mkt.spec$t, vol.mat = vols, jumpTransform = getPointerToJumpTransform('kouExpJumpTransform'), mod.type = 'standard')

# obsData <- cbind(stock.ret,t(dvrg[1,,]),t(skew[1,,])/t(dvrg[1,,]^1.5),t(quart[1,,])/t(dvrg[1,,]^2))
obsData <- matrix(0,nrow = length(stock.ret), ncol = 1+3*prod(dim(dvrg)[1:2]))
obsData[,1] <- stock.ret
for(kk in 1:nrow(obsData)){
  obsData[kk,2:5] <- as.numeric(dvrg[,,kk])/rep(mkt.spec$t,each=2)
  obsData[kk,6:9] <- as.numeric(skew[,,kk])/as.numeric(dvrg[,,kk])^1.5
  obsData[kk,10:13] <- as.numeric(quart[,,kk])/as.numeric(dvrg[,,kk])^2
}

obsDataTrue <- obsData
for(kk in 1:nrow(vols)){
  err <- t(chol(divModelObsNoiseMat(corrs = rep(0.0,3), bpars = c(0.001,0.0020,0.0030), spotVar = vols[kk], matVec = mkt.spec$t, U = dim(dvrg)[1]))) %*% rnorm(ncol(obsData)-1)
  obsData[kk,-1] <- obsData[kk,-1] + err
}

# ---- MODEL RUN ----

ode.solutions <- Re(odeExtSolveWrap(u = matrix(c(0.1,0,0,0.5,0,0),nrow=2, ncol=3,byrow=T), params.Q = parList$Q, mkt = mkt.spec, jumpTransform = getPointerToJumpTransform('kouExpJumpTransform'), N.factors = 2, mod.type = 'standard', rtol = 1e-12, atol = 1e-28))

dvrg_c <- divergenceSwapRateCpp(p = c(0.1,0.5), coeffs = ode.solutions, stateMat = vols)

data.structure <- list(obs.data = obsData, spec.mat = expand.grid(t = mkt.spec$t, p = c(0.1,0.5), type = c("div","skew","quart")))

model.spec <- list(params.P = parList$P, params.Q = parList$Q, jump.type = 'kouExpJumpTransform', dt = 5/252, N.factors  = 2, error = list(cVec=rep(0.0,3), bVec = c(0.001,0.0020,0.0030)), mkt = mkt.spec)

lik.test <- modelLikelihood(data.structure = data.structure, model.spec = model.spec, for.estimation = F, filterFoo = divergenceModelR:::DSQ_sqrtFilter)

# lik.test.simpleFilter <- modelLikelihood(data.structure = data.structure, model.spec = model.spec, for.estimation = F, filterFoo = divergenceModelR:::DSQ_filter)

# ---- NOT RUN ----

# parList$P$jmp$muSc <- 1/parList$P$jmp$muSc
# parList$Q$jmp$muSc <- 1/parList$Q$jmp$muSc

# aa <- modelDynamics(params.P = parList$P, params.Q = parList$Q, dT = 5/252, N.factors = 2, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.points = 2, mod.type = 'standard', rtol = 1e-12, atol = 1e-30)

aa <- modelDynamics(params.P = parList$P, params.Q = parList$Q, dT = 5/252,N.factors = 2, jumpTransform = getPointerToJumpTransform('kouExpJumpTransform')$TF, N.points = 2, mod.type = 'standard', rtol = 1e-14, atol = 1e-30)

state <- matrix(c(8,7),1)

testMean <- meanVecFun(meanListS = aa$mean.vec, currVol = state)
testCovL <- covMatFun(covListS = aa$cov.list, covListDim = c(3,3) ,  currVol = state)
testCovAO <- affineOption::covMatFun(covListS = aa$cov.array[lower.tri(diag(3),diag=T)], U = 0, currVol = state)
testCov <- covMatFun(covListS = aa$cov.array[lower.tri(diag(3),diag=T)], covListDim = c(3,3) ,  currVol = state)

testCov <- testCov - testMean %*% t(testMean)
testBeta <- solve(testCov[-1,-1]) %*% testCov[1,-1]
stockCov <- testCov[1,1] - t(testBeta) %*% testCov[-1,-1] %*% (testBeta)

library(numDeriv)

numMean <- grad(func = function(uu){
  affineCF(u = matrix(c(uu[1],uu[2],uu[3]),nrow=1,ncol=3,byrow=T), params.Q = parList$Q, params.P = parList$P, t.vec = 5/252, v.0 = state, jumpTransform = getPointerToJumpTransform('kouExpJumpTransform')$TF, N.factors = 2, CGF = F)
}, x = c(0,0,0), method = "complex")
numCov <- hessian(func = function(uu){
  affineCF(u = matrix(c(uu[1],uu[2],uu[3]),nrow=1,ncol=3,byrow=T), params.Q = parList$Q, params.P = parList$P, t.vec = 5/252, v.0 = state, jumpTransform = getPointerToJumpTransform('kouExpJumpTransform')$TF, N.factors = 2, CGF = F, rtol = 1e-8)
},x = c(0,0,0), method = "complex", method.args = list(r=4, show.details=T))

numCovMat0 <- numCov - numMean %*% t(numMean)
eigen(numCovMat)$values


numMean <- grad(func = function(uu){
  affineCF(u = matrix(c(uu[1],uu[2],uu[3]),nrow=1,ncol=3,byrow=T), params.Q = parList$Q, params.P = parList$P, t.vec = 5/252, v.0 = state, jumpTransform = getPointerToJumpTransform('kouExpJumpTransform')$TF, N.factors = 2, CGF = F, rtol = 1e-11, atol = 1e-30)
}, x = c(0,0,0), method = "complex")
numCov <- hessian(func = function(uu){
  affineCF(u = matrix(c(uu[1],uu[2],uu[3]),nrow=1,ncol=3,byrow=T), params.Q = parList$Q, params.P = parList$P, t.vec = 5/252, v.0 = state, jumpTransform = getPointerToJumpTransform('kouExpJumpTransform')$TF, N.factors = 2, CGF = F, rtol = 1e-11, atol = 1e-30)
},x = c(0,0,0), method = "complex", method.args = list(r=4,show.details=F))



numCovMat <- numCov - numMean %*% t(numMean)
eigen(numCovMat)$values

numMean <- grad(func = function(uu){
  affineCF(u = matrix(c(uu[1],uu[2],uu[3]),nrow=1,ncol=3,byrow=T), params.Q = parList$Q, params.P = parList$P, t.vec = 5/252, v.0 = state, jumpTransform = getPointerToJumpTransform('kouExpJumpTransform')$TF, N.factors = 2, CGF = F)
}, x = c(0,0,0), method = "Richardson")
numCov <- hessian(func = function(uu){
  affineCF(u = matrix(c(uu[1],uu[2],uu[3]),nrow=1,ncol=3,byrow=T), params.Q = parList$Q, params.P = parList$P, t.vec = 5/252, v.0 = state, jumpTransform = getPointerToJumpTransform('kouExpJumpTransform')$TF, N.factors = 2, CGF = F, rtol=1e-6)
},x = c(0,0,0), method = "Richardson", method.args = list(r=4,show.details=T))

numCovMatR <- numCov - numMean %*% t(numMean)
eigen(numCovMat)$values
