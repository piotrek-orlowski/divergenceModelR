load("D:/git/hfAffine/RpackageDev/affineOption/data/heston.params.RData")

aa <- modelDynamics(params.P = parListHeston$P, params.Q = parListHeston$Q, dT = 5/252, N.factors = 1, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.points = 2, mod.type = 'standard')

testMean <- meanVecFun(meanListS = aa$mean.vec, currVol = matrix(1,1,1))
testCov <- covMatFun(covListS = aa$cov.list, covListDim = c(2,2) ,  currVol = matrix(1,1,1))
testCov <- covMatFun(covListS = aa$cov.array[c(1,2,4)], covListDim = c(2,2) ,  currVol = matrix(1,1,1))

calc.par <- list(mean.vec = aa$mean.vec[-1], cov.array = aa$cov.array[2,2])

meancov <- affineTransitionStateHandler(stateMat = matrix(c(0.5,1,2),nrow=1), modelParameters = calc.par)

#### ---- OBSERVATION ----
mkt.spec <- data.frame(p=1,q=0,r=0,t=c(0.5,1))

sol.check.derivs <- Re(odeExtSolveWrap(u = cbind(0.5,0), params.Q = parListHeston$Q, mkt = mkt.spec, jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, mod.type = 'standard', rtol = 1e-12, atol = 1e-22))

obsParameters <- list(stockParams = list(mean.vec = aa$mean.vec, cov.array = aa$cov.array[lower.tri(diag(2) ,diag = T)]), cfCoeffs = sol.check.derivs, tVec = mkt.spec$t, pVec = 0.5, cVec = rep(0.8,3), bVec = rep(0.05,3))

stobs <- affineObservationStateHandler(stateMat = matrix(c(1,0.5,2), nrow=1), modelParameters = obsParameters)
