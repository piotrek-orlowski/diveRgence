# Test if estimated parameters fit better the empirical means of the portfolios than the true ones
data(bates.params)
data(simulEstParams)

dT <- 5/252

m.true <- NULL
m.est <- NULL
for (tt in c(1/12,1/2)) {
  for (uu in c(-2,3)) {
    port.struct <- hedged.return(u=uu,tau.j=dT/2,delta=1/252/78,ttm=tt)
    m.true <- c(m.true,momentCondition(bates.params.P,bates.params.Q,condition.struct=port.struct,conditional=FALSE))
    m.est <- c(m.est,momentCondition(par.est$P,par.est$Q,condition.struct=port.struct,conditional=FALSE))
  }
}

m.true <- 1/(1/252/78) * dT * m.true
m.est <- 1/(1/252/78) * dT * m.est

print(apply(hedgedRet[,-1],2,mean))
print(m.true)
print(m.est)

# now do simulation and check if conditional + simulation = unconditional
ab.lists.true <- meanAndcovList(expand.grid(u=c(-2,3),t=c(1/12,1/2)),params.P=bates.params.P, params.Q = bates.params.Q,volMean=TRUE,rtol.Q=1e-15,rtol=1e-12)
ab.lists.est <- meanAndcovList(expand.grid(u=c(-2,3),t=c(1/12,1/2)),params.P=par.est$P, params.Q = par.est$Q,volMean=TRUE,rtol.Q=1e-15,rtol=1e-12)

set.seed(1)

v.sim.est <- simulator2fStationaryDistDraw(par.est$P,sample.size=5e4)
v.sim.true <- simulator2fStationaryDistDraw(bates.params.P,sample.size=5e4)

mean.true <- mean.est <- 0

for (vv in 1:(nrow(v.sim.est))) {
  meanCov <- .Call("meanAndCovWrap",ab.lists.true$mean.list,ab.lists.true$cov.list,4,as.numeric(v.sim.true[vv,c("v1","v2")]))
  mean.true <- mean.true + meanCov$mean
  
  meanCov <- .Call("meanAndCovWrap",ab.lists.est$mean.list,ab.lists.est$cov.list,4,as.numeric(v.sim.est[vv,c("v1","v2")]))
  mean.est <- mean.est + meanCov$mean
}

mean.true <- (mean.true/nrow(v.sim.true))[2:5]
mean.est <- (mean.est/nrow(v.sim.est))[2:5]

expect_true(max(abs(c(mean.true/m.true,mean.est/m.est)-1))<1e-2)
