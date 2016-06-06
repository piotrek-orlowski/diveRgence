#### ---- LOAD LIBRARIES ----

library(diveRgence)
library(affineModelR)
Sys.setenv(TZ="UTC")
#### ---- LOAD PARAMETERS ----
data("affine-parameters-bates2006")

#### ---- SIMULATE TIME SERIES ----

sim.paths <- affineSimulate(paramsList = par.bates.svj1, N.factors = 1, t.days = 2*252, t.freq = 1/86400, freq.subdiv = 1, rng.seed = 1221, jumpGeneratorPtr = getPointerToGenerator('expNormJumpTransform'), jumpTransformPtr = getPointerToJumpTransform('expNormJumpTransform')$TF, mod.type = 'standard')

#### ---- AGGREGATE TIME SERIES to 1-MIN OBSERVATIONS ----

calendar <- sim.paths$times
calendar.posix <- as.POSIXct(as.Date("2016-01-01"))
calendar.posix <- calendar.posix + 8.5 * 3600
calendar.posix <- calendar.posix + 86400 * calendar

stock.path <- xts(x = sim.paths$sim.arrays[[1]]$S.array, order.by = calendar.posix)

stock.path <- diveRgence::aggregatePrice(ts = stock.path, on = "minutes", k = 1, marketopen = "08:30:00", marketclose = "15:15:00", pad = TRUE, pad.arg = NA_real_)

stock.path <- do.call(what = rbind.xts, args = lapply(as.character(unique(as.Date(index(stock.path)))), function(x) makeReturns(stock.path[x])[-1]))

vol.path <- xts(x = sim.paths$sim.arrays[[1]]$V.array, order.by = calendar.posix)
vol.path <- diveRgence::aggregatePrice(ts = vol.path, on = "minutes", k = 1, marketopen = "08:30:00", marketclose = "15:15:00", pad = TRUE, pad.arg = NA_real_)

jmp.cal <- sim.paths$jumpSizes[[1]][,1]
jmp.cal <- as.POSIXct(as.POSIXct(as.Date("2016-01-01")) + 8.5*3600 + 86400 * jmp.cal)
jmp.path <- xts(x = sim.paths$jumpSizes[[1]][,-1], order.by = jmp.cal)
jmp.path <- jmp.path["T08:30:00/T15:00:00"]

last.date <- tail(as.character(unique(as.Date(index(stock.path)))))

stock.path <- stock.path[paste0("/",last.date)]
vol.path <- vol.path[paste0("/",last.date)]

#### ---- ESTIMATE INTEGRATED VOLATILITY ----
library(highfrequency)

# jump-robust estimator of Andersen, Dobrev, Schaumburg (2012)
stock.rv <- highfrequency::rCov(rdata = stock.path, align.by = NULL, align.period = NULL, makeReturns = FALSE)

stock.rv <- 252 * 24/6.75 * stock.rv

vol.path.aggr <- apply.daily(vol.path, function(x){
  mean(x)
})

jmp.ind <- which(as.character(as.Date(index(vol.path.aggr))) %in% as.character(as.Date(index(jmp.path))))

vol.path.aggr[jmp.ind] <- vol.path.aggr[jmp.ind] + 406 * as.numeric(jmp.path[,1])^2

plot(sqrt(stock.rv), ylim = range(sqrt(c(vol.path.aggr,stock.rv))))
lines(sqrt(stock.rv), col = "darkblue", lwd = 1.5)
lines(sqrt(vol.path.aggr), col = "darkorange", lwd = 1.5)

#### ---- BOOTSTRAP CONFIDENCE INTERVALS ----
library(boot)
library(foreach)
library(iterators)

stock.dates <- unique(as.Date(index(stock.path)))

system.time(
  rv.ci <- foreach(dt = iter(stock.dates), .combine = rbind, .multicombine = TRUE) %do% {
    loc.ret <- stock.path[as.character(dt)]
    rv.boot <- tsboot(tseries = loc.ret, statistic = function(x) sqrt(252 * 24/6.75 * rCov(rdata = x, align.by = NULL, align.period = NULL, makeReturns = FALSE)), R = 1.5e3, sim = 'fixed', l = 30, orig.t = T)
    rv.ci <- boot.ci(boot.out = rv.boot, conf = 0.95)
    return(xts(x = rv.ci$perc[,4:5,drop=FALSE], order.by = dt))
    # return(xts(x = rv.ci$norm[,2:3,drop=FALSE], order.by = dt))
  }
)

#### ---- IN PARALLEL ----
library(parallel)
library(doParallel)

cl <- makeCluster(4)
registerDoParallel(cl)

clusterEvalQ(cl, library(highfrequency))
clusterEvalQ(cl, library(boot))

clusterExport(cl, 'stock.path')
system.time(
  rv.ci.par <- foreach(dt = iter(stock.dates), .combine = rbind, .multicombine = TRUE) %dopar% {
    loc.ret <- stock.path[as.character(dt)]
    rv.boot <- tsboot(tseries = loc.ret, statistic = function(x) sqrt(252 * 24/6.75 * rCov(rdata = x, align.by = NULL, align.period = NULL, makeReturns = FALSE)), R = 1.5e3, sim = "geom", l = 6, orig.t = T)
    rv.ci <- boot.ci(boot.out = rv.boot, conf = 0.95)
    return(xts(x = rv.ci$perc[,4:5,drop=FALSE], order.by = dt))
  }
)

index(stock.rv) <- as.Date(index(stock.rv))
index(vol.path.aggr) <- as.Date(index(vol.path.aggr))

plot(sqrt(stock.rv), ylim = range(c(sqrt(c(as.numeric(vol.path.aggr),as.numeric(stock.rv))),as.numeric(rv.ci.par))))
lines(sqrt(stock.rv), col = "darkblue", lwd = 1.5)
lines(sqrt(vol.path.aggr), col = "darkorange", lwd = 1.5)
lines(rv.ci.par[,1], col = "red", lwd = 1.5, lty = 2)
lines(rv.ci.par[,2], col = "red", lwd = 1.5, lty = 2)

closeAllConnections()
