#' @name rDivEngine
#' @title Realized divergence inference wrapper.
#' @description A set of functions to calculate realized divergence measures and respective confidence intervals based on semimartingale discretisation theory.
#' @param rdata an \code{xts} object containing 1 price or return series.
#' @param fooStr character, name of base function (realized measure) for inference.
#' @param pow \code{numeric} vector of length \code{P}: functions from the set associated with \code{fooStr} will be evaluated with powers \code{pow}.
#' @param makeReturns boolean, should be \code{TRUE} when price data is supplied. Defaults to \code{FALSE}.
#' @param align.by argument to \code{\link[highfrequency]{aggregatePrice}}
#' @param align.period argument to \code{\link[highfrequency]{aggregatePrice}}
#' @param intradaySeasonFun Function. Allows to control for diurnal patterns in volatility whhen calculating estimators with jump truncation. Accepts 1 argument: time in seconds from start of trading day. The default setting returns 1.
#' @param reference.time Time string in format \code{\%H:\%M:\%S}, the time relative to which inputs to \code{intradaySeasonFun} are calculated.
#' @param ... Arguments passed on to \code{\link[highfrequency]{aggregatePrice}}
#' @details The most important arguments to pass go \code{\link[highfrequency]{aggregatePrice}} are \code{marketopen} and \code{marketclose}, see documentation therein. The default values are different from our test data set. \cr \cr
#' There following divergence types are available (see Khajavi, Orlowski, Trojani 2015):\cr
#' \code{fooStr = 'rDiv'} -- realized power divergence of log returns with \code{p=pow}, \cr
#' \code{fooStr = 'rUDiv'} -- realized power divergence, scaled by value at outset of period, \code{p=pow},\cr
#' \code{fooStr = 'rJSkew'} -- realized skewness divergence of log returns (jump skewness) around power \code{p=pow}\cr
#' \code{fooStr = 'rUSkew'} -- realized skewness divergence, scaled by value at outset of period, around power \code{p=pow}, (similar to signed realized volaility)\cr
#' \code{fooStr = 'rJKurt'} -- realized quarticity divergence of log returns (jump kurtosis) around power \code{p=pow}\cr
#' \code{fooStr = 'rUKurt'} -- realized quarticity divergence, scaled by value at outset of period, around power \code{p=pow}, (similar to realized volaility weighted by divergence of return from outset)\cr
#' 
#' @return \code{rDivEngine} returns an \code{xts} object of dimension \code{num.days x length(pow)} \cr \cr
#' \code{rDivEngineInference} returns a list with fields \code{rDiv} and \code{asy.var}; the former contains the output of \code{rDivEngine}, the latter contains the confidence interval for estimation error, i.e. for \eqn{\hat{D}-D}.
#' 
#' @export rDivEngine

rDivEngine <- function(rdata, fooStr, pow, makeReturns, align.by, align.period, marketopen = "08:30:00", marketclose= "15:15:00" , intradaySeasonFun = function(x) 1 , ...){
  
  if(!(is.null(align.by) & is.null(align.period))){
    marketopen <- as.POSIXlt(x = marketopen, tz = "UTC", format = "%H:%M:%S")
    marketclose <- as.POSIXlt(x = marketclose, tz = "UTC", format = "%H:%M:%S")
    marketopen <- c(marketopen$hour, marketopen$min, marketopen$sec)
    marketclose <- c(marketclose$hour, marketclose$min, marketclose$sec)
    rdata <- aggregatePrice_Xts(rdata = rdata, period_ = align.by, numPeriods_ = align.period, dayStart_ = marketopen, dayEnd_ = marketclose)
  }
  
  if(makeReturns){
    rdata <- makeReturns(ts = rdata)
  }
  
  eval(parse(text = paste0("fooBase <- diveRgence:::",fooStr,"Base")))
  
  # Adjustment for careless coders
  if(hasArg(data)){ rdata <- data }
  
  # Multiday adjustment: 
  multixts <- highfrequency:::.multixts(rdata);
  
  if(multixts){
    
    result <- apply.daily(rdata, fooBase, pow, align.by, align.period, makeReturns, intradaySeasonFun, ...)
    return(result)
    
  } else if(!multixts){
    
    result <- fooBase(rdata, pow, align.by, align.period, makeReturns, intradaySeasonFun, ...)
    return(result)
  }
  
}

#' @export rDivEngineInference
#' @describeIn rDivEngine

rDivEngineInference <- function(rdata, fooStr, pow, align.by, align.period, makeReturns, reference.time, year.days = 365, seconds.per.day = 86400, ...){
  
  if(!(is.null(align.by) | is.null(align.period))){
    rdata <- aggregatePrice(ts = rdata, on = align.by, k = align.period, ...)
  }
  
  if(makeReturns){
    rdata <- makeReturns(ts = rdata)
  }
  
  # Get average volatility during the day -- this is for better truncation
  avg.vol <- rMPVcpp(rdata = rdata, mNum = 4, pPow = 2, yearDays = year.days)
  
  dates <- unique(as.Date(index(rdata)))
  rdata.list <- lapply(dates, function(dd){
    loc.sp <- rdata[as.character(dd)]
    return(loc.sp)
  })
  
  # annualize volatilities
  T.ranges <- t(sapply(rdata.list, function(x) range(index(x))))
  T.ranges <- apply(T.ranges,1,diff)
  T.ranges <- T.ranges / (year.days * seconds.per.day)
  avg.vol <- avg.vol / T.ranges
  
  # volatility seasonality
  rdata.lengths <- sapply(rdata.list,length)
  rdata.unique.lenghts <- unique(rdata.lengths)
  rdata.counts <- sapply(rdata.unique.lenghts, function(x) length(which(rdata.lengths==x)))
  rdata.pick <- rdata.list[which(rdata.lengths == rdata.unique.lenghts[which.max(rdata.counts)])]
  rdata.pick.ranges <- t(sapply(rdata.list, function(x) range(index(x)-index(x)[1])))
  rdata.pick.ranges <- apply(rdata.pick.ranges,2,median)
  rdata.pick.ranges <- as.POSIXct(rdata.pick.ranges, origin = Sys.Date()) + as.integer(as.POSIXct(reference.time,format="%H:%M:%S")) - as.integer(as.POSIXct(Sys.Date()))
  rdata.pick.ranges <- as.character(rdata.pick.ranges, format = "%H:%M:%S")
  rdata.for.seasonality <- do.call(rbind.xts, rdata.pick)
  
  time.alignment.for.seasonality <- floor(median(unlist(sapply(rdata.for.seasonality, function(x) diff(as.numeric(index(x)))))))
  
  seasonVol <- tryCatch(highfrequency::spotvol(data = rdata.for.seasonality, makeReturns = T, method = "detper", on = "seconds", k = time.alignment.for.seasonality, marketopen = rdata.pick.ranges[1], marketclose = rdata.pick.ranges[2])$periodic, error = function(e){
    print(e)
    print("\n")
    print("Error in seasonality calculation, periodic component set to unity.")
    res <- xts(rep(1,rdata.unique.lenghts[which.max(rdata.counts)]), order.by = index(rdata.for.seasonality[as.character(index(rdata.for.seasonality[1]), format = "%Y-%m-%d")]))
    return(res)
  })
  seasonVol <- splinefun(x = as.numeric(sapply(index(seasonVol), function(x) difftime(x,index(seasonVol)[1],units="secs"))), y = as.numeric(seasonVol), method = "natural")
    
  ### Calculate realized divergence estimates
  rDiv.result <- rDivEngine(rdata = rdata, fooStr = fooStr, pow = pow, makeReturns = makeReturns, align.by = align.by, align.period = align.period, intradaySeasonFun = seasonVol, ...)
  
  ### Calculate variance of realized divergence
  rDiv.ci <- rDivEngineVarFun(rdata = rdata, fooStr = fooStr, pow = pow, makeReturns = makeReturns, align.by = align.by, align.period = align.period, avg.vol = avg.vol, intradaySeasonFun = seasonVol, reference.time = as.character(as.POSIXct(median(sapply(rdata.for.seasonality, function(ll) index(head(ll,1)))), origin="1970-01-01"), format= "%H:%M:%S"), ...)
  
#   ### Rescale with frequency
#   if(!is.null(align.period) & !is.null(align.by)){
#     stopifnot(align.by %in% c("minutes","seconds","hours"))
#     if(align.by == "minutes"){
#       
#       mean.dt <- align.period * 60 / (23400 * 365)
#       
#     } else if(align.by == "seconds"){
#       
#       mean.dt <- align.period / (23400 * 365)  
#       
#     } else if(align.by == "hours"){
#       
#       mean.dt <- align.period * 3600 / (23400 * 365)
#       
#     }
#   } else {
#     
#     T.range <- apply.daily(x = rdata, FUN = function(x) range(index(x)))
#     T.range <- T.range[,2] - T.range[,1]
#     mean.dt <- apply.daily(x = rdata, FUN = function(x) mean(difftime(time2 = head(index(rdata),-1), time1 = tail(index(rdata),-1), units = "secs")))
#     mean.dt <- mean.dt/T.range/365
#   }
#   
#   rDiv.variance <- rDiv.variance * mean.dt
  
  ### return
  result <- list(rDiv = rDiv.result, asy.var = rDiv.ci)
  class(result) <- c("list","divinf")
  return(result)
}

rDivEngineVarFun <- function(rdata, fooStr, pow, align.by, align.period, makeReturns, avg.vol, intradaySeasonFun, reference.time, ...){
  # Adjustment for careless coders
  if(hasArg(data)){ rdata <- data }
  
  # Multiday adjustment: 
  multixts <- highfrequency:::.multixts(rdata);
  
  if(multixts){
    
    result <- apply.daily(rdata, rDivEngineVarFunBase, fooStr = fooStr, pow= pow, align.by= align.by, align.period = align.period, makeReturns = makeReturns, avg.vol = avg.vol, intradaySeasonFun = intradaySeasonFun, reference.time = reference.time, ...)
    return(result)
    
  } else if(!multixts){
    
    result <- rDivEngineVarFunBase(rdata, fooStr, pow, align.by, align.period, makeReturns, avg.vol = avg.vol, intradaySeasonFun = intradaySeasonFun, reference.time = reference.time, ...)
    return(result)
  }
}

rDivEngineVarFunBase <- function(rdata, fooStr, pow, align.by, align.period, makeReturns, avg.vol, intradaySeasonFun, reference.time, ...){
  
  if((!is.null(align.by))&&(!is.null(align.period))){
    rdata <- aggregatePrice(rdata, on=align.by, k=align.period, ...);
  }
  
  if(makeReturns){
    rdata <- makeReturns(rdata) 
  }  
  
  result <- apply(matrix(pow),1, rDivEngineVarFoo, fooStr = fooStr, tsMat = rdata, align.by = align.by, align.period = align.period, makeReturns = FALSE, avg.vol = avg.vol, intradaySeasonFun = intradaySeasonFun, reference.time = reference.time, ...) 
  
  return(result)  
}

rDivEngineVarFoo <- function(p, fooStr, tsMat, avg.vol, intradaySeasonFun, reference.time, ...){
  
  eval(parse(text = paste0("fooBaseDeriv <- ",fooStr,"BaseDeriv")))
  eval(parse(text = paste0("fooBaseZDeriv <- ",fooStr,"BaseZDeriv")))
  eval(parse(text = paste0("fooBaseContPart <- ",fooStr,"BaseContPart")))
  
  # where are we in avg vol?
  avg.vol <- avg.vol[as.character(head(index(tsMat),1), format = "%Y-%m-%d")]
  # spot var is hard at the ends of the sample. set k as function of single day length
  k <- ceiling(sqrt(1/3 * nrow(tsMat)))
  spot.var <- spotVol(rdata = tsMat, spot.index = index(tsMat), estFun = spotVolBaseJump, makeReturns = F, align.by = NULL, align.period = NULL, avg.vol = avg.vol, intradaySeasonFun = intradaySeasonFun, reference.time = reference.time)
  
  avg.var <- spot.var$avg
  spot.var <- spot.var$spot[(k+1):(nrow(tsMat)-k)]
  
  # Vol seasonality
  reference.time <- as.POSIXlt(reference.time, format = "%H:%M:%S")
  reference.time$year <- as.POSIXlt(index(tsMat)[1])$year
  reference.time$mon <- as.POSIXlt(index(tsMat)[1])$mon
  reference.time$mday <- as.POSIXlt(index(tsMat)[1])$mday
  time.stamp <- sapply(index(tsMat), function(x) as.numeric(difftime(time1 = x, time2 = reference.time, units = "sec")))
  time.delta <- c(diff(time.stamp), median(diff(time.stamp))) / (86400 * 365)
  
  # divergence derivative squared
  div.deriv <- apply(X = matrix(p), MARGIN = 1, FUN = fooBaseDeriv, tsMat = tsMat[(k+1):(nrow(tsMat)-k)], .sum = FALSE) * as.numeric(abs(tsMat[(k+1):(nrow(tsMat)-k)]) > 3 * sqrt(as.numeric(avg.var) * intradaySeasonFun(time.stamp[(k+1):(nrow(tsMat)-k)])) * time.delta[(k+1):(nrow(tsMat)-k)]^(0.49) )
  
  # divergence derivative squared, with respect to z
  div.deriv.z <- apply(X = matrix(p), MARGIN = 1, FUN = fooBaseZDeriv, tsMat = tsMat[(k+1):(nrow(tsMat)-k)], .sum = FALSE) * as.numeric(abs(tsMat[(k+1):(nrow(tsMat)-k)]) > 3 * sqrt(as.numeric(avg.var) * intradaySeasonFun(time.stamp[(k+1):(nrow(tsMat)-k)]))* time.delta[(k+1):(nrow(tsMat)-k)]^(0.49))
  
  # For the continuous part 
  div.cont.part <- apply(X = matrix(p), MARGIN = 1, FUN = fooBaseContPart, tsMat = tsMat[(k+1):(nrow(tsMat)-k)], .sum = FALSE) * as.numeric(abs(tsMat[(k+1):(nrow(tsMat)-k)]) <= 3 * sqrt(as.numeric(avg.var) * intradaySeasonFun(time.stamp[(k+1):(nrow(tsMat)-k)]))* time.delta[(k+1):(nrow(tsMat)-k)]^(0.49))
  
  div.cont.part <- rxsumCpp(div.cont.part)
  
  # Monte Carlo
  test.size <- 0.05
  num.reps <- 1e3/test.size
  
  # Generate random variates
#   U.plus <- matrix(rnorm(n = (nrow(tsMat)-2*(k)) * num.reps, mean = 0, sd = 1), nrow = nrow(tsMat)-2*k, ncol= num.reps)
#   U.minus <- matrix(rnorm(n = (nrow(tsMat)-2*(k)) * num.reps, mean = 0, sd = 1), nrow = nrow(tsMat)-2*k, ncol= num.reps)
#   K.p <- matrix(runif(n = (nrow(tsMat) - 2*(k)) * num.reps, min = 0, max = 1), nrow = nrow(tsMat)-2*k, ncol= num.reps)
#   U.cont <- rnorm(num.reps, mean = 0, sd = 1)
#   
#   # Calculate limiting variables
#   limiting.variable <- U.plus * sqrt(K.p) * matrix(sqrt(spot.var$minus), nrow = nrow(tsMat) - 2*k, ncol = num.reps) * matrix(div.deriv, nrow = nrow(tsMat) - 2*k, ncol = num.reps)
#   rm(U.plus)
#   limiting.variable <- limiting.variable + U.minus * sqrt(1-K.p) * matrix(sqrt(spot.var$plus), nrow = nrow(tsMat) - 2*k, ncol = num.reps) * matrix(div.deriv, nrow = nrow(tsMat) - 2*k, ncol = num.reps)
#   limiting.variable <- limiting.variable - U.minus * sqrt(K.p) * matrix(div.deriv.z, nrow = nrow(tsMat) - 2*k, ncol = num.reps)
#   rm(U.minus, K.p)
#   limiting.variable <- apply(X = limiting.variable, MARGIN = 2, FUN = rxsumCpp)
#   limiting.variable <- limiting.variable + U.cont * sqrt(div.cont.part)
#   limiting.variable <- limiting.variable * sqrt(median(time.delta))
  limiting.variable <- mcCltInference(rdivDerivX = div.deriv, rdivDerivZ = div.deriv.z, rdivCont = div.cont.part, spotVolPlus = spot.var$plus, spotVolMinus = spot.var$minus, nSampl = num.reps)
  limiting.variable <- limiting.variable * sqrt(median(time.delta))
  
  # Construct confidence interval
  res.ci <- - quantile(limiting.variable, probs = c(1-test.size/2, test.size/2))
  
  # multiply, sum, look what gets
  # spot.var <- matrix(as.numeric(spot.var), nrow = nrow(spot.var), ncol = ncol(div.deriv.sq))
  # result <- apply(spot.var*div.deriv.sq, 2, rxsumCpp)
  
  return(res.ci)
}