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
#' @param spot.vol.series xts object: estimated spot volatility, if you don't like built-in methods; time stamps must correspond to time stamps in \code{rdata} upon input if you don't aggegate, otherwise they must yield a reasonable aggregated series.
#' @param jump.series xts object: estimated jump sizes, if you don't like built-in methods; provide jump times and sizes only, a substitute for the rdata series will be created. 
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
  
  # Adjustment for careless coders
  if(hasArg(data)){ rdata <- data }
  
  if(!(is.null(align.by) & is.null(align.period))){
    marketopen <- as.POSIXlt(x = marketopen, tz = "UTC", format = "%H:%M:%S")
    marketclose <- as.POSIXlt(x = marketclose, tz = "UTC", format = "%H:%M:%S")
    marketopen <- c(marketopen$hour, marketopen$min, marketopen$sec)
    marketclose <- c(marketclose$hour, marketclose$min, marketclose$sec)
    rdata <- aggregatePrice(ts = rdata, on = align.by, k = align.period, ...)
    # rdata <- aggregatePrice_Xts(rdata = rdata, period_ = align.by, numPeriods_ = align.period, dayStart_ = marketopen, dayEnd_ = marketclose)
  }
  
  if(makeReturns){
    rdata <- makeReturns(ts = rdata)
  }
  
  eval(parse(text = paste0("fooBase <- diveRgence:::",fooStr,"Base")))
  
  # Multiday adjustment: 
  multixts <- highfrequency:::.multixts(rdata);
  
  if(multixts){
    
    result <- apply.daily(rdata, fooBase, pow, align.by = NULL, align.period = NULL, intradaySeasonFun, makeReturns = makeReturns)
    return(result)
    
  } else if(!multixts){
    
    result <- fooBase(rdata, pow, align.by = NULL, align.period = NULL, intradaySeasonFun, makeReturns = makeReturns)
    return(result)
  }
  
}

#' @export rDivEngineInference
#' @describeIn rDivEngine

rDivEngineInference <- function(rdata, fooStr, pow, test.size = 0.05, align.by, align.period, makeReturns, reference.time, year.days = 365, seconds.per.day = 86400, cl = NULL, spot.vol.series = NULL, jump.series = NULL, kernel.type = "gaussian", vol.jumping = TRUE, ...){
  
  # Adjustment for careless coders
  if(hasArg(data)){ rdata <- data }
  
  
  if(!(is.null(align.by) | is.null(align.period))){
    rdata <- aggregatePrice(ts = rdata, on = align.by, k = align.period, ...)
    if(!is.null(spot.vol.series)){
      spot.vol.series$minus <- aggregatePrice(ts = spot.vol.series$minus, on = align.by, k = align.period, ...)
      spot.vol.series$plus <- aggregatePrice(ts = spot.vol.series$plus, on = align.by, k = align.period, ...) 
    }
  }
  
  dates <- unique(as.Date(index(rdata)))
  rdata.list <- lapply(dates, function(dd){
    loc.sp <- rdata[as.character(dd)]
    return(loc.sp)
  })
  
  ### Calculate realized divergence estimates
  rDiv.result <- rDivEngine(rdata = rdata, fooStr = fooStr, pow = pow, makeReturns = makeReturns, align.by = NULL, align.period = NULL)
  
  if(makeReturns){
    rdata <- makeReturns(ts = rdata)
  }
  
  # Get average volatility during the day -- this is for better truncation
  avg.vol <- rMPVcpp(rdata = rdata, mNum = 4, pPow = 2, yearDays = year.days)
  
  # annualize volatilities
  T.ranges <- t(sapply(rdata.list, function(x) range(index(x))))
  T.ranges <- apply(T.ranges,1,diff)
  T.ranges <- T.ranges / (year.days * 86400)
  avg.vol <- avg.vol / T.ranges
  
  ### Calculate variance of realized divergence
  rDiv.ci <- rDivEngineVarFun(rdata = rdata, fooStr = fooStr, pow = pow, test.size = test.size, makeReturns = makeReturns, align.by = NULL, align.period = NULL, avg.vol = avg.vol, spot.vol.series = spot.vol.series, jump.series = jump.series, intradaySeasonFun = function(x) 1, reference.time = reference.time, year.days = year.days, kernel.type = kernel.type, cl = cl, vol.jumping = vol.jumping, ...)
  
  ### return
  result <- list(rDiv = rDiv.result, rDiv.clt = rDiv.ci + rep(rDiv.result, ncol(rDiv.ci)))
  class(result) <- c("list","divinf")
  return(result)
}

rDivEngineVarFun <- function(rdata, fooStr, pow, test.size = test.size, align.by, align.period, makeReturns, avg.vol, spot.vol.series, jump.series, intradaySeasonFun, reference.time, year.days, kernel.type, vol.jumping = TRUE, cl = NULL, ...){
  # Adjustment for careless coders
  if(hasArg(data)){ rdata <- data }
  
  # Multiday adjustment: 
  multixts <- highfrequency:::.multixts(rdata);
  
  if(multixts){
    if(is.null(cl)){
      result <- apply.daily(rdata, rDivEngineVarFunBase, fooStr = fooStr, pow= pow, test.size = test.size, align.by= NULL, align.period = NULL, makeReturns = makeReturns, avg.vol = avg.vol, spot.vol.series = spot.vol.series, jump.series = jump.series, intradaySeasonFun = intradaySeasonFun, reference.time = reference.time, year.days, kernel.type = kernel.type, vol.jumping = vol.jumping, ...)  
    } else {
      rdata.dates <- unique(as.Date(index(rdata)))
      rdata.list <- lapply(rdata.dates,function(x) rdata[as.character(x)])
      clusterExport(cl, c("fooStr","pow","test.size","avg.vol","intradaySeasonFun","reference.time","year.days","spot.vol.series","jump.series"), envir = environment())
      if(any(grepl("snow",search()))) {
        res.list <- snow::parLapply(cl = cl, x = rdata.list, fun = rDivEngineVarFunBase, fooStr = fooStr, pow = pow, test.size =test.size, align.by = NULL, align.period = NULL, makeReturns = makeReturns, avg.vol = avg.vol, spot.vol.series = spot.vol.series, jump.series = jump.series, intradaySeasonFun = intradaySeasonFun, reference.time = reference.time, year.days = year.days, kernel.type = kernel.type, vol.jumping = vol.jumping)
      } else if(any(grepl("parallel",search()))){
        res.list <- parallel::parLapply(cl = cl, X = rdata.list, fun = rDivEngineVarFunBase, fooStr = fooStr, pow = pow, test.size =test.size, align.by = NULL, align.period = NULL, makeReturns = makeReturns, avg.vol = avg.vol, spot.vol.series = spot.vol.series, jump.series = jump.series, intradaySeasonFun = intradaySeasonFun, reference.time = reference.time, year.days = year.days, kernel.type = kernel.type, vol.jumping = vol.jumping)
      }
      result <- t(do.call(what = cbind, res.list))
      result <- xts(result, rdata.dates)
    }
    return(result)
  } else if(!multixts){
    result <- rDivEngineVarFunBase(rdata = rdata, fooStr = fooStr, pow = pow, test.size = test.size, align.by = align.by, align.period = align.period, makeReturns = makeReturns, avg.vol = avg.vol, spot.vol.series = spot.vol.series, jump.series = jump.series, intradaySeasonFun = intradaySeasonFun, reference.time = reference.time, year.days = year.days,kernel.type = kernel.type, vol.jumping = vol.jumping, ...)
    return(result)
  }
}

rDivEngineVarFunBase <- function(rdata, fooStr, pow, test.size, align.by, align.period, makeReturns, avg.vol, spot.vol.series, jump.series, intradaySeasonFun, reference.time, year.days, kernel.type, vol.jumping = TRUE, ...){
  
  if((!is.null(align.by))&&(!is.null(align.period))){
    rdata <- aggregatePrice(rdata, on=align.by, k=align.period, ...);
  }
  
  if(makeReturns){
    rdata <- tail(rdata,-1) 
  }  
  
  result <- apply(matrix(pow),1, rDivEngineVarFoo, fooStr = fooStr, tsMat = rdata, test.size = test.size, align.by = align.by, align.period = align.period, makeReturns = FALSE, avg.vol = avg.vol, spot.var = spot.vol.series, jump.series = jump.series, intradaySeasonFun = intradaySeasonFun, reference.time = reference.time, year.days = year.days, kernel.type = kernel.type, vol.jumping = vol.jumping, ...) 
  
  return(result)  
}

rDivEngineVarFoo <- function(p, fooStr, tsMat, avg.vol, spot.var, jump.series, test.size, intradaySeasonFun, reference.time, year.days, kernel.type, vol.jumping = TRUE,  ...){
  
  eval(parse(text = paste0("fooBaseDeriv <- ",fooStr,"BaseDeriv")))
  eval(parse(text = paste0("fooBaseZDeriv <- ",fooStr,"BaseZDeriv")))
  eval(parse(text = paste0("fooBaseContPart <- ",fooStr,"BaseContPart")))
  
  # where are we in avg vol?
  avg.vol <- avg.vol[as.character(head(index(tsMat),1), format = "%Y-%m-%d")]
  # spot var is hard at the ends of the sample. set k as function of single day length
  # k <- ceiling(sqrt(1/3 * nrow(tsMat)))
  k <- 2
  if(is.null(spot.var)){
    spot.var <- spotVol(rdata = tsMat, spot.index = index(tsMat), makeReturns = FALSE, align.by = NULL, align.period = NULL, avg.vol = avg.vol, reference.time = reference.time, vol.jumping = vol.jumping, year.days = year.days, kernel.type = kernel.type) 
    spot.var <- spot.var$spot[(k+1):(nrow(tsMat)-k)]
  } else {
    if(nrow(spot.var$minus) > nrow(tsMat)){
      spot.var$minus <- spot.var$minus[as.character(as.Date(index(tsMat)[1]))]
      spot.var$plus <- spot.var$plus[as.character(as.Date(index(tsMat)[1]))] 
      spot.var$minus <- head(spot.var$minus, -(nrow(spot.var$minus) - nrow(tsMat)))
      spot.var$plus <- tail(spot.var$plus, -(nrow(spot.var$plus) - nrow(tsMat)))
    } else {
      spot.var$minus <- spot.var$minus[as.character(as.Date(index(tsMat)[1]))]
      spot.var$plus <- spot.var$plus[as.character(as.Date(index(tsMat)[1]))] 
    }
  }
  
  if(!is.null(jump.series)){
    jump.series <- jump.series[as.character(as.Date(index(tsMat)[1]))]
    if(dim(jump.series)[1]==0){
      jump.series <- xts(x = rep(0,nrow(tsMat)), order.by = index(tsMat))
    } else {
      jump.series <- aggregatePrice(ts = jump.series, on = "hours", k = 7, aggr.vec = index(tsMat), pad = TRUE, pad.arg = 0, marketopen = format(head(index(tsMat),1), "%H:%M:%S"), marketclose = format(tail(index(tsMat)), "%H:%M:%S"))
    }
  }
  
  avg.var <- avg.vol
  
  # Vol seasonality
  reference.time <- as.POSIXlt(reference.time, format = "%H:%M:%S")
  reference.time$year <- as.POSIXlt(index(tsMat)[1])$year
  reference.time$mon <- as.POSIXlt(index(tsMat)[1])$mon
  reference.time$mday <- as.POSIXlt(index(tsMat)[1])$mday
  time.stamp <- sapply(index(tsMat), function(x) as.numeric(difftime(time1 = x, time2 = reference.time, units = "sec")))
  time.delta <- c(diff(time.stamp), median(diff(time.stamp))) / (86400 * year.days)
  
  # divergence derivative squared
  if(is.null(jump.series)){
    div.deriv <- apply(X = matrix(p), MARGIN = 1, FUN = fooBaseDeriv, tsMat = tsMat[(k+1):(nrow(tsMat)-k)], cumTsMat = cumsum(tsMat[(k+1):(nrow(tsMat)-k)]),.sum = FALSE) * as.numeric(abs(tsMat[(k+1):(nrow(tsMat)-k)]) > 2.5 * sqrt(as.numeric(avg.var) * intradaySeasonFun(time.stamp[(k+1):(nrow(tsMat)-k)])) * time.delta[(k+1):(nrow(tsMat)-k)]^(0.4999))    
  } else {
    div.deriv <- apply(X = matrix(p), MARGIN = 1, FUN = fooBaseDeriv, tsMat = as.numeric(jump.series), cumTsMat = cumsum(tsMat), .sum = FALSE)
  }

  
  # divergence derivative squared, with respect to z
  if(is.null(jump.series)){
    div.deriv.z <- apply(X = matrix(p), MARGIN = 1, FUN = fooBaseZDeriv, tsMat = tsMat[(k+1):(nrow(tsMat)-k)], cumTsMat = cumsum(tsMat[(k+1):(nrow(tsMat)-k)]), .sum = FALSE) * as.numeric(abs(tsMat[(k+1):(nrow(tsMat)-k)]) > 2.5 * sqrt(as.numeric(avg.var) * intradaySeasonFun(time.stamp[(k+1):(nrow(tsMat)-k)]))* time.delta[(k+1):(nrow(tsMat)-k)]^(0.4999))
  } else {
    div.deriv.z <- apply(X = matrix(p), MARGIN = 1, FUN = fooBaseZDeriv, tsMat = as.numeric(jump.series), cumTsMat = cumsum(tsMat), .sum = FALSE)
    }
  
  # For the continuous part 
  if(is.null(jump.series)){
    div.cont.part <- apply(X = matrix(p), MARGIN = 1, FUN = fooBaseContPart, tsMat = tsMat[(k+1):(nrow(tsMat)-k)], .sum = FALSE) * as.numeric(abs(tsMat[(k+1):(nrow(tsMat)-k)]) <= 2.5 * sqrt(as.numeric(avg.var) * intradaySeasonFun(time.stamp[(k+1):(nrow(tsMat)-k)]))* time.delta[(k+1):(nrow(tsMat)-k)]^(0.4999))
  } else {
    div.cont.part <- apply(X = matrix(p), MARGIN = 1, FUN = fooBaseContPart, tsMat = tsMat - as.numeric(jump.series), .sum = FALSE)
  }
  
  div.cont.part <- sum(div.cont.part)
  
  # Monte Carlo
  if(length(test.size) == 1){
    num.reps <- 1e3/test.size  
    test.size <- c(1-test.size/2, test.size/2)
  } else {
    num.reps <- min(1e4,1e3/min(min(test.size),diff(test.size)))
    test.size <- rev(test.size)
  }
  
  limiting.variable <- mcCltInference(rdivDerivX = div.deriv, rdivDerivZ = div.deriv.z, rdivCont = div.cont.part, spotVolPlus = as.numeric(spot.var$plus), spotVolMinus = as.numeric(spot.var$minus), nSampl = num.reps)
  limiting.variable <- limiting.variable * sqrt(median(time.delta))
  
  # Construct confidence interval
  res.ci <- - quantile(limiting.variable, probs = test.size)
  
  # multiply, sum, look what gets
  # spot.var <- matrix(as.numeric(spot.var), nrow = nrow(spot.var), ncol = ncol(div.deriv.sq))
  # result <- apply(spot.var*div.deriv.sq, 2, sum)
  
  return(res.ci)
}