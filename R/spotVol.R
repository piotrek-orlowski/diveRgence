#' @name spotVolatility
#' @title Spot volatility estimation
#' @description Function for calculating spot volatility under the assumptions of (a) continuous underlying process, (b) noisy observations of a continuous underlying process, (c) underlying process with jumps. The \code{apply.at} and \code{apply.period} parameters indicate at which time stamps the estimator should be applied. If you specify to \code{align.by='seconds'}, \code{align.period=5}, then pick, for example, \code{apply.at='minutes'}, \code{apply.by=10}. For more details, see Mancini, Matiussi and Reno (2012).
#' @param rdata an \code{xts} object containing 1 return series.
#' @param spot.index POSIXct vector of times when volatility should be estimated. If this argument is given, \code{apply.at} and \code{apply.period} are ignored.
#' @param makeReturns boolean, should be \code{TRUE} when price data is supplied. Defaults to \code{FALSE}.
#' @param align.by
#' @param align.period
#' @param apply.at at what time unit intervals should the estimator be applied?
#' @param apply.period every how many time unit intervals should the estimator be applied?
#' @param estFun local volatility estimating function, one of \code{spotVolBase*} specified in this document.
#' @param nbar two-scale (two-frequency) parameter for eliminating microstructure noise in estimation.
#' @param ... Arguments passed on to \code{\link[highfrequency]{aggregatePrice}}
#' @details The most important arguments to pass go \code{\link[highfrequency]{aggregatePrice}} are \code{marketopen} and \code{marketclose}, see documentation therein. The default values are different from our test data set.
#' @export spotVol

spotVol <-  function(rdata, spot.index = NULL, makeReturns, align.by, align.period, apply.at, apply.period, avg.vol, reference.time= "07:30:00", vol.jumping = TRUE, year.days = 365, seconds.per.day = 86400, kernel.type = "gaussian", ...){
  # Adjustment for careless coders
  if(hasArg(data)){ rdata <- data }
  
  # Make default spot.index
  if(is.null(spot.index)){
    spot.index <- index(rdata)[tail(endpoints(x = rdata, on = apply.at, k = apply.period),-1)]
  }
  
  # Aggregate price time series
  if(!(is.null(align.by) | is.null(align.period))){
    rdata <- diveRgence::aggregatePrice(ts = rdata, on = align.by, k = align.period, ...)
  }
  
  if(makeReturns){
    spot.index <- tail(spot.index,-1)
    rdata <- makeReturns(ts = rdata)
  }
    
  # We only use data within a single day for each spot estimation
  spot.days <- as.character(unique(as.Date(spot.index)))
  rdata <- rdata[spot.days]
  
  # preparatory calculations
  T.range <- index(rdata[c(1,nrow(rdata))])
  T.range <- as.numeric(difftime(T.range[2],T.range[1], units = "secs"))
  delta.bar <- T.range / (length(rdata)-1)
  delta.bar.years <- (delta.bar/(86400*year.days))
  
  # Time stamps
  reference.time <- as.POSIXlt(reference.time, format = "%H:%M:%S")
  reference.time$year <- as.POSIXlt(spot.index[1])$year
  reference.time$mon <- as.POSIXlt(spot.index[1])$mon
  reference.time$mday <- as.POSIXlt(spot.index[1])$mday
  time.stamp <- as.numeric(index(rdata)) - as.numeric(reference.time)
  time.stamp.years <- time.stamp / (86400*year.days)
  
  # kernel time for bandwidth
  delta.bandwidth <- delta.bar.years
  
  # rdata transformations
  rdataSq <- as.numeric(rdata)^2
  rdataAbs <- abs(as.numeric(rdata))
  
  # Apply estimator
  # result <- sapply(X = spot.index, FUN = estFun, rdataSq = rdataSq, rdataAbs = rdataAbs, rdataInd = index(rdata), T.range = T.range, time.stamp.years = time.stamp.years, avg.vol = avg.vol, intradaySeasonFun = intradaySeasonFun, reference.time = reference.time, sep.lr = TRUE, time.delta = delta.bandwidth)
  result <- sapply(X = as.numeric(spot.index), FUN = diveRgence:::spotVolBaseJump_cpp, rdataSq = rdataSq, rdataAbs = rdataAbs, rdataInd = as.numeric(index(rdata)), tRange = as.numeric(T.range), timeStampYears = time.stamp.years, avgVol = as.numeric(avg.vol), referenceTime = as.numeric(reference.time), sepLR = vol.jumping, timeDelta = delta.bandwidth, yearLength = year.days, kernelType = kernel.type)
  
  result <- t(result)
  
  result <- xts(x = result, order.by = spot.index)
  colnames(result) <- c("minus","plus")
  result <- list(spot = result, avg = avg.vol)
  return(result)
}
# 
# #' @export spotVolBaseContinuous
# 
# spotVolBaseContinuous <- function(spot.point, rdata, align.by, align.period, makeReturns, ...){
# #   
# #   if(makeReturns){
# #     rdata <- makeReturns(rdata) 
# #   }  
# #   
#   # Get range of day, multiply by 3600 to get seconds
#   #T.range <- as.numeric(diff(range(index(rdata)))*3600)
#   T.range <- range(index(rdata))
#   T.range <- as.numeric(difftime(T.range[2],T.range[1], units = "secs"))
#   delta.bar <- T.range / (length(rdata)-1)
#   
#   # Find a kernel
#   delta.bandwidth <- exp(1.4)*length(rdata)^0.62 * delta.bar/3
#   
#   # Calculate estimator
#   time.kernel <- expDoubleKernel(x0 = spot.point, cnv.data = index(rdata), delta.bandwidth = delta.bandwidth)
#   time.kernel <- time.kernel / sum(time.kernel * delta.bar)
#   result <- sum(time.kernel * rdata^2* delta.bar)
#   result <- result / (delta.bar * 1/365 * 1/T.range)
#   
#     result <- xts(x = result, order.by = spot.point)
#   
#   return(result)
# }
# 
# #' @export spotVolBaseJump
# #' 
# spotVolBaseJump <- function(spot.point, rdataSq, rdataAbs, rdataInd, T.range, time.stamp.years, avg.vol, intradaySeasonFun, reference.time, sep.lr, time.delta){
#   
#   # Calculate kernel
#   time.kernel <- diveRgence:::kernel_gaussian(time.stamp.years, (as.numeric(spot.point)-as.numeric(reference.time))/(86400*365), time.delta/10)
#   
#   if(sep.lr){
#     # Separate left and right sides if necessary
#     time.kernel.minus <- time.kernel * as.numeric(rdataInd < spot.point)
#     time.kernel.plus <- time.kernel * as.numeric(rdataInd >= spot.point)
#     
#     result.minus <- sum(time.kernel.minus * rdataSq * as.numeric((rdataAbs < (3 * sqrt(as.numeric(avg.vol)*intradaySeasonFun(time.stamp)) * time.delta^0.49)))) / sum(head(time.kernel.minus,-1) * diff(time.stamp.years))
#     result.plus <- sum(time.kernel.plus * rdataSq * as.numeric((rdataAbs < (3 * sqrt(as.numeric(avg.vol)*intradaySeasonFun(time.stamp)) * time.delta^0.49)))) / sum(tail(time.kernel.plus,-1) * diff(time.stamp.years))
#   } else {
#     result.minus <- sum(time.kernel * rdataSq * as.numeric((rdataAbs < (3 * sqrt(as.numeric(avg.vol)*intradaySeasonFun(time.stamp)) * time.delta^0.49)))) / sum(head(time.kernel,-1) * diff(time.stamp.years))
#     result.plus <- result.minus
#   }
#   
#   result.plus <- xts(x = result.plus, order.by = spot.point)
#   result.minus <- xts(x = result.minus, order.by = spot.point)  
#   
#   return(rbind(minus = result.minus, plus = result.plus))
# }
# 
# #' @export spotVolBaseNoise
# 
# spotVolBaseNoise <- function(spot.point, rdata, nbar, align.by, align.period, makeReturns, ...){
#   
#   if(makeReturns){
#     rdata.nbar <- diff(x = log(rdata), lag = nbar)
#     rdata.nbar <- xts(x = as.numeric(rdata.nbar[-(1:nbar)]), order.by = index(head(rdata.nbar,-nbar)))
#     rdata <- makeReturns(rdata) 
#     rdata <- rdata[index(rdata.nbar)]
#   } else {
#     rdata.nbar <- rollsum(x = rdata, k = nbar)
#   }
#   
#   # Get range of day, multiply by 3600 to get seconds
#   # T.range <- as.numeric(diff(range(index(rdata)))*3600)
#   T.range <- range(index(rdata))
#   T.range <- as.numeric(difftime(T.range[2],T.range[1], units = "secs"))
#   delta.bar <- T.range / (length(rdata)-1)
#   delta.bar.years <- (delta.bar/(T.range*365))^0.95 # in business time
#   
#   # Find a kernel
#   delta.bandwidth <- exp(1.4)*length(rdata)^0.62 * delta.bar/3
#   
#   # Calculate estimator
#   time.kernel <- expDoubleKernel(x0 = spot.point, cnv.data = index(rdata), delta.bandwidth = delta.bandwidth)
#   time.kernel <- time.kernel / sum(time.kernel * delta.bar)
#   result <- 1/nbar * sum(time.kernel * (rdata.nbar^2-rdata^2) * delta.bar)
#   result <- result / (delta.bar * 1/365 * 1/T.range)
#   
#   result <- xts(x = result, order.by = spot.point)
#   
#   return(result)
# }