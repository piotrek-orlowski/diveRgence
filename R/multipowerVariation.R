#' @name multipowerVariation
#' @title Multipower variation
#' @description Calculate multipower variation, an estimator of integrated volatiltiy. See summary in Veraart (2010) and Jacod (2006).
#' @param rdata an \code{xts} object containing 1 return series.
#' @param makeReturns boolean, should be \code{TRUE} when price data is supplied. Defaults to \code{FALSE}.
#' @param align.by
#' @param align.period
#' @param r.tot which power-variation to estimate.
#' @param num.int how many neighbouring intervals are taken into consideration.
#' @param ... Arguments passed on to \code{\link[highfrequency]{aggregatePrice}}
#' @details The most important arguments to pass go \code{\link[highfrequency]{aggregatePrice}} are \code{marketopen} and \code{marketclose}, see documentation therein. The default values are different from our test data set. \cr \cr
#' Choose \code{r.tot} and \code{num.int} such that \code{r.tot/num.int < 2}.

#' @export rMPV

rMPV <- function(rdata, r.tot, num.int, makeReturns, year.days, align.by, align.period, ...){
  # Adjustment for careless coders
  if(hasArg(data)){ rdata <- data }
  
  if(!(is.null(align.by) | is.null(align.period))){
    rdata <- aggregatePrice(ts = rdata, on = align.by, k = align.period, ...)
  }
  
  if(makeReturns){
    rdata <- makeReturns(ts = rdata)
  }
  
  rmpv <- rMPVcpp(rdata = rdata, mNum = num.int, pPow = r.tot, yearDays = year.days)
  
  index(rmpv) <- as.Date(index(rmpv))
  return(rmpv)
}

rMPV_legacy <- function(rdata, r.tot, num.int, makeReturns, align.by, align.period, ...){
  # Adjustment for careless coders
  if(hasArg(data)){ rdata <- data }
  
  # Multiday adjustment: 
  multixts <- highfrequency:::.multixts(rdata);
  
  if(multixts){
    
    result <- apply.daily(rdata, rMPVbase, r.tot = r.tot, num.int = num.int, align.by = align.by, align.period = align.period, makeReturns = makeReturns, ...)
    return(result)
    
  } else if(!multixts){
    
    result <- rMPVbase(rdata, r.tot, num.int, align.by, align.period, makeReturns, ...)
    return(result)
  }
}

rMPVbase <- function(rdata, r.tot, num.int, align.by, align.period, makeReturns, ...){
  
  if((!is.null(align.by))&&(!is.null(align.period))){
    rdata <- aggregatePrice(rdata, on=align.by, k=align.period, ...);
  }
  
  if(makeReturns){
    rdata <- makeReturns(rdata) 
  }
  
  # Constants for scaling
  r.loc <- r.tot/num.int
  mu.scale <- pi^(-0.5) * 2^(r.loc/2) * gamma((r.loc+1)/2)
  
  # delta-t
  if(!is.null(align.period) & !is.null(align.by)){
    stopifnot(align.by %in% c("minutes","seconds","hours"))
    if(align.by == "minutes"){
      
      mean.dt <- align.period * 60 / (23400 * 252)
      
    } else if(align.by == "seconds"){
      
      mean.dt <- align.period / (23400 * 252)  
      
    } else if(align.by == "hours"){
      
      mean.dt <- align.period * 3600 / (23400 * 252)
      
    }
  } else {
    
    T.range <- apply.daily(x = rdata, FUN = function(x) range(index(x)))
    T.range <- T.range[,2] - T.range[,1]
    mean.dt <- apply.daily(x = rdata, FUN = function(x) mean(difftime(time2 = head(index(rdata),-1), time1 = tail(index(rdata),-1), units = "secs")))
    mean.dt <- mean.dt/T.range/252
  }
  mean.dt.r <- mean.dt^(1-r.tot/2)
  
  mpvFoo <- function(ts, r){
    ts <- as.numeric(ts)
    ts <- prod(abs(ts)^r)
    return(ts)
  }
  
  rmpv <- rollapply(data = rdata, width = num.int, FUN = mpvFoo, r = r.loc, fill = 0)
  rmpv <- length(rdata)/(length(rdata) - num.int + 1) * mean.dt.r * mu.scale^(-num.int) * sum(as.numeric(rmpv))
  
  return(rmpv)
}