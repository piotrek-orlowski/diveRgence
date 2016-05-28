#' @name realizedUnscaledDivergence
#' @title Realized Unscaled Divergence
#' @description Functions for calculating realized univariate divergences with power \code{pow}, and the variance of the estimates. Use string \code{"rUDiv"} to call via \code{\link{rDivEngine}} functions
#' @param rdata an \code{xts} object containing 1 return series.
#' @param pow \code{numeric} vector of length \code{P}: divergences will be calculated at all powers.
#' @param makeReturns boolean, should be \code{TRUE} when price data is supplied. Defaults to \code{FALSE}.
#' @param align.by
#' @param align.period
#' @param ... Arguments passed on to \code{\link[highfrequency]{aggregatePrice}}
#' @details The most important arguments to pass go \code{\link[highfrequency]{aggregatePrice}} are \code{marketopen} and \code{marketclose}, see documentation therein. The default values are different from our test data set. \cr
#' \cr
#' \code{rDivInference} calculates the estimator and asymptotic variance based on the feasible estimator in Khajavi, Or{\l}owski and Trojani (2015). \cr \cr
#' \code{rDivTrueInference} calculates the estimator and the infeasible asymptotic variance based on jump data from simulation.
#' @export rUDiv

rUDiv <- function(rdata, pow, makeReturns, align.by, align.period, intradaySeasonFun, ...){
  
  # Adjustment for careless coders
  if(hasArg(data)){ rdata <- data }
  
  # Multiday adjustment: 
  multixts <- highfrequency:::.multixts(rdata);
  
  if(multixts){
    
    result <- apply.daily(rdata, rDivBase, pow, align.by, align.period, makeReturns, ...)
    return(result)
    
  } else if(!multixts){
    
    result <- rDivBase(rdata, pow, align.by, align.period, makeReturns, ...)
    return(result)
  }
}

rUDivBase <- function(rdata, pow, align.by, align.period, makeReturns, intradaySeasonFun, ...){
  
  if((!is.null(align.by))&&(!is.null(align.period))){
    rdata <- aggregatePrice(rdata, on=align.by, k=align.period, ...);
  }
  
  if(makeReturns){
    rdata <- tail(rdata,-1)
  }  
  
  rdata <- as.numeric(rdata)
  
  result <- apply(matrix(pow),1,uDivFoo, tsMat = rdata) 
  return(result)
}

uDivFoo <- function(p, tsMat){
  z <- c(0,head(cumsum(tsMat),-1))
  if(p == 0){
    res <- - tsMat + (exp(tsMat) - 1)
  } else if(p == 1){
    res <- exp(tsMat)*tsMat - (exp(tsMat) - 1)
  } else {
    cnst.1 <- 1/(p*(p-1))
    cnst.2 <- 1/(p-1)
    res <- cnst.1 * (exp(p*tsMat) - 1)
    res <- res - cnst.2*(exp(tsMat) - 1)
  }
  res <- exp(p*z)*res
  res <- sum(res)
  return(res)
}

rUDivBaseDeriv <- function(p, tsMat, .sum = FALSE){
  z <- c(0,head(cumsum(tsMat),-1))
  if(p == 0){
    res <- exp(tsMat)-1
  } else if(p == 1){
    res <- exp(tsMat)*tsMat
  } else {
    cnst.1 <- 1/(p-1)
    res <- cnst.1 * (exp(p*tsMat) - exp(tsMat))
  }
  res <- exp(p*z)*res
  if(.sum){
    res <- sum(res)  
  }
  return(res)
}

rUDivBaseZDeriv <- function(p, tsMat, .sum = FALSE){
  z <- c(0,head(cumsum(tsMat),-1))
  if(p == 0){
    res <- - tsMat + (exp(tsMat) - 1)
  } else if(p == 1){
    res <- exp(tsMat)*tsMat - (exp(tsMat) - 1)
  } else {
    cnst.1 <- 1/(p*(p-1))
    cnst.2 <- 1/(p-1)
    res <- cnst.1 * (exp(p*tsMat) - 1)
    res <- res - cnst.2*(exp(tsMat) - 1)
  }
  res <- p*exp(p*z)*res
  res <- sum(res)
  return(res)
}

rUDivBaseContPart <- function(p, tsMat, .sum = TRUE){
  
  z <- c(0,head(cumsum(tsMat),-1))
  
  # The continuous part here is integral of sigma_s^4 ds
  time.incr <- diff(as.integer(index(tsMat)))
  time.incr <- c(time.incr,median(time.incr))
  time.incr <- time.incr / (365*86400)
  res <- 1/3 * 1/2 * tsMat^4 / time.incr * exp(2*p*z)
  if(!.sum){
    return(res)
  } else {
    return(sum(res))
  }
  
}