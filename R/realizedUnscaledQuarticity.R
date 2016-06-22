#' @name realizedUnscaledKurticity
#' @title Realized Unscaled Kurticity
#' @description Functions for calculating realized univariate divergences with power \code{pow}, and the variance of the estimates. Use string \code{"rUKurt"} to call via \code{\link{rDivEngine}} functions
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
#' @export rUKurt
#' 

rUKurt <-  function(rdata, pow, makeReturns, align.by, align.period, intradaySeasonFun, ...){
  # Adjustment for careless coders
  if(hasArg(data)){ rdata <- data }
  
  # Multiday adjustment: 
  multixts <- highfrequency:::.multixts(rdata);
  
  if(multixts){
    
    result <- apply.daily(rdata, rJKurtBase, pow, align.by, align.period, makeReturns, ...)
    return(result)
    
  } else if(!multixts){
    
    result <- rJKurtBase(rdata, pow, align.by, align.period, makeReturns, ...)
    return(result)
  }
}

rUKurtBase <- function(rdata, pow, align.by, align.period, makeReturns, intradaySeasonFun, ...){
  
  if((!is.null(align.by))&&(!is.null(align.period))){
    rdata <- aggregatePrice(rdata, on=align.by, k=align.period, ...);
  }
  
  if(makeReturns){
    rdata <- tail(rdata,-1)
  }  
  
  rdata <- as.numeric(rdata)
  
  result <- apply(matrix(pow),1,divUKurtFoo, tsMat = rdata) 
  
}

divUKurtFoo <- function(p, tsMat){
  z <- c(0,head(cumsum(tsMat),-1))
  if(p == 0){ 
    res <- -2 - tsMat^3/3 - 2*z - z^2 - tsMat^2*(1 + z) + exp(tsMat)* (2 + 2*z + z^2) - tsMat * (2 + 2*z + z^2)
  } else if(p == 1){
    res <- 1/3*exp(z)*(3*(2 - 2*z + z^2) + exp(tsMat)*(tsMat^3 + 3*tsMat^2*(-1 + z) - 3*(2 - 2*z + z^2) + 3*tsMat*(2 - 2*z + z^2)))
  } else {
    res <- -exp(tsMat)*p^3*(2 + (-1 + p)*z*(-2 + (-1 + p)*z))
    res <- res + (-1 + p)^3*(2 + p*z*(-2 + p*z))
    res <- res + exp(p*tsMat)*(2 + (-1 + p)*p*(6 + (-1 + p)*p*tsMat^2 + 2*tsMat*(1 - 2*p + (-1 + p)*p*z) + z*(2 + p*(-4 + (-1 + p)*z))))
    res <- res * exp(p*z)/(p^3 * (p-1)^3)
    # sm.tsMat <- which(abs(tsMat) < 1e-5)
    # res[sm.tsMat] <- tsMat[sm.tsMat]^3/6 + tsMat[sm.tsMat]^4/24*(1+2*p)
  }
  res <- sum(res)
  return(res)
}

divUKurtFoo_true <- function(p, tsMat, z){
  
  if(p == 0){ 
    res <- -2 - tsMat^3/3 - 2*z - z^2 - tsMat^2*(1 + z) + exp(tsMat)* (2 + 2*z + z^2) - tsMat * (2 + 2*z + z^2)
  } else if(p == 1){
    res <- 1/3*exp(z)*(3*(2 - 2*z + z^2) + exp(tsMat)*(tsMat^3 + 3*tsMat^2*(-1 + z) - 3*(2 - 2*z + z^2) + 3*tsMat*(2 - 2*z + z^2)))
  } else {
    res <- -exp(tsMat)*p^3*(2 + (-1 + p)*z*(-2 + (-1 + p)*z))
    res <- res + (-1 + p)^3*(2 + p*z*(-2 + p*z))
    res <- res + exp(p*tsMat)*(2 + (-1 + p)*p*(6 + (-1 + p)*p*tsMat^2 + 2*tsMat*(1 - 2*p + (-1 + p)*p*z) + z*(2 + p*(-4 + (-1 + p)*z))))
    res <- res * exp(p*z)/(p^3 * (p-1)^3)
    # sm.tsMat <- which(abs(tsMat) < 1e-5)
    # res[sm.tsMat] <- tsMat[sm.tsMat]^3/6 + tsMat[sm.tsMat]^4/24*(1+2*p)
  }
  res <- sum(res)
  return(res)
}

rUKurtBaseDeriv <- function(p, tsMat, .sum = FALSE){
  z <- c(0,head(cumsum(tsMat),-1))
  if(p == 0){
    res <- -2 - tsMat^2 - 2*z - z^2 - 2*tsMat*(1 + z) + exp(tsMat)*(2 + 2*z + z^2)
  } else if(p == 1){
    res <- 1/3*exp(tsMat + z)*tsMat*(tsMat^2 + 3*tsMat*z + 3*z^2)
  } else {
    res <- -exp(tsMat)*(2 + (-1 + p)*z*(-2 + (-1 + p)*z))
    res <- res + exp(p*tsMat)*(2 + (-1 + p)^2*tsMat^2 + (-1 + p)*z*(-2 + (-1 + p)*z) + 2*(-1 + p)*tsMat*(-1 + (-1 + p)*z))
  }
  if(.sum){
    res <- sum(res)  
  }
  return(res)
}

rUKurtBaseZDeriv <- function(p, tsMat, .sum = FALSE){
  z <- c(0,head(cumsum(tsMat),-1))
  if(p == 0){
    res <- -tsMat^2 - 2*(1 + z) + 2*exp(tsMat)*(1 + z) - 2*tsMat*(1 + z)
  } else if(p == 1){
    res <- 1/3*exp(z)*(3*z^2 + exp(tsMat)*(tsMat^3 + 3*tsMat^2*z - 3*z^2 + 3*tsMat* z^2))
  } else {
    res <- (-1 + p)^3*z^2
    res <- res + exp(p*tsMat)*(2 + (-1 + p)^2*tsMat^2 + (-1 + p)*z*(-2 + (-1 + p)*z) + 2*(-1 + p)*tsMat*(-1 + (-1 + p)*z))
    res <- res - exp(tsMat)*(2 + (-1 + p)*z*(-2 + (-1 + p)*p*z))
    res <- res * 1/(-1 + p)^3*exp(p*z)
  }
  if(.sum){
    res <- sum(res)  
  }
  return(res)
}

rUKurtBaseContPart <- function(p, tsMat, .sum = TRUE){
  
  z <- c(0,head(cumsum(tsMat),-1))
  
  # The continuous part here is integral of sigma_s^4 ds
  time.incr <- diff(as.integer(index(tsMat)))
  time.incr <- c(time.incr,median(time.incr))
  time.incr <- time.incr / (365*86400)
  res <- 1/3 * 1/2 * tsMat^4 / time.incr * z^4 * exp(2*p*z)
  if(!.sum){
    return(res)
  } else {
    return(sum(res))
  }
  
}