#' @name realizedUnscaledSkew
#' @title Realized Unscaled Skewness
#' @description Functions for calculating realized univariate divergences with power \code{pow}, and the variance of the estimates. Use string \code{"rUSkew"} to call via \code{\link{rDivEngine}} functions
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
#' @export rUSkew
#' 

rUSkew <-  function(rdata, pow, makeReturns, align.by, align.period, intradaySeasonFun, ...){
# Adjustment for careless coders
if(hasArg(data)){ rdata <- data }

# Multiday adjustment: 
multixts <- highfrequency:::.multixts(rdata);

if(multixts){
  
  result <- apply.daily(rdata, rJSkewBase, pow, align.by, align.period, makeReturns, ...)
  return(result)
  
} else if(!multixts){
  
  result <- rJSkewBase(rdata, pow, align.by, align.period, makeReturns, ...)
  return(result)
}
}

rUSkewBase <- function(rdata, pow, align.by, align.period, makeReturns, intradaySeasonFun, ...){
  
  if((!is.null(align.by))&&(!is.null(align.period))){
    rdata <- aggregatePrice(rdata, on=align.by, k=align.period, ...);
  }
  
  if(makeReturns){
    rdata <- tail(rdata,-1)
  }  
  
  rdata <- as.numeric(rdata)
  
  result <- apply(matrix(pow),1,divUSkewFoo, tsMat = rdata) 
  
}

divUSkewFoo <- function(p, tsMat){
  z <- c(0,head(cumsum(tsMat),-1))
  if(p == 0){ 
    res <- -1 - tsMat^2/2 - z + exp(tsMat) * (1 + z) - tsMat * (1 + z)
  } else if(p == 1){
    res <- tsMat^2/2*exp(tsMat + z) + tsMat * exp(tsMat + z)*(z - 1) + expm1(tsMat)*(exp(z) - z*exp(z))
  } else {
    res <- 1/(p^2*(p-1)^2)
    res <- res * exp(p*z)
    res <- res * (-p^2*exp(tsMat)*(z*(p-1)-1) + (p*z-1)*(p-1)^2 + exp(p*tsMat)*(1 + p^2*(z+tsMat)-p*(2+z+tsMat)))
    # sm.tsMat <- which(abs(tsMat) < 1e-5)
    # res[sm.tsMat] <- tsMat[sm.tsMat]^3/6 + tsMat[sm.tsMat]^4/24*(1+2*p)
  }
  res <- sum(res)
  return(res)
}

divUSkewFoo_true <- function(p, tsMat, z){
  if(p == 0){ 
    res <- -1 - tsMat^2/2 - z + exp(tsMat) * (1 + z) - tsMat * (1 + z)
  } else if(p == 1){
    res <- tsMat^2/2*exp(tsMat + z) + tsMat * exp(tsMat + z)*(z - 1) + expm1(tsMat)*(exp(z) - z*exp(z))
  } else {
    res <- 1/(p^2*(p-1)^2)
    res <- res * exp(p*z)
    res <- res * (-p^2*exp(tsMat)*(z*(p-1)-1) + (p*z-1)*(p-1)^2 + exp(p*tsMat)*(1 + p^2*(z+tsMat)-p*(2+z+tsMat)))
    # sm.tsMat <- which(abs(tsMat) < 1e-5)
    # res[sm.tsMat] <- tsMat[sm.tsMat]^3/6 + tsMat[sm.tsMat]^4/24*(1+2*p)
  }
  res <- sum(res)
  return(res)
}


rUSkewBaseDeriv <- function(p, tsMat, .sum = FALSE){
  z <- c(0,head(cumsum(tsMat),-1))
  if(p == 0){
    res <- -1 - tsMat - z + exp(tsMat)*(1 + z)
  } else if(p == 1){
    res <- 1/2*exp(tsMat + z) *  tsMat *  (tsMat + 2*z)
  } else {
    res <- exp(p*z)/(p-1)^2
    res <- res * (exp(p*tsMat)*(-1 + tsMat*(p-1) + z * (p-1)) + exp(tsMat)*(1+z-p*z))
  }
  if(.sum){
    res <- sum(res)  
  }
  return(res)
}

rUSkewBaseZDeriv <- function(p, tsMat, .sum = FALSE){
  z <- c(0,head(cumsum(tsMat),-1))
  if(p == 0){
    res <- expm1(tsMat) - tsMat
  } else if(p == 1){
    res <- 1/2*exp(z)*(2*z + exp(tsMat)*(tsMat^2 - 2*z + 2*tsMat*z))
  } else {
    d_left_outer <- 1/(p^2*(p-1)^2) * p * exp(p*z)
    d_right_outer <- d_left_outer / p
    d_left_inner <- (-p^2 * exp(tsMat)*((p-1)*z-1) + (p*z-1)*(p-1)^2 + exp(p*tsMat)*(1 + p^2*(z+tsMat) - p*(z + tsMat + 2)))
    d_right_inner <- (-p^2 * exp(tsMat)*(p-1) + p*(p-1)^2 + exp(p*tsMat)*(p^2 - p))
    res <- d_left_outer * d_left_inner + d_right_outer * d_right_inner
  }
  if(.sum){
    res <- sum(res)  
  }
  return(res)
}

rUSkewBaseContPart <- function(p, tsMat, .sum = TRUE){
  
  z <- c(0,head(cumsum(tsMat),-1))
  
  # The continuous part here is integral of sigma_s^4 ds
  time.incr <- diff(as.integer(index(tsMat)))
  time.incr <- c(time.incr,median(time.incr))
  time.incr <- time.incr / (365*86400)
  res <- 1/3 * 1/2 * tsMat^4 / time.incr * z^2 * exp(2*p*z)
  if(!.sum){
    return(res)
  } else {
    return(sum(res))
  }
  
}