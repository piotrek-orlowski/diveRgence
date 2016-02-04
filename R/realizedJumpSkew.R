#' @name realizedJumpSkew
#' @title Realized jump skew divergence
#' @description Functions for calculating realized weighted univariate jump skew divergences around divergence with power \code{pow}, and the variance of the estimates.
#' @param rdata an \code{xts} object containing 1 return series.
#' @param pow \code{numeric} vector of length \code{P}: jump skew divergences will be calculated around all powers.
#' @param makeReturns boolean, should be \code{TRUE} when price data is supplied. Defaults to \code{FALSE}.
#' @param align.by
#' @param align.period
#' @param ... Arguments passed on to \code{\link[highfrequency]{aggregatePrice}}
#' @details The most important arguments to pass go \code{\link[highfrequency]{aggregatePrice}} are \code{marketopen} and \code{marketclose}, see documentation therein. The default values are different from our test data set. \cr \cr
#' \code{rJSkew} calculates the realized skew jump divergence. \cr \cr
#' \code{rJSkewInference} calls \code{rJSkew} and calculates the asymptotic variance of the estimates with the use of the feasible estimator. \cr \cr
#' \code{rJSkewTrueInference} calculates the realized measure and true asymptotic variance from simulated data, including jump data.

#' @export rJSkew

rJSkew <-  function(rdata, pow, makeReturns, align.by, align.period, intradaySeasonFun, ...){
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

rJSkewBase <- function(rdata, pow, align.by, align.period, makeReturns, intradaySeasonFun, ...){
  
  if((!is.null(align.by))&&(!is.null(align.period))){
    rdata <- aggregatePrice(rdata, on=align.by, k=align.period, ...);
  }
  
  if(makeReturns){
    rdata <- tail(makeReturns(rdata),-1)
  }  
  
  rdata <- as.numeric(rdata)
  
  result <- apply(matrix(pow),1,divSkewFoo, tsMat = rdata) 
  
}

divSkewFoo <- function(p, tsMat){
  if(p == 0){
    res <- -(tsMat + 0.5*tsMat^2) + (exp(tsMat) - 1)
  } else if(p == 1){
    res <- 0.5*exp(tsMat)*(2 - 2*tsMat + tsMat^2) - 1
  } else {
    cnst <- 1/(p^2*(p-1)^2)
    res <- 2*p-1
    res <- res + p^2*(exp(tsMat)-1)
    res <- res + exp(p*tsMat)*(1 - 2*p + (p - 1)*p*tsMat)
    res <- cnst * res
    sm.tsMat <- which(abs(tsMat) < 1e-5)
    res[sm.tsMat] <- tsMat[sm.tsMat]^3/6 + tsMat[sm.tsMat]^4/24*(1+2*p)
  }
  res <- rxsumCpp(res)
  return(res)
}

rJSkewBaseDeriv <- function(p, tsMat, .sum = FALSE){
  if(p == 0){
    res <- expm1(tsMat) - tsMat
  } else if(p == 1){
    res <- exp(tsMat)*tsMat^2/2
  } else {
    cnst.1 <- 1/(p-1)^2
    cnst.2 <- 1/(p-1)
    res <- exp(tsMat)*cnst.1 - exp(p*tsMat)*cnst.1 + exp(p*tsMat)*tsMat*cnst.2
  }
  if(.sum){
    res <- rxsumCpp(res)  
  }
  return(res)
}

rJSkewBaseZDeriv <- function(p, tsMat, .sum = FALSE){
  res <- rep(0, nrow(tsMat))
  if(.sum){
    res <- 0
  }
  return(res)
}

rJSkewBaseContPart <- function(p, tsMat, .sum = TRUE){
  
  res <- rep(0, nrow(tsMat))
  if(.sum){
    res <- 0
  }
  return(res)
  
}

# @export rJSkewTrueInference
# @describeIn realizedJumpSkew

# rJSkewTrueInference <- function(rdata, pow, align.by, align.period, jump.data, var.data, makeReturns, ...){
#   
#   ### Calculate realized divergence estimates
#   rDiv.result <- rJSkew(rdata = rdata, pow = pow, makeReturns = makeReturns, align.by = align.by, align.period = align.period, ...)
#   
#   ### Calculate true CLT components
#   spot.var.at.jumps <- var.data[index(jump.data)]
#   
#   ### True realized quarticity
#   ## First determine dt, assuming a 6.5-hr trading day, i.e. 23400 seconds, 252 trading days per year
#   dt.base <- 1/23400/252
#   
#   ### True jump divergences at jump times
#   jmp.divs.true <- apply(X = matrix(pow), MARGIN = 1, FUN = divSkewDerivFoo, tsMat = as.numeric(jump.data))
#   jmp.divs.true <- jmp.divs.true^2
#   jmp.divs.true <- xts(x = jmp.divs.true, order.by = index(jump.data))
#   
#   ### compile inference results: for every day take quarticity plus eventual jump cotributions
#   jmp.var <- 1/2 * matrix(as.numeric(jmp.divs.true), nrow = nrow(jmp.divs.true), ncol = ncol(jmp.divs.true)) * matrix(as.numeric(spot.var.at.jumps), nrow = nrow(spot.var.at.jumps), ncol = ncol(jmp.divs.true)) * 2
#   jmp.var <- xts(x = jmp.var, order.by = index(jmp.divs.true))
#   
#   ## Aggregate jmp.var within days
#   jmp.var <- apply.daily(x = jmp.var, FUN = colSums)
#   
#   cnt.var <- xts(x = rep(0,nrow(rDiv.result)), order.by = index(rDiv.result))
#   ## sum
#   inf.result <- cnt.var
#   for(nn in 1:ncol(inf.result)){
#     inf.result[as.character(as.Date(index(jmp.var))),nn] <- inf.result[as.character(as.Date(index(jmp.var))),nn] + as.numeric(jmp.var[,nn])
#   }
#   
#   ### calculate asymptotic variance with true spot vol and all increments (semi-true asy variance)
#   if(makeReturns){
#     rdata <- makeReturns(rdata)
#   }
#   divs.all <- apply(X = matrix(pow), MARGIN = 1, FUN = divSkewDerivFoo, tsMat = as.numeric(rdata))
#   divs.all <- divs.all^2
#   jmp.var.all <- matrix(as.numeric(divs.all), nrow = nrow(divs.all), ncol = ncol(divs.all)) * matrix(as.numeric(var.data), nrow = nrow(var.data), ncol = ncol(divs.all))
#   jmp.var.all <- xts(x = jmp.var.all, order.by = index(var.data))
#   
#   jmp.var.all <- apply.daily(x = jmp.var.all, FUN = colSums)
#   
#   ### Rescale with frequency
#   if(!is.null(align.period) & !is.null(align.by)){
#     stopifnot(align.by %in% c("minutes","seconds","hours"))
#     if(align.by == "minutes"){
#       
#       mean.dt <- align.period * 60 / (23400 * 252)
#       
#     } else if(align.by == "seconds"){
#       
#       mean.dt <- align.period / (23400 * 252)  
#       
#     } else if(align.by == "hours"){
#       
#       mean.dt <- align.period * 3600 / (23400 * 252)
#       
#     }
#   } else {
#     
#     T.range <- apply.daily(x = rdata, FUN = function(x) range(index(x)))
#     T.range <- T.range[,2] - T.range[,1]
#     mean.dt <- apply.daily(x = rdata, FUN = function(x) mean(difftime(time2 = head(index(rdata),-1), time1 = tail(index(rdata),-1), units = "secs")))
#     mean.dt <- mean.dt/T.range/252
#   }
#   
#   inf.result <- inf.result * mean.dt
#   jmp.var.all <- jmp.var.all * mean.dt
#   
#   ### return list
#   result.list <- list(rDiv = rDiv.result, asy.var = inf.result, asy.var.semi = jmp.var.all)
#   class(result.list) <- c("list","divinf")
#   return(result.list)
# }
# 
# #' @export rJSkewInference
# #' @describeIn realizedJumpSkew
# 
# rJSkewInference <- function(rdata, pow, align.by, align.period, makeReturns, ...){
#   
#   ### Calculate realized divergence estimates
#   rDiv.result <- rJSkew(rdata = rdata, pow = pow, makeReturns = makeReturns, align.by = align.by, align.period = align.period, ...)
#   
#   ### Calculate variance of realized divergence
#   rDiv.variance <- rJSkewVarFun(rdata = rdata, pow = pow, makeReturns = makeReturns, align.by = align.by, align.period = align.period, ...)
#   
#   ### Rescale with frequency
#   if(!is.null(align.period) & !is.null(align.by)){
#     stopifnot(align.by %in% c("minutes","seconds","hours"))
#     if(align.by == "minutes"){
#       
#       mean.dt <- align.period * 60 / (23400 * 252)
#       
#     } else if(align.by == "seconds"){
#       
#       mean.dt <- align.period / (23400 * 252)  
#       
#     } else if(align.by == "hours"){
#       
#       mean.dt <- align.period * 3600 / (23400 * 252)
#       
#     }
#   } else {
#     
#     T.range <- apply.daily(x = rdata, FUN = function(x) range(index(x)))
#     T.range <- T.range[,2] - T.range[,1]
#     mean.dt <- apply.daily(x = rdata, FUN = function(x) mean(difftime(time2 = head(index(rdata),-1), time1 = tail(index(rdata),-1), units = "secs")))
#     mean.dt <- mean.dt/T.range/252
#   }
#   
#   rDiv.variance <- rDiv.variance * mean.dt
#   
#   ### return
#   result <- list(rDiv = rDiv.result, asy.var = rDiv.variance)
#   class(result) <- c("list","divinf")
#   return(result)
# }
# 
# rJSkewVarFun <- function(rdata, pow, align.by, align.period, makeReturns, ...){
#   # Adjustment for careless coders
#   if(hasArg(data)){ rdata <- data }
#   
#   # Multiday adjustment: 
#   multixts <- highfrequency:::.multixts(rdata);
#   
#   if(multixts){
#     
#     result <- apply.daily(rdata, rJSkewVarFunBase, pow= pow, align.by= align.by, align.period = align.period, makeReturns = makeReturns, ...)
#     return(result)
#     
#   } else if(!multixts){
#     
#     result <- rJSkewVarFunBase(rdata, pow, align.by, align.period, makeReturns, ...)
#     return(result)
#   }
# }
# 
# rJSkewVarFunBase <- function(rdata, pow, align.by, align.period, makeReturns, ...){
#   
#   if((!is.null(align.by))&&(!is.null(align.period))){
#     rdata <- aggregatePrice(rdata, on=align.by, k=align.period, ...);
#   }
#   
#   if(makeReturns){
#     rdata <- makeReturns(rdata) 
#   }  
#   
#   result <- apply(matrix(pow),1, rJSkewVarFoo, tsMat = rdata, align.by = align.by, align.period = align.period, makeReturns = FALSE) 
#   
#   return(result)  
# }
# 
# rJSkewVarFoo <- function(p, tsMat, ...){
#   
#   # spot var
#   spot.var <- spotVol(rdata = tsMat, spot.index = index(tsMat), estFun = spotVolBaseJump, ...)
#   
#   # divergence derivative squared
#   div.deriv.sq <- apply(X = matrix(p), MARGIN = 1, FUN = divSkewDerivFoo, tsMat = tsMat, .sum = FALSE)^2
#   
#   # multiply, sum, look what gets
#   spot.var <- matrix(as.numeric(spot.var), nrow = nrow(spot.var), ncol = ncol(div.deriv.sq))
#   result <- apply(spot.var*div.deriv.sq, 2, rxsumCpp)
#   
#   return(result)
# }