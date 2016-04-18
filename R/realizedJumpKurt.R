#' @name realizedJumpKurt
#' @title Realized jump kurtosis divergence
#' @description Functions for calculating realized weighted univariate jump kurtosis divergences around divergence with power \code{pow}, and the variance of the estimates.
#' @param rdata an \code{xts} object containing 1 return series.
#' @param pow \code{numeric} vector of length \code{P}: jump skew divergences will be calculated around all powers.
#' @param makeReturns boolean, should be \code{TRUE} when price data is supplied. Defaults to \code{FALSE}.
#' @param align.by
#' @param align.period
#' @param ... Arguments passed on to \code{\link[highfrequency]{aggregatePrice}}
#' @details The most important arguments to pass go \code{\link[highfrequency]{aggregatePrice}} are \code{marketopen} and \code{marketclose}, see documentation therein. The default values are different from our test data set. \cr \cr
#' \code{rJKurt} calculates the realized kurtosis jump divergence. \cr \cr
#' \code{rJKurtInference} calls \code{rJKurt} and calculates the asymptotic variance of the estimates with the use of the feasible estimator. \cr \cr
#' \code{rJKurtTrueInference} calculates the realized measure and true asymptotic variance from simulated data, including jump data.
 
#' @export rJKurt

rJKurt <-  function(rdata, pow, makeReturns, align.by, align.period, intradaySeasonFun, ...){
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

rJKurtBase <- function(rdata, pow, align.by, align.period, makeReturns, intradaySeasonFun, ...){
  
  if((!is.null(align.by))&&(!is.null(align.period))){
    rdata <- aggregatePrice(rdata, on=align.by, k=align.period, ...);
  }
  
  if(makeReturns){
    rdata <- tail(makeReturns(rdata),-1)
  }  
  
  rdata <- as.numeric(rdata)
  
  result <- apply(matrix(pow),1,divKurtFoo, tsMat = rdata) 
  
}

divKurtFoo <- function(p, tsMat){
  if(p == 0){
    res <- 2*expm1(tsMat) - 2*tsMat - tsMat^2 - 1/3*tsMat^3
  } else if(p == 1){
    res <- 2 + 1/3*(expm1(tsMat)+1)*(-6 + 6*tsMat - 3*tsMat^2 + tsMat^3)
  } else {
    cnst <- 1/(p^3*(p-1)^3)
    res <- 2*(p-1)^3
    res <- res - 2*p^3*(expm1(tsMat)+1)
    res <- res + (expm1(p*tsMat)+1)*(2+(p-1)*p*(6+tsMat*(2+p*(-4+(p-1)*tsMat))))
    res <- cnst * res
    sm.tsMat <- which(abs(tsMat) < 1e-3)
    res[sm.tsMat] <- tsMat[sm.tsMat]^4/12 + tsMat[sm.tsMat]^5/60*(1+3*p)
  }
  res <- sum(res)
  return(res)
}

rJKurtBaseDeriv <- function(p, tsMat, .sum = FALSE){
  if(p == 0){
    res <- 2*expm1(tsMat) - 2*tsMat - tsMat^2
  } else if(p == 1){
    res <- exp(tsMat)*tsMat^3/3
  } else {
    cnst.1 <- 1/(p-1)^3
    res <- -2*cnst.1*exp(tsMat) + 2*cnst.1*exp(p*tsMat) + cnst.1*exp(p*tsMat)*tsMat*(2 - 2*p) + cnst.1*exp(p*tsMat)*tsMat^2*(1-2*p+p^2)
  }
  if(.sum){
    res <- sum(res)  
  }
  return(res)
}

rJKurtBaseZDeriv <- function(p, tsMat, .sum = FALSE){
  res <- rep(0, nrow(tsMat))
  if(.sum){
    res <- 0
  }
  return(res)
}

rJKurtBaseContPart <- function(p, tsMat, .sum = TRUE){
  
  res <- rep(0, nrow(tsMat))
  if(.sum){
    res <- 0
  }
  return(res)
  
}


# #' @export rJKurtTrueInference
# #' @describeIn realizedJumpKurt
# 
# rJKurtTrueInference <- function(rdata, pow, align.by, align.period, jump.data, var.data, makeReturns, ...){
#   
#   ### Calculate realized divergence estimates
#   rDiv.result <- rJKurt(rdata = rdata, pow = pow, makeReturns = makeReturns, align.by = align.by, align.period = align.period, ...)
#   
#   ### Calculate true CLT components
#   spot.var.at.jumps <- var.data[index(jump.data)]
#   
#   ### True realized quarticity
#   ## First determine dt, assuming a 6.5-hr trading day, i.e. 23400 seconds, 252 trading days per year
#   dt.base <- 1/23400/252
#   
#   ### True jump divergences at jump times
#   jmp.divs.true <- apply(X = matrix(pow), MARGIN = 1, FUN = divKurtDerivFoo, tsMat = as.numeric(jump.data))
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
#   divs.all <- apply(X = matrix(pow), MARGIN = 1, FUN = divKurtDerivFoo, tsMat = as.numeric(rdata))
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
# #' @export rJKurtInference
# #' @describeIn realizedJumpKurt
# 
# rJKurtInference <- function(rdata, pow, align.by, align.period, makeReturns, ...){
#   
#   ### Calculate realized divergence estimates
#   rDiv.result <- rJKurt(rdata = rdata, pow = pow, makeReturns = makeReturns, align.by = align.by, align.period = align.period, ...)
#   
#   ### Calculate variance of realized divergence
#   rDiv.variance <- rJKurtVarFun(rdata = rdata, pow = pow, makeReturns = makeReturns, align.by = align.by, align.period = align.period, ...)
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
# rJKurtVarFun <- function(rdata, pow, align.by, align.period, makeReturns, ...){
#   # Adjustment for careless coders
#   if(hasArg(data)){ rdata <- data }
#   
#   # Multiday adjustment: 
#   multixts <- highfrequency:::.multixts(rdata);
#   
#   if(multixts){
#     
#     result <- apply.daily(rdata, rJKurtVarFunBase, pow= pow, align.by= align.by, align.period = align.period, makeReturns = makeReturns, ...)
#     return(result)
#     
#   } else if(!multixts){
#     
#     result <- rJKurtVarFunBase(rdata, pow, align.by, align.period, makeReturns, ...)
#     return(result)
#   }
# }
# 
# rJKurtVarFunBase <- function(rdata, pow, align.by, align.period, makeReturns, ...){
#   
#   if((!is.null(align.by))&&(!is.null(align.period))){
#     rdata <- aggregatePrice(rdata, on=align.by, k=align.period, ...);
#   }
#   
#   if(makeReturns){
#     rdata <- makeReturns(rdata) 
#   }  
#   
#   result <- apply(matrix(pow),1, rJKurtVarFoo, tsMat = rdata, align.by = align.by, align.period = align.period, makeReturns = FALSE) 
#   
#   return(result)  
# }
# 
# rJKurtVarFoo <- function(p, tsMat, ...){
#   
#   # spot var
#   spot.var <- spotVol(rdata = tsMat, spot.index = index(tsMat), estFun = spotVolBaseJump, ...)
#   
#   # divergence derivative squared
#   div.deriv.sq <- apply(X = matrix(p), MARGIN = 1, FUN = divKurtDerivFoo, tsMat = tsMat, .sum = FALSE)^2
#   
#   # multiply, sum, look what gets
#   spot.var <- matrix(as.numeric(spot.var), nrow = nrow(spot.var), ncol = ncol(div.deriv.sq))
#   result <- apply(spot.var*div.deriv.sq, 2, sum)
#   
#   return(result)
# }