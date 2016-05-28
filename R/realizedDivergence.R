#' @name realizedDivergence
#' @title Realized Divergence
#' @description Functions for calculating realized weighted univariate divergences with power \code{pow}, and the variance of the estimates.
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
#' @export rDiv

rDiv <- function(rdata, pow, makeReturns, align.by, align.period, intradaySeasonFun, ...){
  
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

rDivBase <- function(rdata, pow, align.by, align.period, makeReturns, intradaySeasonFun, ...){
  
  if((!is.null(align.by))&&(!is.null(align.period))){
    rdata <- aggregatePrice(rdata, on=align.by, k=align.period, ...);
  }
  
  if(makeReturns){
    rdata <- tail(rdata,-1)
  }  
  
  rdata <- as.numeric(rdata)
  
  result <- apply(matrix(pow),1,divFoo, tsMat = rdata) 
  return(result)
}

divFoo <- function(p, tsMat){
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
  res <- sum(res)
  return(res)
}

rDivBaseDeriv <- function(p, tsMat, .sum = FALSE){
  if(p == 0){
    res <- exp(tsMat)-1
  } else if(p == 1){
    res <- exp(tsMat)*tsMat
  } else {
    cnst.1 <- 1/(p-1)
    res <- cnst.1 * (exp(p*tsMat) - exp(tsMat))
  }
  if(.sum){
    res <- sum(res)  
  }
  return(res)
}

rDivBaseZDeriv <- function(p, tsMat, .sum = FALSE){
  res <- rep(0, nrow(tsMat))
  if(.sum){
    res <- 0
  }
  return(res)
}

rDivBaseContPart <- function(p, tsMat, .sum = TRUE){
  
  # The continuous part here is integral of sigma_s^4 ds
  time.incr <- diff(as.integer(index(tsMat)))
  time.incr <- c(time.incr,median(time.incr))
  time.incr <- time.incr / (365*86400)
  res <- 1/3 * 1/2 * tsMat^4 / time.incr
  if(!.sum){
    return(res)
  } else {
    return(sum(res))
  }
  
}

# #' @export rDivTrueInference
# #' @describeIn realizedDivergence
# 
# rDivTrueInference <- function(rdata, pow, align.by, align.period, jump.data, var.data, makeReturns, ...){
#   
#   ### Calculate realized divergence estimates
#   rDiv.result <- rDiv(rdata = rdata, pow = pow, makeReturns = makeReturns, align.by = align.by, align.period = align.period, ...)
#   
#   ### Calculate true CLT components
#   spot.var.at.jumps <- var.data[index(jump.data)]
#   
#   ### True realized quarticity
#   ## First determine dt, assuming a 6.5-hr trading day, i.e. 23400 seconds, 252 trading days per year
#   dt.base <- 1/23400/252
#   RQ.true <- apply.daily(x = var.data, FUN = function(v.ts){
#     dt.loc <- as.numeric(diff(index(v.ts))) * dt.base
#     v.ts <- head(v.ts,-1)
#     rq.loc <- sum(v.ts^2 * dt.loc)
#     return(rq.loc)
#   })
#   
#   ### True jump divergences at jump times
#   jmp.divs.true <- apply(X = matrix(pow), MARGIN = 1, FUN = divDerivFoo, tsMat = as.numeric(jump.data))
#   jmp.divs.true <- jmp.divs.true^2
#   jmp.divs.true <- xts(x = jmp.divs.true, order.by = index(jump.data))
#   
#   ### compile inference results: for every day take quarticity plus eventual jump cotributions
#   jmp.var <- 1/2 * matrix(as.numeric(jmp.divs.true), nrow = nrow(jmp.divs.true), ncol = ncol(jmp.divs.true)) * matrix(as.numeric(spot.var.at.jumps), nrow = nrow(spot.var.at.jumps), ncol = ncol(jmp.divs.true)) * 2
#   jmp.var <- xts(x = jmp.var, order.by = index(jmp.divs.true))
#   cnt.var <- 1/2 * RQ.true
#   cnt.var <- matrix(cnt.var, nrow = nrow(cnt.var), ncol = ncol(jmp.divs.true))
#   cnt.var <- xts(cnt.var, order.by = index(RQ.true))
#   
#   ## Aggregate jmp.var within days
#   jmp.var <- apply.daily(x = jmp.var, FUN = colSums)
#   
#   ## sum
#   inf.result <- cnt.var
#   for(nn in 1:ncol(inf.result)){
#     inf.result[as.character(as.Date(index(jmp.var))),nn] <- inf.result[as.character(as.Date(index(jmp.var))),nn] + as.numeric(jmp.var[,nn])
#   }
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
#    
#   ### return list
#   result.list <- list(rDiv = rDiv.result, asy.var = inf.result)
#   
#   class(result.list) <- c("list","divinf")
#   
#   return(result.list)
# }
# 
# #' @export rDivInference
# #' @describeIn realizedDivergence
# 
# rDivInference <- function(rdata, pow, align.by, align.period, makeReturns, ...){
# 
#   ### Calculate realized divergence estimates
#   rDiv.result <- rDiv(rdata = rdata, pow = pow, makeReturns = makeReturns, align.by = align.by, align.period = align.period, ...)
#   
#   ### Calculate variance of realized divergence
#   rDiv.variance <- rDivVarFun(rdata = rdata, pow = pow, makeReturns = makeReturns, align.by = align.by, align.period = align.period, ...)
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
# rDivVarFun <- function(rdata, pow, align.by, align.period, makeReturns, ...){
#   # Adjustment for careless coders
#   if(hasArg(data)){ rdata <- data }
#   
#   # Multiday adjustment: 
#   multixts <- highfrequency:::.multixts(rdata);
#   
#   if(multixts){
#     
#     result <- apply.daily(rdata, rDivVarFunBase, pow= pow, align.by= align.by, align.period = align.period, makeReturns = makeReturns, ...)
#     return(result)
#     
#   } else if(!multixts){
#     
#     result <- rDivVarFunBase(rdata, pow, align.by, align.period, makeReturns, ...)
#     return(result)
#   }
# }
# 
# rDivVarFunBase <- function(rdata, pow, align.by, align.period, makeReturns, ...){
#  
#   if((!is.null(align.by))&&(!is.null(align.period))){
#     rdata <- aggregatePrice(rdata, on=align.by, k=align.period, ...);
#   }
#   
#   if(makeReturns){
#     rdata <- makeReturns(rdata) 
#   }  
#   
#   result <- apply(matrix(pow),1,divVarFoo, tsMat = rdata, align.by = align.by, align.period = align.period, makeReturns = FALSE) 
#   
#   return(result)  
# }
# 
# divVarFoo <- function(p, tsMat, ...){
#   
#   # spot var
#   spot.var <- spotVol(rdata = tsMat, spot.index = index(tsMat), estFun = spotVolBaseJump, ...)
# 
#   # divergence derivative squared
#   div.deriv.sq <- apply(X = matrix(p), MARGIN = 1, FUN = divDerivFoo, tsMat = tsMat, .sum = FALSE)^2
#   
#   # multiply, sum, look what gets
#   spot.var <- matrix(as.numeric(spot.var), nrow = nrow(spot.var), ncol = ncol(div.deriv.sq))
#   result <- apply(spot.var*div.deriv.sq, 2, sum)
#   
#   # the estimator doesn't exactly yield what we need: correct with quarticity
#   RQ <- as.numeric(rMPV(rdata = tsMat, r.tot = 4, num.int = 10, makeReturns = FALSE, align.by = NULL, align.period = NULL))
#   result <- max(result - 0.5*RQ, 0.5*RQ)
#   
#   return(result)
# }