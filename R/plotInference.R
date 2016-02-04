#' @name plotInference
#' @title Plot realized divergence inference results
#' @description This replaces the plot interface for objects of class \code{divinf}, output from the divergence estimating functions with asymptotic variance estimates.
#' @param x list with fields \code{rDiv} and \code{asy.var}.
#' @param true xts time series.
#' @param quantile for confidence bounds
#' @param ... graphical parameters
#' @export plot.divinf

plot.divinf <- function(x, true = NULL, quantile = 0.95, ...){
  
  estimated <- x
  
  est.ts <- xts(as.numeric(estimated$rDiv), order.by = as.Date(index(estimated$rDiv)))
  est.var <- xts(as.numeric(estimated$asy.var), order.by = as.Date(index(estimated$asy.var)))
  est.upper.pos <- est.ts + qnorm(p = 1-(1-quantile)/2)*sqrt(est.var)
  est.upper <- est.upper.pos
  est.upper.pos <- est.upper.pos[which(est.upper.pos >= 0)]
  est.lower.pos <- est.ts - qnorm(p = 1-(1-quantile)/2)*sqrt(est.var)
  est.lower <- est.lower.pos
  est.lower.pos <- est.lower.pos[which(est.lower.pos >= 0)]
  est.upper.neg <- est.ts + qnorm(p = 1-(1-quantile)/2)*sqrt(est.var)
  est.upper.neg <- est.upper.neg[which(est.upper.neg < 0)]
  est.lower.neg <- est.ts - qnorm(p = 1-(1-quantile)/2)*sqrt(est.var)
  est.lower.neg <- est.lower.neg[which(est.lower.neg < 0)]
  
  plot.range <- range(c(est.ts, est.upper.pos, est.upper.neg, est.lower.pos, est.lower.neg))
  
  if(!is.null(true)){
    true <- xts(as.numeric(true), order.by = as.Date(index(true)))
    plot.range <- range(c(plot.range, true))
  }
  
  plot(est.ts, type="h", ylim = plot.range, main="", xlab = "", ylab = "")
  par(new = TRUE)
  plot(index(est.upper.pos), as.numeric(est.upper.pos), col= "red", lty = 1, lwd = 2, ylim = plot.range, type="h", xlim = range(index(est.ts)), ylab = "Realized", xlab = "Time", main = "")
  par(new = TRUE)
  plot(index(est.lower.pos), as.numeric(est.lower.pos), col= "white", lty = 1, lwd = 2.5, ylim = plot.range, type="h", xlim = range(index(est.ts)), ylab = "", xlab = "", main = "")
  par(new = TRUE)
  plot(index(est.lower.neg), as.numeric(est.lower.neg), col= "red", lty = 1, lwd = 2, ylim = plot.range, type="h",xlim = range(index(est.ts)), main="", xlab = "", ylab = "")
  par(new = TRUE)
  plot(index(est.upper.neg), as.numeric(est.upper.neg), col= "white", lty = 1, lwd = 2.5, ylim = plot.range, type="h",xlim = range(index(est.ts)), main="", xlab = "", ylab = "")
  par(new = TRUE)
  plot(index(est.ts), as.numeric(est.ts), pch=20, lwd = 2, ylim = plot.range, col = addAlpha("yellow",0.8),xlim = range(index(est.ts)),main="", xlab = "", ylab = "")
  

  if(!is.null(true)){
    par(new=TRUE) 
    plot(index(true), as.numeric(true), pch = 20, col = addAlpha("purple",1), ylim = plot.range,xlim = range(index(est.ts)),main="", xlab = "", ylab = "")
  }
}