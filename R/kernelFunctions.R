#' @name deltaKernels
#' @title Kernel functions for delta sequences
#' @description Functions for constructing delta-sequences, based on kernels, for spot estimation from HF data.
#' @param x0 kernel evaluation argument
#' @param cnv.data data for kernel-based convolution
#' @param delta.bandwidth kernel bandwiwdth
#' @param normalize if \code{TRUE}, the resulting weights will be normalized so that \code{sum(expDoubleKernel(x0,cnv.data,0.01)*0.001)}
#' @export expDoubleKernel

expDoubleKernel <- function(x0, cnv.data, delta.bandwidth){
  
  cnv.data.0 <- cnv.data
  cnv.data <- cnv.data - x0
  cnv.data <- as.numeric(abs(cnv.data))
  result <- 1/delta.bandwidth * exp(-2*cnv.data/delta.bandwidth)
  
  result <- xts(x = result, order.by = cnv.data.0)
  
  return(result)
}

#' @export gaussianKernel
#' 

gaussianKernel <- function(x0, cnv.data, delta.bandwidth){
  
  cnv.data <- cnv.data - x0
  renorm <- mean(diff(cnv.data))
  result <- 1/delta.bandwidth * dnorm(x = cnv.data/delta.bandwidth)
  #result <- result/sum(result * renorm)
  
  return(result)
  
}