% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/spotVol.R
\name{spotVolatility}
\alias{spotVol}
\alias{spotVolatility}
\title{Spot volatility estimation}
\usage{
spotVol(rdata, spot.index = NULL, makeReturns, align.by, align.period,
  apply.at, apply.period, avg.vol, reference.time = "07:30:00",
  vol.jumping = TRUE, year.days = 365, seconds.per.day = 86400,
  kernel.type = "gaussian", ...)
}
\arguments{
\item{rdata}{an \code{xts} object containing 1 return series.}

\item{spot.index}{POSIXct vector of times when volatility should be estimated. If this argument is given, \code{apply.at} and \code{apply.period} are ignored.}

\item{makeReturns}{boolean, should be \code{TRUE} when price data is supplied. Defaults to \code{FALSE}.}

\item{align.by}{}

\item{align.period}{}

\item{apply.at}{at what time unit intervals should the estimator be applied?}

\item{apply.period}{every how many time unit intervals should the estimator be applied?}

\item{...}{Arguments passed on to \code{\link[highfrequency]{aggregatePrice}}}

\item{estFun}{local volatility estimating function, one of \code{spotVolBase*} specified in this document.}

\item{nbar}{two-scale (two-frequency) parameter for eliminating microstructure noise in estimation.}
}
\description{
Function for calculating spot volatility under the assumptions of (a) continuous underlying process, (b) noisy observations of a continuous underlying process, (c) underlying process with jumps. The \code{apply.at} and \code{apply.period} parameters indicate at which time stamps the estimator should be applied. If you specify to \code{align.by='seconds'}, \code{align.period=5}, then pick, for example, \code{apply.at='minutes'}, \code{apply.by=10}. For more details, see Mancini, Matiussi and Reno (2012).
}
\details{
The most important arguments to pass go \code{\link[highfrequency]{aggregatePrice}} are \code{marketopen} and \code{marketclose}, see documentation therein. The default values are different from our test data set.
}

