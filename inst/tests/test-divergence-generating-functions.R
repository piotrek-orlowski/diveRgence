
####
# Tests of correctness of unweighted divergence generating functions, and their derivatives
####
library(diveRgence)
library(testthat)
library(numDeriv)

p.range <- seq(-0.5,1,by=0.1)

foo.list <- c("divFoo","divSkewFoo","divKurtFoo","uDivFoo","divUSkewFoo","divUKurtFoo")
foo.list <- paste0("diveRgence:::",foo.list)

#### Is uDivFoo well defined ? ####
ts <- c(0.2,-0.22, -0.3)

divs <- sapply(p.range, diveRgence:::divFoo, tsMat = ts, .sum = FALSE)

divs.u.check <- sapply(p.range, function(pp) exp(pp * head(cumsum(c(0,ts)),-1)))
divs.u.check <- divs.u.check * divs
divs.u <- sapply(p.range, diveRgence:::uDivFoo, tsMat = ts, .sum = FALSE)


#### are rUSkewFoo and rUQuartFoo well defined? ####
skew.check <- sapply(p.range, function(pp) {
      jacobian(diveRgence:::uDivFoo, x = pp, tsMat = ts, .sum = FALSE, method = "Richardson", method.args = list(eps=1e-3, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))}
   )
skew <- sapply(p.range, diveRgence:::divUSkewFoo, tsMat = ts, .sum = FALSE)

quart.check <- sapply(p.range, function(pp) {
  jacobian(diveRgence:::divUSkewFoo, x = pp, tsMat = ts, .sum = FALSE, method = "complex")}
)

quart <- sapply(p.range, diveRgence:::divUKurtFoo, tsMat = ts, .sum = FALSE)

#### is the Z-deriv of uDivFoo correct? ####
divs.u.z.deriv <- sapply(p.range, FUN = diveRgence:::rUDivBaseZDeriv, tsMat = ts, .sum= F)

divs.u.z.deriv.check <- sapply(p.range, FUN = function(pp, tsMat, .sum){
    res <- diveRgence:::divFoo(p = pp, tsMat = tsMat, .sum = .sum)
    res <- pp * res * exp(pp * head(cumsum(c(0,tsMat)),-1))
    return(res)
}, tsMat= ts, .sum = FALSE)

#### is the Y-deriv of uDivFoo correct? ####
divs.u.y.deriv <- sapply(p.range, FUN = diveRgence:::rUDivBaseDeriv, tsMat = ts, .sum= F)

divs.u.y.deriv.check <- sapply(p.range, FUN = function(pp, tsMat, .sum){
  res <- diveRgence:::rDivBaseDeriv(p = pp, tsMat = tsMat, .sum = .sum)
  res <- res * exp(pp * head(cumsum(c(0,tsMat)),-1))
  return(res)
}, tsMat= ts, .sum = FALSE)

#### how about the Z-derivs of skew and quart ####
skew.z.deriv <- sapply(p.range, FUN = diveRgence:::rUSkewBaseZDeriv, tsMat = ts, .sum = F)

skew.z.deriv.check <- sapply(p.range, function(pp) {
  jacobian(diveRgence:::rUDivBaseZDeriv, x = pp, tsMat = ts, .sum = FALSE, method = "complex", method.args = list(eps=1e-3, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))}
)

quart.z.deriv <- sapply(p.range, FUN = diveRgence:::rUKurtBaseZDeriv, tsMat = ts, .sum = F)

quart.z.deriv.check <- sapply(p.range, function(pp) {
  jacobian(diveRgence:::rUSkewBaseZDeriv, x = pp, tsMat = ts, .sum = FALSE, method = "complex", method.args = list(eps=1e-3, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))}
)

#### how about the Y-derivs of skew and quart ####
skew.y.deriv <- sapply(p.range, FUN = diveRgence:::rUSkewBaseDeriv, tsMat = ts, .sum = F)

skew.y.deriv.check <- sapply(p.range, function(pp) {
  jacobian(diveRgence:::rUDivBaseDeriv, x = pp, tsMat = ts, .sum = FALSE, method = "complex", method.args = list(eps=1e-3, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))}
)

quart.y.deriv <- sapply(p.range, FUN = diveRgence:::rUKurtBaseDeriv, tsMat = ts, .sum = F)

quart.y.deriv.check <- sapply(p.range, function(pp) {
  jacobian(diveRgence:::rUSkewBaseDeriv, x = pp, tsMat = ts, .sum = FALSE, method = "complex", method.args = list(eps=1e-3, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))}
)
