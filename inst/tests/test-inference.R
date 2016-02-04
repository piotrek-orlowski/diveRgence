### TEST FILE FOR PACKAGE

# This tests whether the functions work, not whether they're correct!

### ---- SETUP ----
library(highfrequency)
library(affinesim)

data(sp500)

# To aggregate data it's not enough to use aggregatePrice, because (a) it has some unwanted behaviour, as filling in the start and end of the day, if trading starts late and finishes early, and (b) it sometimes churns out two records per timestamp, which is totally not nice.

dates <- unique(as.Date(index(sp500)))
sp500.list <- lapply(dates, function(dd){
  print(paste0("Aggregating SP500 data for ", dd))
  loc.sp <- sp500[as.character(dd)]
  loc.hrspan <- range(index(loc.sp))
  if(as.POSIXlt(loc.hrspan[1])$hour < 7){
    loc.index <- as.POSIXlt(index(loc.sp))
    loc.index$hour <- loc.index$hour + rep(1, length(loc.index$hour))
    index(loc.sp) <- as.POSIXct(loc.index)
  }
  loc.hrspan <- range(index(loc.sp))
  loc.hrspan <- as.character(loc.hrspan, format = "%H:%M:%S")
  substr(loc.hrspan[1],start = 7,stop = 8) <- "00"
  substr(loc.hrspan[2],start = 7,stop = 8) <- "00"
  substr(loc.hrspan[2],start = 4,stop = 5) <- sprintf("%1.2d",as.numeric(substr(loc.hrspan[2],start = 4,stop = 5)) + (60 - as.numeric(substr(loc.hrspan[2],start = 4,stop = 5))) %% 15)
  agg <- aggregatePrice(ts = loc.sp, on = "minutes", k = 1, marketopen = loc.hrspan[1], marketclose = loc.hrspan[2])
  
  return(agg)
})
sp500 <- do.call(rbind.xts, sp500.list)

plot(sp500)

### --- Define a plottin function ----

plotDivInf <- function(rdiv, title = ""){
  plot(rdiv$rDiv, ylim = range(rdiv$rDiv[,c(1,1)] + rdiv$asy.var), main = title)
  lines(index(rdiv$rDiv), as.numeric(rdiv$rDiv + rdiv$asy.var[,1]), col = "red", lty = 2)
  lines(index(rdiv$rDiv), as.numeric(rdiv$rDiv + rdiv$asy.var[,2]), col = "red", lty = 2)  
}

### ---- Realized Divergence----
print("Testing realized divergences")
rDiv.sp500 <- rDivEngine(rdata = sp500, fooStr = "rDiv", pow = c(1/2), makeReturns = T, align.by = NULL, align.period = NULL)

plot(rDiv.sp500, type = "h", title = "RDiv of return")
system.time(
rDiv.sp500.inf <- rDivEngineInference(rdata = sp500["1996-07-01/1996-09-30"], fooStr = "rDiv", pow = c(1/2), align.by = NULL, align.period = NULL, makeReturns = T, reference.time = "07:30:00")
)

plotDivInf(rDiv.sp500.inf, title = "RDiv of return")

rUDiv.sp500 <- rDivEngine(rdata = sp500, fooStr = "rUDiv", pow = c(1/2,0), makeReturns = T, align.by = NULL, align.period = NULL)

rUDiv.sp500.inf <- rDivEngineInference(rdata = sp500["1996-07-01/1996-09-30"], fooStr = "rUDiv", pow = c(1/2), align.by = NULL, align.period = NULL, makeReturns = T, reference.time = "07:30:00")

plotDivInf(rUDiv.sp500.inf, title = "rUDiv")


### ---- Realized Skewness ----

rJSkew.sp500 <- rDivEngine(rdata = sp500, fooStr = "rJSkew", pow = c(1/2), makeReturns = T, align.by = NULL, align.period = NULL)

rJSkew.sp500.inf <- rDivEngineInference(rdata = sp500["1996-07-01/1996-09-30"], fooStr = "rJSkew", pow = c(1/2), align.by = NULL, align.period = NULL, makeReturns = T, reference.time = "07:30:00")

plotDivInf(rJSkew.sp500.inf, title = "rJSkew")

rUSkew.sp500 <- rDivEngine(rdata = sp500, fooStr = "rUSkew", pow = c(1/2), makeReturns = T, align.by = NULL, align.period = NULL)

rUSkew.sp500.inf <- rDivEngineInference(rdata = sp500["1996-07-01/1996-09-30"], fooStr = "rUSkew", pow = c(1/2), align.by = NULL, align.period = NULL, makeReturns = T, reference.time = "07:30:00")

plotDivInf(rUSkew.sp500.inf, title = "rUSkew")

### ---- Realized kurtosis ----

rJKurt.sp500 <- rDivEngine(rdata = sp500, fooStr = "rJKurt", pow = c(1/2), makeReturns = T, align.by = NULL, align.period = NULL)

rJKurt.sp500.inf <- rDivEngineInference(rdata = sp500["1996-07-01/1996-09-30"], fooStr = "rJKurt", pow = c(1/2), align.by = NULL, align.period = NULL, makeReturns = T, reference.time = "07:30:00")

plotDivInf(rJKurt.sp500.inf, title = "rJKurt")

rUKurt.sp500 <- rDivEngine(rdata = sp500, fooStr = "rUKurt", pow = c(1/2), makeReturns = T, align.by = NULL, align.period = NULL)

rUKurt.sp500.inf <- rDivEngineInference(rdata = sp500["1996-07-01/1996-09-30"], fooStr = "rUKurt", pow = c(1/2), align.by = NULL, align.period = NULL, makeReturns = T, reference.time = "07:30:00")

plotDivInf(rUKurt.sp500.inf, title = "rUKurt")
