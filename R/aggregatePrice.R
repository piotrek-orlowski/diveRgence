#' @title Aggregate a time series to within open/close bounds
#' @description Function returns new time series as an \code{xts} object. The observations are always previous-tick observations. If \code{marketopen} is set to before the start of the time series \code{ts}, the timestamp is discarded and the first returned timestamp is the first in-fact time stamp. The \code{ts} object can contain multiple days of data.
#' @param ts \code{xts} object to be aggregated; \code{index(ts)} is \code{POSIXct}. It contains intraday observations of a security price, potentially for multiple days.
#' @param on character, indicating the time scale at which the series will be aggregated: \code{"seconds"}, \code{"minutes"} or \code{"hours"}.
#' @param k positive integer, aggregation is performed every k \code{on}s
#' @param marketopen string in format \code{"HH:MM:SS"}, start of trading time.
#' @param marketclose string in format \code{"HH:MM:SS"}, end of trading time.
#' @param pad.na logical, if \code{TRUE}, any time when the previous tick is earlier than \code{k} times \code{on}, the latter time stamp is populated with \code{NA_real}.
#' @export
aggregatePrice <- function(ts, on = 'minutes', k = 1, marketopen = "08:30:00", marketclose = "15:15:00", pad = TRUE, pad.arg = NA_real_, aggr.vec = numeric(0)){
  
  marketopen <- as.POSIXlt(x = marketopen, tz = "UTC", format = "%H:%M:%S")
  marketclose <- as.POSIXlt(x = marketclose, tz = "UTC", format = "%H:%M:%S")
  marketopen <- c(marketopen$hour, marketopen$min, marketopen$sec)
  marketclose <- c(marketclose$hour, marketclose$min, marketclose$sec)
  
  rdata <- aggregatePrice_Xts(rdata = ts, period_ = on, numPeriods_ = k, dayStart_ = marketopen, dayEnd_ = marketclose, aggr_vec = as.numeric(aggr.vec), pad = pad, pad_arg = pad.arg)
  
  return(rdata)
}