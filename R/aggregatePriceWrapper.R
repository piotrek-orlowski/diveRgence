#' @name aggregatePrice
#' @title Time-series price aggregation
#' @description This is a wrapper to \code{\link{highfrequency::aggregatePrice}} which checks for double time stamps in the results.
#' @details all argument and outputs as in \code{\link{highfrequency::aggregatePrice}}, except you can't specify the aggregation function because \code{\link{highfrequency}} has been coded up bad.
#' @export aggregatePrice

aggregatePrice <- function(ts, on = "minutes", k=1, marketopen = "07:30:00", marketclose = "14:15:00", tz = "GMT"){
  
  res <- highfrequency::aggregatePrice(ts, on = on, k=k, marketopen = marketopen, marketclose = marketclose, tz = tz)
  
  if(length(index(res)) == length(unique(index(res)))){
    return(res)
  } else {
    res <- lapply(unique(index(res)), function(xx){
      record <- res[xx]
      if(nrow(record) > 1){
        record <- xts(x = median(as.numeric(record)), order.by = unique(index(record)))
      }
      return(record)
    })
    res <- do.call(what = rbind.xts, args = res)
    return(res)
  }
}