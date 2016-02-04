#' @title Add alpha
#' @description Take any colour and make it transparent by setting alpha < 1
#' @export addAlpha


addAlpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  return(apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  )
}