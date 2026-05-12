

#' Named repititions
#' 
#' @description Wrapper on rep() to input names instead of length.
#'   
#' @param x Value to be repeated.
#' @param names List of names.
#'   
#' @export
repWithNames <- function(x, names){
  repnum <- rep(x,length(names))
  names(repnum) <- names
  return(repnum)
}

