#' Batched Matrix Multiplication
#'
#' Efficient batched matrix multiplication with BLAS-backend.
#'
#' @section Supported Contractions:
#' \describe{
#'   \item{\code{[M,K,B] x [K,N] -> [M,N,B]}}{Left batched}
#'   \item{\code{[M,K] x [K,N,B] -> [M,N,B]}}{Right batched}
#'   \item{\code{[M,K,B] x [K,N,B] -> [M,N,B]}}{Both batched}
#' }
#'
#' @param A Numeric matrix or 3D array
#' @param B Numeric matrix or 3D array
#' @return Numeric 3D array with dim \code{[M, N, B]}
#' @export
`%bmm%` <- function(A, B) {
  da <- dim(A)
  db <- dim(B)
  
  if (length(da) == 3 && length(db) == 2) {
    # [M,K,B] x [K,N] -> [M,N,B]
    stopifnot(da[2] == db[1])
    bmm_lb(A, B, da[3], da[1], da[2], db[2])
    
  } else if (length(da) == 2 && length(db) == 3) {
    # [M,K] x [K,N,B] -> [M,N,B]
    stopifnot(da[2] == db[1])
    bmm_rb(A, B, db[3], da[1], da[2], db[2])
    
  } else if (length(da) == 3 && length(db) == 3) {
    # [M,K,B] x [K,N,B] -> [M,N,B]
    stopifnot(da[3] == db[3], da[2] == db[1])
    bmm_bb(A, B, da[3], da[1], da[2], db[2])
    
  } else {
    stop("Invalid dimensions for %bmm%")
  }
}