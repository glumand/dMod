#' Batched Matrix Multiplication
#'
#' Efficient batched matrix multiplication with BLAS-backend.
#' Batch index is the FIRST dimension.
#'
#' @section Supported Contractions:
#' \describe{
#'   \item{\code{[B,M,K] x [K,N] -> [B,M,N]}}{Left batched}
#'   \item{\code{[M,K] x [B,K,N] -> [B,M,N]}}{Right batched}
#'   \item{\code{[B,M,K] x [B,K,N] -> [B,M,N]}}{Both batched}
#' }
#'
#' @param A Numeric matrix or 3D array
#' @param B Numeric matrix or 3D array
#' @return Numeric 3D array with dim \code{[B, M, N]}
#' @export
`%bmm%` <- function(A, B) {
  da <- dim(A)
  db <- dim(B)

  if (length(da) == 3 && length(db) == 2) {
    # [B,M,K] x [K,N] -> [B,M,N]
    stopifnot(da[3] == db[1])
    bmm_lb(A, B, da[1], da[2], da[3], db[2])

  } else if (length(da) == 2 && length(db) == 3) {
    # [M,K] x [B,K,N] -> [B,M,N]
    stopifnot(da[2] == db[2])
    bmm_rb(A, B, db[1], da[1], da[2], db[3])

  } else if (length(da) == 3 && length(db) == 3) {
    # [B,M,K] x [B,K,N] -> [B,M,N]
    stopifnot(da[1] == db[1], da[3] == db[2])
    bmm_bb(A, B, da[1], da[2], da[3], db[3])

  } else {
    stop("Invalid dimensions for %bmm%")
  }
}
