#' Build per-subject quadrature nodes for adaptive NLME marginal integration
#'
#' @description
#' Given the posterior mode `etaHat` and negative-Hessian `Hi` of the joint
#' log-density at the mode (= inverse covariance of the local Gaussian
#' approximation), constructs sparse-grid Gauss-Hermite quadrature nodes in
#' eta-space along with the log-magnitudes of their augmented weights (with
#' signs tracked separately).
#'
#' The augmented weight `W_b` for the integral `int g(eta) deta` approximated
#' as `sum_b W_b * g(eta_b)` includes: the raw Smolyak-Gauss-Hermite weight
#' `w_b_GH`, the change-of-variable Jacobian `2^(K/2) / |det L_H|`, and the
#' `exp(z_b' z_b)` factor that un-does the implicit `exp(-z'z)` weight of the
#' physicists' GH rule. Combined: `log W_b = (K/2)*log(2) - log|det L_H| +
#' log|w_b_GH| + z_b' z_b`, with sign carried in `weightSigns`.
#'
#' @param etaHat Length-K numeric, posterior mode of the random effects.
#' @param Hi K x K positive-definite matrix, the joint's negative-Hessian
#'   at `etaHat` (= inverse covariance of the local Gaussian).
#' @param level Integer >= K, Smolyak depth. K+1 to K+3 is the useful range
#'   for typical NLME problems; higher levels exponentially raise node count.
#'
#' @return A list with components:
#' \describe{
#'   \item{`etaNodes`}{B x K matrix of quadrature nodes in eta-space
#'     (batch-first; row `b` is the b-th node).}
#'   \item{`logAbsWeights`}{Length-B numeric, `log|W_b|` of the augmented
#'     weights including the change-of-variable Jacobian. Use with
#'     `weightSigns` for downstream signed log-sum-exp.}
#'   \item{`weightSigns`}{Length-B numeric `{-1, +1}`, sign of `W_b`.}
#'   \item{`K`, `level`}{Echoed.}
#' }
#'
#' @seealso `sparse_grid_gh` (Rcpp-exported, in `src/quadrature.cpp`).
#' @keywords internal
makeSubjectNodes <- function(etaHat, Hi, level) {
  K <- length(etaHat)
  if (!is.matrix(Hi) || nrow(Hi) != K || ncol(Hi) != K)
    stop("`Hi` must be a K x K matrix matching length(etaHat).")
  if (!is.numeric(level) || length(level) != 1L || level < 1L)
    stop("`level` must be a positive integer scalar.")

  raw     <- sparse_grid_gh(K, as.integer(level))
  z_nodes <- raw$nodes
  w_raw   <- raw$weights
  B       <- nrow(z_nodes)

  R         <- chol(Hi)
  log_det_L <- sum(log(abs(diag(R))))

  z_solved  <- backsolve(R, t(z_nodes))
  etaNodes <- t(sqrt(2) * z_solved) +
               matrix(etaHat, B, K, byrow = TRUE)
  if (!is.null(names(etaHat))) colnames(etaNodes) <- names(etaHat)

  z2_sum   <- rowSums(z_nodes^2)
  log_abs  <- (K / 2) * log(2) - log_det_L + log(abs(w_raw)) + z2_sum
  signs    <- sign(w_raw)

  list(etaNodes      = etaNodes,
       logAbsWeights = log_abs,
       weightSigns   = signs,
       K             = K,
       level         = as.integer(level))
}
