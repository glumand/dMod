#' @export
#' @rdname trust
#' @param mu Named numeric vector of reference values for the L1-penalised
#'   parameters. Names must be a subset of \code{names(parinit)}; only the
#'   named parameters receive a penalty. Defaults to a zero vector covering
#'   all of \code{parinit}.
#' @param one.sided Logical. If \code{TRUE}, the penalty is one-sided and
#'   acts as a lower wall at \code{mu}: \code{lambda * max(0, mu - p)}.
#'   Otherwise it is the symmetric \code{lambda * |p - mu|}.
#' @param lambda Strength of the L1 penalty. Either a scalar (broadcast to
#'   all entries of \code{mu}) or a named numeric aligned with \code{mu}.
trustL1 <- function(objfun, parinit, mu = 0 * parinit, one.sided = FALSE, lambda = 1,
                    rinit, rmax,
                    parscale  = NULL,
                    iterlim   = 100L,
                    fterm     = 1e-6,
                    mterm     = 1e-6,
                    minimize  = TRUE,
                    blather   = FALSE,
                    parupper  = NULL,
                    parlower  = NULL,
                    printIter = FALSE,
                    traceFile = NULL,
                    ...) {
  sanePars <- sanitizePars(parinit, list(...)$fixed)
  parinit  <- sanePars$pars

  if (is.null(names(parinit)))
    stop("trustL1: parinit must be a named numeric vector")
  if (length(mu) > 0L && is.null(names(mu)))
    stop("trustL1: mu must be a named numeric vector")

  unknown <- setdiff(names(mu), names(parinit))
  if (length(unknown) > 0L)
    stop("trustL1: mu has names not present in parinit: ",
         paste(unknown, collapse = ", "))

  if (length(lambda) == 1L) {
    lambda <- structure(rep(as.numeric(lambda), length(mu)),
                        names = names(mu))
  } else {
    if (is.null(names(lambda)))
      stop("trustL1: lambda must be scalar or a named numeric vector")
    if (!setequal(names(lambda), names(mu)))
      stop("trustL1: names(lambda) must equal names(mu)")
    lambda <- lambda[names(mu)]
  }

  dots <- list(...)
  fn <- if (length(dots) > 0L) {
    function(x) do.call(objfun, c(list(x), dots))
  } else {
    objfun
  }

  mu <- structure(as.numeric(mu), names = names(mu))
  lambda <- structure(as.numeric(lambda), names = names(lambda))

  trustL1_impl(fn, parinit, mu, lambda,
               as.logical(one.sided)[1L], rinit, rmax,
               parscale, as.integer(iterlim), fterm, mterm,
               minimize, blather, parupper, parlower, printIter, traceFile)
}
