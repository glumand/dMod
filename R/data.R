#' Compute residuals between data and model prediction
#'
#' Matches data to predictions by time and observable, computes (weighted)
#' residuals, and propagates parameter derivatives. Values below `lloq` are
#' censored via `pmax(value, lloq)`.
#'
#' @md
#' @param data Data frame with columns `time`, `name`, `value`, `sigma`, `lloq`.
#'   Rows with `sigma = NA` are filled from `err`.
#' @param out Prediction matrix (first column = time, remaining = observables).
#'   Optional `"deriv"` attribute: `[name, param, time]` array.
#' @param err Optional error-model matrix (same layout as `out`).
#'   Optional `"deriv"` attribute: `[name, param, time]` array.
#'
#' @details
#' The returned `"deriv"` and `"deriv.err"` matrices have shape
#' \eqn{n \times p}{n x p} (residuals x parameters), extracted from the
#' `[name, param, time]` arrays on `out` and `err`.
#'
#' @return An [objframe()] with columns `time`, `name`, `value`, `prediction`,
#'   `sigma`, `residual`, `weighted.residual`, `bloq`, `weighted.0` and
#'   attributes `"deriv"` and `"deriv.err"`.
#'
#' @seealso [objframe()]
#' @export
res <- function(data, out, err = NULL) {
  
  data$name <- as.character(data$name)
  n <- nrow(data)
  times <- sort(unique(data$time))
  names <- unique(data$name)
  
  ti <- match.num(times, out[, 1])[match.num(data$time, times)]
  ni <- match(names, colnames(out))[match(data$name, names)]
  if (anyNA(ni))
    stop("Observable not found: ",
         paste(setdiff(names, colnames(out)), collapse = ", "))
  if (anyNA(ti)) stop("Some data$time not found in out[,1]")
  
  pred <- out[cbind(ti, ni)]
  
  deriv <- NULL
  if (!is.null(d <- attr(out, "deriv"))) {
    oi <- match(data$name, dimnames(d)[[1]])
    np <- dim(d)[2]
    deriv <- matrix(
      d[cbind(rep(oi, np), rep(seq_len(np), each = n), rep(ti, np))],
      n, np, dimnames = list(NULL, dimnames(d)[[2]]))
  }
  
  sig  <- data$sigma
  sNA  <- is.na(sig)
  derr <- NULL
  
  if (any(sNA)) {
    if (is.null(err)) stop("NA sigmas but no errmodel")
    ti_e <- match.num(times, err[, 1])[match.num(data$time, times)]
    ni_e <- match(names, colnames(err))[match(data$name, names)]
    sig[sNA] <- err[cbind(ti_e, ni_e)][sNA]
    
    if (!is.null(de <- attr(err, "deriv"))) {
      oi <- match(data$name, dimnames(de)[[1]])
      np <- dim(de)[2]
      ns <- sum(sNA)
      derr <- matrix(0, n, np, dimnames = list(NULL, dimnames(de)[[2]]))
      derr[sNA, ] <- matrix(
        de[cbind(rep(oi[sNA], np), rep(seq_len(np), each = ns), rep(ti_e[sNA], np))],
        ns, np)
    }
  }
  
  val  <- pmax(data$value, data$lloq)
  resi <- pred - val
  inv  <- 1 / sig
  
  objframe(
    data.table::data.table(
      time = data$time, name = data$name, value = val,
      prediction = pred, sigma = sig, residual = resi,
      weighted.residual = resi * inv,
      bloq = val <= data$lloq, weighted.0 = pred * inv),
    deriv = deriv, deriv.err = derr)
}



#' Time-course data for the JAK-STAT cell signaling pathway
#'
#' Phosphorylated Epo receptor (pEpoR), phosphorylated STAT in the
#' cytoplasm (tpSTAT) and total STAT (tSTAT) in the cytoplasmhave been 
#' measured at times 0, ..., 60.
#'
#' @name jakstat
#' @docType data
#' @keywords data
NULL


#' Time-course data for the Bile-Acid demonstration model
#'
#' @name badata
#' @docType data
#' @keywords data
NULL


# Match with numeric tolerance 
match.num <- function(x, y, tol = 1e-8) {
  sapply(x, function(xi) {
    d <- abs(y - xi)
    if (min(d) > tol) return(NA_integer_)
    which.min(d)
  })
}


