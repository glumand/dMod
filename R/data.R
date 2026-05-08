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
    oi <- match(data$name, dimnames(d)[[2]])
    np <- dim(d)[3]
    deriv <- matrix(
      d[cbind(rep(ti, np), rep(oi, np), rep(seq_len(np), each = n))],
      n, np, dimnames = list(NULL, dimnames(d)[[3]]))
  }

  deriv2 <- NULL
  if (!is.null(d2 <- attr(out, "deriv2"))) {
    oi2 <- match(data$name, dimnames(d2)[[2]])
    np2 <- dim(d2)[3]
    # Build [n*np*np x 4] index matrix; outermost loop = k, then j, then i.
    idx <- cbind(
      rep(ti,  np2 * np2),
      rep(oi2, np2 * np2),
      rep(rep(seq_len(np2), each = n), np2),
      rep(seq_len(np2), each = n * np2)
    )
    deriv2 <- array(d2[idx], c(n, np2, np2),
                    dimnames = list(NULL, dimnames(d2)[[3]], dimnames(d2)[[4]]))
  }

  sig  <- data$sigma
  sNA  <- is.na(sig)
  derr <- NULL
  derr2 <- NULL

  if (any(sNA)) {
    if (is.null(err)) stop("NA sigmas but no errmodel")
    ti_e <- match.num(times, err[, 1])[match.num(data$time, times)]
    ni_e <- match(names, colnames(err))[match(data$name, names)]
    sig[sNA] <- err[cbind(ti_e, ni_e)][sNA]

    if (!is.null(de <- attr(err, "deriv"))) {
      oi <- match(data$name, dimnames(de)[[2]])
      np <- dim(de)[3]
      ns <- sum(sNA)
      derr <- matrix(0, n, np, dimnames = list(NULL, dimnames(de)[[3]]))
      derr[sNA, ] <- matrix(
        de[cbind(rep(ti_e[sNA], np), rep(oi[sNA], np), rep(seq_len(np), each = ns))],
        ns, np)
    }

    if (!is.null(de2 <- attr(err, "deriv2"))) {
      oi <- match(data$name, dimnames(de2)[[2]])
      np2 <- dim(de2)[3]
      ns <- sum(sNA)
      derr2 <- array(0, c(n, np2, np2),
                     dimnames = list(NULL, dimnames(de2)[[3]], dimnames(de2)[[4]]))
      idx <- cbind(
        rep(ti_e[sNA],  np2 * np2),
        rep(oi[sNA],    np2 * np2),
        rep(rep(seq_len(np2), each = ns), np2),
        rep(seq_len(np2), each = ns * np2)
      )
      derr2[sNA, , ] <- array(de2[idx], c(ns, np2, np2))
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
    deriv = deriv, deriv.err = derr,
    deriv2 = deriv2, deriv2.err = derr2)
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


