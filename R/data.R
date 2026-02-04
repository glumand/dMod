#' Compare data and model prediction by computing residuals
#'
#' @param data data.frame with columns: time, name, value, sigma, lloq
#' @param out matrix with predictions (first column = time), optional "deriv" attribute
#' @param err optional matrix with error model predictions, optional "deriv" attribute
#' @return objframe with attributes "deriv", "deriv.err"
#' @export
res <- function(data, out, err = NULL) {
  
  data$name <- as.character(data$name)
  n <- nrow(data)
  
  times <- sort(unique(data$time))
  names <- unique(data$name)
  
  dt <- match.num(data$time, times)
  dn <- match(data$name, names)
  
  ot <- match.num(times, out[,1])
  on <- match(names, colnames(out))
  if (anyNA(on))
    stop("Observable not found: ",
         paste(setdiff(names, colnames(out)), collapse = ", "))
  
  ti <- ot[dt]; ni <- on[dn]
  if (anyNA(ti)) stop("Some data$time not found in out[,1]")
  
  pred <- out[cbind(ti, ni)]
  
  deriv <- NULL
  if (!is.null(d <- attr(out, "deriv"))) {
    oi <- match(data$name, dimnames(d)[[2]])
    np <- dim(d)[3]
    deriv <- matrix(
      d[cbind(rep(ti, np), rep(oi, np), rep(seq_len(np), each = n))],
      n, np, dimnames = list(NULL, dimnames(d)[[3]])
    )
  }
  
  sig <- data$sigma
  sNA <- is.na(sig)
  derr <- NULL
  
  if (any(sNA)) {
    if (is.null(err)) stop("NA sigmas but no errmodel")
    et <- match.num(times, err[,1])
    en <- match(names, colnames(err))
    ti_e <- et[dt]; ni_e <- en[dn]
    sig[sNA] <- err[cbind(ti_e, ni_e)][sNA]
    
    if (!is.null(de <- attr(err, "deriv"))) {
      oi <- match(data$name, dimnames(de)[[2]])
      np <- dim(de)[3]
      derr <- matrix(0, n, np, dimnames = list(NULL, dimnames(de)[[3]]))
      derr[sNA,] <- matrix(
        de[cbind(rep(ti_e[sNA], np),
                 rep(oi[sNA], np),
                 rep(seq_len(np), each = sum(sNA)))],
        sum(sNA), np
      )
    }
  }
  
  val <- pmax(data$value, data$lloq)
  resi <- pred - val
  inv <- 1 / sig
  
  objframe(
    data.frame(data[c("time","name")], value = val,
               prediction = pred, sigma = sig,
               residual = resi, weighted.residual = resi * inv,
               bloq = val <= data$lloq, weighted.0 = pred * inv),
    deriv = deriv, deriv.err = derr
  )
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


