## ---------------------------------------------------------------------------
## Shared setup for the scripts in this folder: loads the installed dMod2 and
## defines the local helpers that the old in-repo dMod fork used to provide
## (mstrust_fim()). Sourced by every script via
##   source(file.path(script_dir, "setup.R"))
## Fixes needed on top of dMod2 live here, never in dMod2 itself.
## ---------------------------------------------------------------------------

# The old dMod is still installed for other projects, but it must never be
# loaded next to dMod2: both register S3 methods for the same classes
# ("[.parvec", "*.fn", ...), so whichever loads last silently replaces the
# other's. Unloading doesn't help: R's S3 table keeps pointing into the old
# namespace and reloads it as soon as dMod2 registers its methods. So refuse
# to start in a session that has seen dMod, and warn if anything (library(),
# dMod::f, dMod:::f, a package importing it) loads it later. The hook can
# only warn: R wraps load hooks in try(), which swallows errors.
if ("dMod" %in% loadedNamespaces())
  stop("The old dMod is loaded in this R session; restart R before running ",
       "this script (it clashes with dMod2).", call. = FALSE)
setHook(packageEvent("dMod", "onLoad"), function(...)
  warning("The old dMod was just loaded -- results from here on may be wrong ",
          "(it clashes with dMod2). Find the dMod:: call and restart R.",
          call. = FALSE, immediate. = TRUE), action = "replace")

suppressMessages(library(dMod2))
stopifnot(packageVersion("dMod2") >= "0.12.7")
cat("Using dMod2", as.character(packageVersion("dMod2")),
    "from", find.package("dMod2"), "\n\n")

# Indices where a sorted objective-value vector jumps by more than tol (the
# start of each plateau in a waterfall). Copied from the old dMod fork, where
# it was the internal dMod:::stepDetect(); dMod2 no longer has it.
stepDetect <- function(x, tol) {
  jumps <- 1
  while (TRUE) {
    i <- which(x - x[1] > tol)[1]
    if (is.na(i)) break
    jumps <- c(jumps, tail(jumps, 1) - 1 + i)
    x <- x[-seq(1, i - 1, 1)]
  }
  jumps
}

# Proxies objfun and, on its FIRST call only, prints the singular values and
# condition number of the flattened [time*observable, parameter] sensitivity
# matrix per condition. Never evaluates objfun an extra time, so
# call-count-sensitive objectives are unaffected.
.mstrust_sensitivity_wrapper <- function(objfun, i) {
  printed <- FALSE
  wrapped <- function(...) {
    out <- objfun(...)
    if (!printed) {
      printed <<- TRUE
      env        <- attr(out, "env")
      prediction <- if (!is.null(env)) env$prediction else NULL
      cat("\n===== Fit", i, ": sensitivity singular values at initial evaluation =====\n")
      if (is.null(prediction)) {
        cat("(not available -- objfun did not attach a prediction with sensitivities)\n")
      } else {
        for (cond in names(prediction)) {
          S <- attr(prediction[[cond]], "deriv")
          cat("-- condition:", cond, "--\n")
          if (is.null(S)) {
            cat("(no sensitivities attached for this condition)\n")
          } else {
            J <- S
            dim(J) <- c(dim(S)[1] * dim(S)[2], dim(S)[3])
            colnames(J) <- dimnames(S)[[3]]
            sv <- sort(svd(J, nu = 0, nv = 0)$d, decreasing = TRUE)
            cat("singular values:", paste(signif(sv, 6), collapse = ", "), "\n")
            cat("condition number (max/min):", signif(max(sv) / min(sv), 6), "\n")
          }
        }
      }
    }
    out
  }
  # mstrust() dispatches on the objective's class/attributes (loadDLL(),
  # resetWarmStarts(), ...), so the proxy has to carry them.
  attributes(wrapped) <- attributes(objfun)
  wrapped
}

#' mstrust() with per-fit FIM diagnostics: identical to dMod2::mstrust(),
#' except that each fit's first objfun evaluation prints the singular values
#' of the sensitivity matrix (see .mstrust_sensitivity_wrapper()). Implemented
#' as an optmethod that wraps objfun freshly for every fit (and retry), then
#' hands over to the real optimiser.
mstrust_fim <- function(objfun, center, ..., optmethod = "trust") {
  optfun <- if (is.function(optmethod)) optmethod else match.fun(optmethod)
  fit_counter <- 0L
  fim_opt <- function() {
    args <- mget(setdiff(names(formals(sys.function())), "..."), environment())
    fit_counter <<- fit_counter + 1L
    args$objfun <- .mstrust_sensitivity_wrapper(args$objfun, fit_counter)
    do.call(optfun, c(args, list(...)))
  }
  # Same formals as the optimiser, so mstrust() routes trust() arguments to it.
  formals(fim_opt) <- formals(optfun)
  mstrust(objfun, center, ..., optmethod = fim_opt)
}
