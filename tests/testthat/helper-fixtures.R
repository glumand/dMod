# Compiled-model fixtures for behavioral tests.
#
# These functions build small dMod prediction / observation / trafo chains
# once per R process and cache the result via a globalenv-based store
# (testthat sources helper-*.R fresh per test_file, so a function-local
# cache would not survive across files; globalenv does).
#
# Convention: every compiled artifact lives in <tempdir>/dmod_fx and is
# named with a stable prefix so dyn.load can reuse it across calls.


## ---- Cache + workdir machinery ------------------------------------------

.dmod_fx_cache <- function() {
  if (!exists("..dmod_fx_cache..", envir = globalenv(), inherits = FALSE))
    assign("..dmod_fx_cache..", new.env(parent = emptyenv()),
           envir = globalenv())
  get("..dmod_fx_cache..", envir = globalenv())
}

.dmod_fx_workdir <- function() {
  cache <- .dmod_fx_cache()
  if (is.null(cache$workdir)) {
    cache$workdir <- file.path(tempdir(), "dmod_fx")
    dir.create(cache$workdir, recursive = TRUE, showWarnings = FALSE)
  }
  cache$workdir
}

# Internal: setwd to the fixture workdir for the duration of `expr`.
.dmod_with_fx_workdir <- function(expr) {
  oldwd <- setwd(.dmod_fx_workdir())
  on.exit(setwd(oldwd), add = TRUE)
  force(expr)
}


## ---- Linear decay fixture -----------------------------------------------

# Build (and cache) the compiled chain for a one-state linear-decay model.
#   ODE:        d A / dt = -k * A
#   Observable: y = A
#   Trafo:      identity (A, k) or log (A = exp(A_log), k = exp(k_log))
#
# Returns a list with elements:
#   m           odemodel
#   xfn         prediction function from Xs(m)
#   gfn         observation function Y(y = A, ...)
#   pfn_id      identity parameter trafo, single condition C1
#   pfn_log     log-trafo, single condition C1
#   prd_id      gfn * xfn * pfn_id
#   prd_log     gfn * xfn * pfn_log
#   outerpars_id, outerpars_log:  default named parameter vectors
fx_decay_compiled <- function() {
  cache <- .dmod_fx_cache()
  if (!is.null(cache$decay)) return(cache$decay)

  .dmod_with_fx_workdir({
    reactions <- addReaction(eqnlist(), from = "A", to = "",
                             rate = "k*A",
                             description = "linear decay")
    # Compile via dMod's `compile()` instead of cOde's internal compileAndLoad
    # so the build flags include -w; otherwise R's default -Wall surfaces a
    # constant pile of unused-variable noise from cOde-generated code.
    m <- odemodel(reactions, modelname = "fx_decay", compile = FALSE)
    xfn <- Xs(m)

    gfn <- Y(c(y = "A"), f = xfn, condition = NULL, attach.input = FALSE,
             modelname = "fx_decay_obs", compile = FALSE)

    trafo_id  <- eqnvec(A = "A",         k = "k")
    trafo_log <- eqnvec(A = "exp(A_log)", k = "exp(k_log)")

    pfn_id <- P(trafo_id,  condition = "C1",
                modelname = "fx_decay_p_id",  compile = FALSE)
    pfn_log <- P(trafo_log, condition = "C1",
                 modelname = "fx_decay_p_log", compile = FALSE)

    compile(xfn, gfn, pfn_id, pfn_log, cores = 1)

    cache$decay <- list(
      m           = m,
      xfn         = xfn,
      gfn         = gfn,
      pfn_id      = pfn_id,
      pfn_log     = pfn_log,
      prd_id      = gfn * xfn * pfn_id,
      prd_log     = gfn * xfn * pfn_log,
      outerpars_id  = c(A = 1.0, k = 0.5),
      outerpars_log = c(A_log = 0, k_log = log(0.5))
    )
  })

  cache$decay
}

# Build a noisy single-condition decay dataset using closed-form A(t) plus
# Gaussian noise. Returns a datalist. Default times / sigma chosen so the
# fit is well-conditioned and the gradient at truth is small but not zero
# (good for numDeriv comparisons).
fx_decay_data <- function(pars  = c(A = 1.0, k = 0.5),
                          times = seq(0, 10, by = 1),
                          sigma = 0.05,
                          condition = "C1",
                          seed  = 1L) {
  df <- make_noisy_data(
    truth_fn  = function(t, p) p["A"] * exp(-p["k"] * t),
    pars      = pars,
    times     = times,
    name      = "y",
    sigma     = sigma,
    condition = condition,
    seed      = seed)
  as.datalist(df, split.by = "condition")
}

# Two-condition decay dataset: same model, two different "true" k values
# encoded in conditions C1 and C2. Used to validate multi-condition
# aggregation in normL2 / objfn composition.
fx_decay_data_multi <- function(parslist = list(C1 = c(A = 1.0, k = 0.5),
                                                C2 = c(A = 1.0, k = 1.0)),
                                times = seq(0, 10, by = 1),
                                sigma = 0.05,
                                seed = 1L) {
  out <- do.call(rbind, lapply(seq_along(parslist), function(i) {
    set.seed(seed + i - 1L)
    p <- parslist[[i]]
    cn <- names(parslist)[i]
    data.frame(name = "y", time = times,
               value = p["A"] * exp(-p["k"] * times) +
                 rnorm(length(times), 0, sigma),
               sigma = sigma, condition = cn,
               stringsAsFactors = FALSE)
  }))
  as.datalist(out, split.by = "condition")
}


## ---- BLOQ-augmented decay dataset ---------------------------------------

# Same decay simulation but with rows below a chosen LLOQ marked BLOQ.
# The LLOQ defaults to a value that censors ~30% of late-time observations
# (decay to small values), so all four BLOQ modes (M1 drop / M3 / M4NM /
# M4BEAL) have non-trivial work to do.
fx_decay_data_bloq <- function(pars  = c(A = 1.0, k = 0.5),
                               times = seq(0, 10, by = 0.5),
                               sigma = 0.05,
                               lloq  = 0.1,
                               condition = "C1",
                               seed  = 1L) {
  df <- make_noisy_data(
    truth_fn  = function(t, p) p["A"] * exp(-p["k"] * t),
    pars      = pars,
    times     = times,
    name      = "y",
    sigma     = sigma,
    condition = condition,
    lloq      = lloq,
    seed      = seed)
  as.datalist(df, split.by = "condition")
}


## ---- Backend parametrisation --------------------------------------------

# The objective functions now have a single C++ backend; these shims keep
# legacy call sites working by running each block once under that backend.
for_each_backend <- function(fn) {
  fn(TRUE)
}

with_cpp_backend <- function(enabled, code) {
  force(code)
}
