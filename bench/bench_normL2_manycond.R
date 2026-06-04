## Scaling sweep: how does the R-glue share in normL2 evolve as the
## condition count grows? Each condition adds:
##   - one res() call (data.table indexing, attribute slicing)
##   - one objlist accumulation in the Reduce(`+`, ...)
##   - one prediction-function dispatch in g*x*p
##
## We replicate the BA "closed" data N times under fresh condition names,
## with one shared trafo (per-condition keys via branch()), measure
## prd, normL2, glue share, and Rprof a 20-restart trust loop.

suppressPackageStartupMessages({
  library(dMod)
  library(dplyr)
  library(microbenchmark)
})

setwd(tempdir())
set.seed(11)

data(badata)
single <- subset(badata, condition == "closed")

reactions <- eqnlist() |>
  addReaction("TCA_buffer", "TCA_cell",  rate = "k_import*TCA_buffer") |>
  addReaction("TCA_cell",   "TCA_buffer", rate = "k_export_sinus*TCA_cell") |>
  addReaction("TCA_cell",   "TCA_cana",   rate = "k_export_cana*TCA_cell") |>
  addReaction("TCA_cana",   "TCA_buffer", rate = "k_reflux*TCA_cana")

mymodel <- odemodel(reactions, modelname = "ba_manycond",
                    solver = "CppODE", compile = FALSE)
x <- Xs(mymodel,
        optionsOde  = list(atol = 1e-8, rtol = 1e-8),
        optionsSens = list(atol = 1e-8, rtol = 1e-8))

observables <- eqnvec(buffer   = "s*TCA_buffer",
                      cellular = "s*(TCA_cana + TCA_cell)")
g <- Y(observables, f = x, condition = NULL,
       compile = FALSE, modelname = "obs_manycond", attach.input = FALSE)

innerpars <- getParameters(x, g)
base_trafo <- NULL |>
  define("x~x", x = innerpars) |>
  define("TCA_buffer~0") |>
  insert("x~exp10(y)", x = .currentSymbols, y = toupper(.currentSymbols))

make_problem <- function(N_cond) {
  cond_names <- sprintf("c%02d", seq_len(N_cond))
  df <- do.call(rbind, lapply(cond_names, function(cn) {
    d <- single; d$condition <- cn
    d$value <- d$value + rnorm(nrow(d), 0, d$sigma * 0.2)
    d
  }))
  data_N <- as.datalist(df)

  ## one P per condition, sharing the same outer parameter set:
  ## the trafo only differs by condition name in branch().
  trafo_branched <- branch(base_trafo, conditions = cond_names)
  p_N <- P(trafo_branched, compile = FALSE,
           modelname = paste0("p_branched_N", N_cond))
  list(data = data_N, p = p_N, conds = cond_names)
}

## build all sizes up-front (compile once)
sizes <- c(2, 8, 32, 64)
probs <- lapply(sizes, make_problem)
names(probs) <- as.character(sizes)

cat("Compiling g, x ... ")
t0 <- Sys.time()
compile(g, x, output = "ba_manycond_gx", cores = 4)
cat(sprintf("(%.1fs)\n", as.numeric(Sys.time() - t0, units = "secs")))
for (N in sizes) {
  cat(sprintf("Compiling P (N=%d) ... ", N))
  t1 <- Sys.time()
  compile(probs[[as.character(N)]]$p,
          output = paste0("ba_manycond_p_N", N), cores = 4)
  cat(sprintf("(%.1fs)\n", as.numeric(Sys.time() - t1, units = "secs")))
}
cat(sprintf("total compile: %.1fs\n\n",
            as.numeric(Sys.time() - t0, units = "secs")))

bench_one <- function(N_cond) {
  pp <- probs[[as.character(N_cond)]]
  prd <- g * x * pp$p
  outerpars <- getParameters(pp$p)
  pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
  obj <- normL2(pp$data, prd)
  times <- sort(unique(c(0, unlist(lapply(pp$data, `[[`, "time")))))
  invisible(prd(times, pouter, deriv = TRUE))
  invisible(obj(pouter))
  mb <- microbenchmark(
    prd_d1 = prd(times, pouter, deriv = TRUE),
    normL2 = obj(pouter),
    times = 20L, unit = "ms"
  )
  s <- summary(mb)
  list(N = N_cond,
       prd_ms    = s$median[s$expr == "prd_d1"],
       normL2_ms = s$median[s$expr == "normL2"])
}

cat("=== microbenchmark sweep ===\n")
res <- do.call(rbind, lapply(sizes, function(N) {
  r <- bench_one(N)
  data.frame(N_cond = r$N, prd_ms = r$prd_ms, normL2_ms = r$normL2_ms,
             glue_ms = r$normL2_ms - r$prd_ms,
             glue_pct = 100 * (r$normL2_ms - r$prd_ms) / r$normL2_ms)
}))
print(res, row.names = FALSE)

## N=64 detailed analysis: backend comparison + Rprof on the glue path
big <- probs[["64"]]
prd_big <- g * x * big$p
obj_big <- normL2(big$data, prd_big)
outerpars <- getParameters(big$p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
times_b <- sort(unique(c(0, unlist(lapply(big$data, `[[`, "time")))))
invisible(obj_big(pouter))

## cores scaling: the C++ residual kernel parallelises over conditions via
## OpenMP, controlled by normL2(..., cores = ).
cat("\n=== cores scaling (N = 64, deriv = TRUE) ===\n")
for (nc in c(1L, 2L, 4L)) {
  obj_nc <- normL2(big$data, prd_big, cores = nc)
  invisible(obj_nc(pouter))
  mb <- microbenchmark(obj_nc(pouter), times = 20L, unit = "ms")
  cat(sprintf("cores = %d: normL2 median = %.2f ms\n",
              nc, median(mb$time) / 1e6))
}

cat("\n=== Rprof: 200 obj_big() calls (N_cond = 64, default backend) ===\n")
pf <- tempfile(fileext = ".out")
Rprof(pf, interval = 0.005)
for (i in 1:200) invisible(obj_big(pouter))
Rprof(NULL)
prof <- summaryRprof(pf)$by.self
total <- sum(prof$self.time)
cat(sprintf("total: %.2fs\n", total))

classify <- function(name) {
  if (grepl("CppODE|solveODE|\\.Call", name)) return("ode_solve")
  if (grepl("res\\b|evalConditionResidual|objlist|constraintL2|normL2|Reduce|lapply|match\\.num|\\.bmm|data\\.table|objframe|parvec|prdframe|getDerivs|\\+\\.objfn|\\*\\.fn|\\+\\.objlist|\\+\\.fn", name))
    return("r_glue")
  "other"
}
prof$family <- vapply(rownames(prof), classify, "")
agg <- aggregate(self.time ~ family, prof, sum)
agg$pct <- 100 * agg$self.time / total
agg <- agg[order(-agg$self.time), ]
print(agg, row.names = FALSE)

cat("\n--- top 20 (N=64) ---\n")
top <- prof[order(-prof$self.time),
            c("self.time","self.pct","family")][1:20, ]
top$fun <- rownames(top)
print(top[, c("fun","self.time","self.pct","family")], row.names = FALSE)

## profvis: same Rprof samples, but flame-graph UI with explicit .Call time
## attribution. Lets us see how much wall time goes into CppODE's .Call
## entries vs the per-condition R glue around it.
if (requireNamespace("profvis", quietly = TRUE)) {
  cat("\n=== profvis: 200 obj_big() calls (N_cond = 64) ===\n")
  pv_N64 <- profvis::profvis(
    for (i in 1:200) invisible(obj_big(pouter)),
    interval = 0.005
  )
  html_path <- file.path(tempdir(), "profvis_normL2_N64.html")
  ## selfcontained = FALSE: HTML + lib/ sidecar, no pandoc needed.
  htmlwidgets::saveWidget(pv_N64, html_path, selfcontained = FALSE)
  assign("pv_N64", pv_N64, envir = .GlobalEnv)
  cat(sprintf("flame graph saved: %s\n", html_path))
  cat("interactive view:  print(pv_N64)   (or browseURL(html_path))\n")

  ## --- aggregate profvis samples by .Call vs R-glue ---
  ## profvis stores parsed prof data; we walk its sample frames and
  ## attribute each sample to the topmost .Call (= native C++ time) or
  ## to the topmost R frame above it.
  prof_df <- pv_N64$x$message$prof
  if (!is.null(prof_df)) {
    interval_ms <- pv_N64$x$message$interval %||% 5
    by_sample <- split(prof_df, prof_df$time)
    classify_stack <- function(s) {
      labs <- s$label
      if (any(grepl("\\.Call|CppODE::solveODE|CppODE_", labs))) return("cpp_call")
      if (any(grepl("tryCatch", labs))) return("tryCatch_glue")
      if (any(grepl("intersect|setdiff|union|unique|\\.set_ops", labs))) return("set_ops")
      if (any(grepl("parvec|c\\.parvec|as\\.parvec|\\[\\.parvec", labs))) return("parvec")
      if (any(grepl("lapply|do\\.call|\\*\\.fn|Reduce|FUN", labs))) return("compose_glue")
      "other"
    }
    cats <- vapply(by_sample, classify_stack, "")
    tab <- sort(table(cats), decreasing = TRUE)
    tot <- sum(tab) * interval_ms
    cat(sprintf("\nprofvis sample breakdown (interval = %g ms, %d samples = %.2fs):\n",
                interval_ms, sum(tab), tot/1000))
    for (k in names(tab))
      cat(sprintf("  %-15s %5d samples  %6.0f ms  (%5.1f%%)\n",
                  k, tab[[k]], tab[[k]] * interval_ms,
                  100 * tab[[k]] / sum(tab)))
  }
} else {
  cat("\n(profvis not installed; skip flame graph. install.packages('profvis'))\n")
}

cat("\nDONE\n")
