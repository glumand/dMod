## Same breakdown for the Boehm benchmark model (PEtab),
## a STAT5 phosphorylation model — 8 dynamic states, ~17 fit pars,
## 3 observables, 16 timepoints. Bigger than BA so we can see how
## the R-glue share changes as ODE/sens work grows.

suppressPackageStartupMessages({
  library(dMod)
  library(microbenchmark)
})

setwd(tempdir())
set.seed(5555)

yaml <- file.path("/home/simon/Documents/Projects/dMod",
                  "BenchmarkModels/Boehm_JProteomeRes2014",
                  "Boehm_JProteomeRes2014.yaml")
cat("Importing Boehm PEtab ... ")
t0 <- Sys.time()
pp <- importPEtab(yaml, solver = "CppODE", compile = TRUE, cores = 4)
cat(sprintf("done (%.1fs)\n", as.numeric(Sys.time() - t0, units = "secs")))

obj <- pp$obj
prd <- pp$prd
pouter <- pp$bestfit
times <- sort(unique(unlist(lapply(pp$dataList, `[[`, "time"))))
times <- sort(unique(c(0, times)))

cat(sprintf("dim: pars=%d  conds=%d  ndata=%d\n",
            length(pouter), length(pp$dataList),
            sum(vapply(pp$dataList, nrow, 0L))))

## warmup
invisible(prd(times, pouter, deriv = TRUE))
invisible(obj(pouter))

cat("\n=== microbenchmark ===\n")
mb <- microbenchmark(
  prd_value_only = prd(times, pouter, deriv = FALSE),
  prd_deriv1     = prd(times, pouter, deriv = TRUE,  deriv2 = FALSE),
  normL2_d1      = obj(pouter),
  times = 30L, unit = "ms"
)
print(summary(mb)[, c("expr", "min", "median", "max")], row.names = FALSE)

s <- summary(mb); med <- setNames(s$median, s$expr)
cat(sprintf("\nglue share (d1) = %.1f%%\n",
            100 * (med["normL2_d1"] - med["prd_deriv1"]) / med["normL2_d1"]))

cat("\n=== Rprof over 20 trust starts ===\n")
prof_file <- tempfile(fileext = ".out")
sds <- runif(20, 0.2, 0.5)
Rprof(prof_file, interval = 0.005)
for (i in seq_along(sds)) {
  pinit <- pouter + rnorm(length(pouter), sd = sds[i])
  try(trust(obj, pinit, rinit = 0.1, rmax = 10,
            iterlim = 80, printIter = FALSE), silent = TRUE)
}
Rprof(NULL)

prof <- summaryRprof(prof_file)$by.self
total <- sum(prof$self.time)
cat(sprintf("total sampled time: %.2fs\n", total))

classify <- function(name) {
  if (grepl("CppODE|solveODE|\\.Call", name)) return("ode_solve")
  if (grepl("res\\b|evalConditionResidual|objlist|constraintL2|normL2|Reduce|lapply|match\\.num|\\.bmm|data\\.table|objframe|parvec|prdframe|getDerivs|\\+\\.objfn|\\*\\.fn", name))
    return("r_glue")
  if (grepl("^trust\\b|trustOptim|solve\\b|qr\\.", name)) return("trust_step")
  "other"
}
prof$family <- vapply(rownames(prof), classify, "")
agg <- aggregate(self.time ~ family, prof, sum)
agg$pct <- 100 * agg$self.time / total
agg <- agg[order(-agg$self.time), ]
cat("\n--- by family ---\n")
print(agg, row.names = FALSE)

cat("\n--- top 15 ---\n")
top <- prof[order(-prof$self.time), c("self.time", "self.pct", "family")][1:15, ]
top$fun <- rownames(top)
print(top[, c("fun", "self.time", "self.pct", "family")], row.names = FALSE)

cat("\nDONE\n")
