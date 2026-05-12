## Profile the classical dMod fit pipeline
##   trust(normL2(data, g*x*p) + constraintL2(...) [+ datapointL2(...)])
##
## Goal: how much wallclock is ODE solve + sensitivities (inside the
## prediction function x), and how much is R-glue (res() / chain-rule /
## objlist + assembly in normL2, constraintL2, datapointL2, trust())?
##
## If the R-glue is < 10% the rewrite is not worth it.
## If 25%+ a shared C++ residual kernel + a C++ trust would help.

suppressPackageStartupMessages({
  library(dMod)
  library(dplyr)
  library(microbenchmark)
})

setwd(tempdir())
set.seed(5555)

cat("R:", R.version.string, "\n")
cat("dMod:", as.character(packageVersion("dMod")), "\n")
cat("CppODE:", as.character(packageVersion("CppODE")), "\n\n")

## ---- BA_transport model (small, realistic) ----
data(badata)
data_full <- as.datalist(badata)

reactions <- eqnlist() |>
  addReaction("TCA_buffer", "TCA_cell",  rate = "k_import*TCA_buffer") |>
  addReaction("TCA_cell",   "TCA_buffer", rate = "k_export_sinus*TCA_cell") |>
  addReaction("TCA_cell",   "TCA_cana",   rate = "k_export_cana*TCA_cell") |>
  addReaction("TCA_cana",   "TCA_buffer", rate = "k_reflux*TCA_cana")

mymodel <- odemodel(reactions, modelname = "bamodel_bench",
                    solver = "CppODE", deriv2 = TRUE, compile = FALSE)
x <- Xs(mymodel,
        optionsOde  = list(atol = 1e-8, rtol = 1e-8),
        optionsSens = list(atol = 1e-8, rtol = 1e-8))

observables <- eqnvec(buffer   = "s*TCA_buffer",
                      cellular = "s*(TCA_cana + TCA_cell)")
g <- Y(observables, f = x, condition = NULL, deriv2 = TRUE,
       compile = FALSE, modelname = "obsfn_bamodel_bench", attach.input = FALSE)

innerpars <- getParameters(x, g)
trafo <- NULL |>
  define("x~x", x = innerpars) |>
  define("TCA_buffer~0") |>
  insert("x~exp10(y)", x = .currentSymbols, y = toupper(.currentSymbols))

p_closed <- P(trafo, condition = "closed", compile = FALSE,
              modelname = "p_closed_bench", deriv2 = TRUE)
trafo_open <- getEquations(p_closed, conditions = "closed") |>
  insert("K_REFLUX~K_REFLUX_OPEN")
p_open <- P(trafo_open, condition = "open", compile = FALSE,
            modelname = "p_open_bench", deriv2 = TRUE)
p <- p_closed + p_open

cat("Compiling ... ")
t0 <- Sys.time()
compile(g, x, p, output = "bamodel_bench", cores = 4)
cat(sprintf("done (%.1fs)\n", as.numeric(Sys.time() - t0, units = "secs")))

prd <- g * x * p
outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)

times <- sort(unique(c(0, unlist(lapply(data_full, `[[`, "time")))))

obj_data   <- normL2(data_full, g * x * p)
obj_prior  <- constraintL2(pouter, sigma = 4)
obj_data_prior <- obj_data + obj_prior

## warmup
invisible(prd(times, pouter, deriv = TRUE,  deriv2 = FALSE))
invisible(prd(times, pouter, deriv = TRUE,  deriv2 = TRUE))
invisible(obj_data(pouter))
invisible(obj_data_prior(pouter))
invisible(obj_data_prior(pouter, deriv2 = TRUE))

## ---- microbenchmark each layer ----
cat("\n=== microbenchmark (median wallclock per call) ===\n")
mb <- microbenchmark(
  prd_value_only       = prd(times, pouter, deriv = FALSE),
  prd_deriv1           = prd(times, pouter, deriv = TRUE,  deriv2 = FALSE),
  prd_deriv2           = prd(times, pouter, deriv = TRUE,  deriv2 = TRUE),
  normL2_d1            = obj_data(pouter),
  normL2_d2            = obj_data(pouter, deriv2 = TRUE),
  normL2_plus_prior_d1 = obj_data_prior(pouter),
  normL2_plus_prior_d2 = obj_data_prior(pouter, deriv2 = TRUE),
  times = 40L,
  unit  = "ms"
)
print(summary(mb)[, c("expr", "min", "median", "max")], row.names = FALSE)

s <- summary(mb)
med <- setNames(s$median, s$expr)
cat(sprintf(
  "\nglue share, d1:  (normL2 - prd_d1) / normL2 = %.1f%%\n",
  100 * (med["normL2_d1"] - med["prd_deriv1"]) / med["normL2_d1"]))
cat(sprintf(
  "glue share, d2:  (normL2 - prd_d2) / normL2 = %.1f%%\n",
  100 * (med["normL2_d2"] - med["prd_deriv2"]) / med["normL2_d2"]))

## --- repeat with INFLATED data + dense time grid to simulate a bigger study
cat("\n=== inflated data (8x replicated, 200 timepoint grid) ===\n")
big <- do.call(rbind, lapply(1:8, function(rep) {
  d <- badata
  d$value <- d$value + rnorm(nrow(d), 0, d$sigma * 0.2)
  d
}))
big_data <- as.datalist(big)
obj_big  <- normL2(big_data, g * x * p) + constraintL2(pouter, sigma = 4)
times_big <- sort(unique(c(0, unlist(lapply(big_data, `[[`, "time")))))
invisible(obj_big(pouter))
mb_big <- microbenchmark(
  prd_big_d1   = prd(times_big, pouter, deriv = TRUE),
  normL2_big   = obj_big(pouter),
  times = 30L, unit = "ms"
)
print(summary(mb_big)[, c("expr", "min", "median", "max")], row.names = FALSE)
s2 <- summary(mb_big); m2 <- setNames(s2$median, s2$expr)
cat(sprintf("\nglue share inflated: %.1f%%\n",
            100 * (m2["normL2_big"] - m2["prd_big_d1"]) / m2["normL2_big"]))

## ---- Rprof over a full trust() run ----
cat("\n=== Rprof over many trust(obj_data_prior, ...) starts ===\n")
prof_file <- tempfile(fileext = ".out")
sds <- runif(40, 0.3, 0.6)
Rprof(prof_file, interval = 0.005, line.profiling = FALSE)
for (i in seq_along(sds)) {
  pinit <- pouter + rnorm(length(pouter), sd = sds[i])
  fit <- try(trust(obj_data_prior, pinit, rinit = 0.1, rmax = 10,
                   iterlim = 200, printIter = FALSE), silent = TRUE)
}
Rprof(NULL)
cat(sprintf("last trust: iter=%d  value=%.4f\n",
            fit$iterations, fit$value))

prof <- summaryRprof(prof_file)$by.self
total <- sum(prof$self.time)
cat(sprintf("total sampled time: %.2fs\n", total))

## Classify samples by family
classify <- function(name) {
  ode_pattern <- c("\\.Call", "lsoda", "lsodes", "ode_", "Xs\\.", "CppODE",
                   "\\.C\\b", "solveODE", "compileAndLoad")
  glue_pattern <- c("res\\b", "evalConditionResidual", "objlist",
                    "constraintL2", "datapointL2", "normL2",
                    "Reduce", "lapply", "match\\.num", "\\.bmm",
                    "match\\b", "data\\.table", "as\\.data\\.frame",
                    "\\[\\.objframe", "objframe", "as\\.parvec",
                    "\\+\\.objfn", "\\*\\.fn",  "\\+\\.fn",
                    "getDerivs", "modify",  "prdframe")
  trust_pattern <- c("^trust\\b", "trustOptim", "newton", "solve\\b", "qr\\.")
  if (any(grepl(paste(ode_pattern,  collapse = "|"), name))) return("ode_solve")
  if (any(grepl(paste(trust_pattern,collapse = "|"), name))) return("trust_step")
  if (any(grepl(paste(glue_pattern, collapse = "|"), name))) return("r_glue")
  "other"
}
prof$family <- vapply(rownames(prof), classify, "")
agg <- aggregate(self.time ~ family, prof, sum)
agg$pct <- 100 * agg$self.time / total
agg <- agg[order(-agg$self.time), ]
cat("\n--- by family ---\n")
print(agg, row.names = FALSE)

cat("\n--- top 20 functions ---\n")
top <- prof[order(-prof$self.time), c("self.time", "self.pct", "family")][1:20, ]
top$fun <- rownames(top)
print(top[, c("fun", "self.time", "self.pct", "family")], row.names = FALSE)

cat("\nDONE\n")
