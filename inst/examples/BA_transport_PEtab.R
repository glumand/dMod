## Round-trip the BileAcid (TCA transport) model through PEtab v2
## ---------------------------------------------------------------
## 1. Build the dMod-native model (reactions, observables, steady-state
##    trafo, log10 outer parameters, branched closed/open conditions).
## 2. Fit it on `badata`.
## 3. Export the fitted problem to PEtab v2 — `exportPEtab` symbolically
##    decomposes the trafo `p` into parameters / conditions / SBML
##    initialAssignments. The `10^(...)` wraps stay intact in
##    conditions.tsv:targetValue / SBML initialAssignment, so PEtab
##    parameters live on the linear log10-domain (parameterScale = "lin").
## 4. Re-import the YAML manifest and verify the objective at the same
##    *outer* parameter vector matches the dMod-native one.

rm(list = ls())
library(dMod)
library(dplyr)

setwd(tempdir())  # generated artefacts (.cpp/.so/.tsv/.xml) land here
set.seed(5555)

## --- 1. dMod-native model ------------------------------------------------

reactions <- eqnlist() %>%
  addReaction("TCA_buffer", "TCA_cell", rate = "k_import*TCA_buffer",       description = "Uptake") %>%
  addReaction("TCA_cell",   "TCA_buffer", rate = "k_export_sinus*TCA_cell", description = "Sinusoidal export") %>%
  addReaction("TCA_cell",   "TCA_cana",   rate = "k_export_cana*TCA_cell",  description = "Canalicular export") %>%
  addReaction("TCA_cana",   "TCA_buffer", rate = "k_reflux*TCA_cana",       description = "Reflux into the buffer")

mymodel <- odemodel(reactions, modelname = "bamodel",
                    compile = FALSE, solver = "Sundials")
x <- Xs(mymodel)

observables <- eqnvec(buffer   = "s*TCA_buffer",
                      cellular = "s*(TCA_cana + TCA_cell)")
g <- Y(observables, f = x, condition = NULL,
       compile = FALSE, modelname = "obsfn_bamodel", attach.input = TRUE)

innerpars  <- getParameters(x, g)
mysteadies <- steadyStates(reactions, forcings = "k_import")

trafo <- eqnvec() %>%
  define("x~x", x = innerpars) %>%
  define("TCA_buffer~0") %>%
  define("x~y", x = names(mysteadies), y = mysteadies) %>%
  insert("x~10^y", x = .currentSymbols, y = toupper(.currentSymbols)) %>%
  branch(conditions = c("closed", "open")) %>%
  define("k_reflux~10^K_REFLUX_OPEN", conditionMatch = "open") %>%
  insert("S~0")

p <- P(trafo, modelname = "parfn", compile = FALSE)

compile(g, x, p, cores = 5)

## --- 2. Data + fit -------------------------------------------------------

data(badata)
data <- as.datalist(badata)
## badata carries per-row sigmas in the `sigma` column. exportPEtab
## auto-encodes them via PEtab's `noiseParameter1_<obsId>` placeholder,
## writing per-row sigma values into `noiseParameters` of measurements.tsv
## — the round-trip preserves the data-side noise model exactly.

outerpars <- getParameters(p)
pouter    <- structure(rep(-1, length(outerpars)), names = outerpars)

obj <- normL2(data, g * x * p)
fit <- trust(obj, pouter, rinit = 0.1, rmax = 5,
             iterlim = 500, printIter = TRUE)
phat <- fit$argument

cat("dMod-native obj(phat) =", obj(phat)$value, "\n")

times <- seq(0, 45, len = 300)
plot((g*x*p)(times, phat), data)

## --- 3. Export to PEtab v2 ----------------------------------------------
## `exportPEtab` reads `getEquations(p)` and decomposes each per-condition
## eqnvec, classifying each LHS:
##   - state with constant numeric RHS  -> SBML <initialConcentration>
##   - state with constant symbolic RHS -> SBML <initialAssignment>
##   - inner_par varying across conds   -> conditions.tsv parameter column
##   - inner_par with constant RHS      -> conditions.tsv (or SBML default)
## In v2 the `10^(...)` chain-rule wraps from `p` are kept verbatim in
## conditions.tsv:targetValue / SBML initialAssignment (PEtab v2 dropped
## parameterScale, so the trafo lives entirely in the model).
## parameters.tsv contains only estimated parameters (matching the v2
## reference test suite); inner symbols referenced by conditions.tsv
## targetIds (k_import, k_export_*, k_reflux) and identity-mapped outer
## pars (s) live in the SBML <parameter> table and are recovered by the
## importer's SBML-default-as-fixed gap.
out_dir  <- file.path(tempdir(), "BA_PEtab_export")
yaml_out <- exportPEtab(
  data           = data,
  reactions      = reactions,
  observables    = observables,
  p              = p,
  pouter         = phat,
  parameterScale = "lin",
  modelID        = "BileAcid",
  dir            = out_dir,
  overwrite      = TRUE
)

cat("\nWritten files:\n  ", paste(list.files(out_dir), collapse = "\n  "), "\n", sep = "")

## --- 4. Re-import and verify --------------------------------------------

petab <- importPEtab(yaml_out, solver = "CppODE", cores = 4)

cat("\nReimported bestfit names:       ", paste(names(petab$bestfit), collapse = ", "), "\n")
cat("Reimported parameterScales:     ", paste(attr(petab$bestfit, "petab_scales"), collapse = ", "), "\n")
cat("Native pouter names:            ", paste(names(phat), collapse = ", "), "\n")

stopifnot(setequal(names(petab$bestfit), names(phat)))

## Compare objective values at the same outer point. `petab$obj` has the
## PEtab-fixed parameters baked in, so a single argument suffices.
val_native <- obj(phat)$value
val_petab  <- petab$obj(petab$bestfit)$value
cat("\nNative    obj(phat) =", val_native, "\n")
cat("Re-import obj(phat) =", val_petab, "\n")
cat("Difference          =", val_petab - val_native, "\n")

stopifnot(abs(val_petab - val_native) < 1e-3)

## And compare predictions across conditions. `petab$prd` is the raw
## composite g * x * p; pass the estimated bestfit plus any fixed outer
## parameters of `p` (here `s = 1`, an identity outer that the native
## trafo collapsed via `insert("S~0")`). `petab_meta$fixed` holds exactly
## those — inner symbols pulled from SBML defaults are stashed in
## `petab_meta$sbml_only_pars` and not needed at the prd boundary.
fixed_petab <- attr(petab, "petab_meta")$fixed
times       <- seq(0, 45, length.out = 200)
pred_native <- (g * x * p)(times, phat)
pred_petab  <- petab$prd(times, c(petab$bestfit, fixed_petab))
plot(pred_petab, data)
