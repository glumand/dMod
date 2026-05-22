## ============================================================================
## BA_transport Bayes: posterior sampling on a mechanistic kinetic model
## ============================================================================
##
## Companion to inst/examples/BA_transport.R.
##
##     mcmc(flatTarget(normL2(data, prd), priorObj, priorSample),
##          sequenceType = "sequential", ...)
##
## Three deliberate choices in this example:
##
##   1. `deriv2 = FALSE` everywhere in the model build. The flat-target +
##      Langevin SMC path uses only first-order sensitivities — the
##      Gauss-Newton Hessian `J^T W J` (with J = dpred/dpar from the
##      first-order sensitivity ODEs) is both the trust() Hessian and the
##      Langevin metric. The `_s2` codegen pipeline is skipped entirely.
##      Compile cost drops 2-3x vs the NLME example.
##
##   2. Weakly-informative Gaussian prior `constraintL2(0, sigma = 1)`
##      attached to the objective. The closed-condition simple trafo has
##      a known structural non-identifiability between (S, TCA_CELL,
##      TCA_CANA) — BA_transport.R flags it at line 144 as "the triple
##      compensate by structure". A truly improper flat prior would give
##      an improper posterior on that ridge (the chain drifts to infinity
##      along the flat direction); the weakly-informative prior keeps the
##      problem proper without meaningfully constraining the identifiable
##      directions. SD = 1 on log10 scale means 68% mass within factor-of-10
##      per parameter. Loose enough that identifiable directions are barely
##      moved; tight enough that the Hessian at the MAP is well-conditioned.
##
##      This is the dMod-idiomatic Bayes prior. constraintL2 plays the
##      double role of (a) regularising mstrust during the MAP search and
##      (b) being the priorObj for mcmc(). Identical reference point,
##      identical sigma — no double counting because we add only ONE prior
##      contribution to the joint objective.
##
##   3. SMC with a MAP-Laplace cloud, fixed-MAP Langevin preconditioner,
##      adaptive-ESS tempering. SMC is the right tool here because
##      tempering at low beta lets particles cover the non-identifiable
##      ridge directions, then the move kernel tightens them onto the
##      data-supported subspace as beta -> 1. Returns particles + free
##      log-evidence logZ. Single-chain MALA from the MAP does NOT mix
##      on this posterior (the sharp identifiable directions + the soft
##      ridge break the local-curvature assumption; acceptance collapses).
##      Use SMC unless you have an identifiable reparameterisation in
##      hand — see the optional block at the bottom for that route.
##
## Wallclock guidance: compile + mstrust + SMC ~1-2 min.
##
## Returns map_frame, mapfit, chain_smc and the ggplot panels
## p_marginals / p_trace / p_pairs in the global env.
## ============================================================================

rm(list = ls(all.names = TRUE))
setwd(tempdir())
set.seed(1)

library(dMod)
library(dplyr)
library(ggplot2)

ncores <- max(1L, parallel::detectCores() - 1L)

## ----------------------------------------------------------------------------
## 1. Data + 4-reaction kinetic model (closed condition only)
## ----------------------------------------------------------------------------
data(badata)
data <- badata %>% subset(condition == "closed") %>% as.datalist()

reactions <- eqnlist() %>%
  addReaction("TCA_buffer", "TCA_cell",  rate = "k_import * TCA_buffer", description = "Uptake") %>%
  addReaction("TCA_cell",   "TCA_buffer", rate = "k_export_sinus * TCA_cell", description = "Sinusoidal export") %>%
  addReaction("TCA_cell",   "TCA_cana",   rate = "k_export_cana * TCA_cell", description = "Canalicular export") %>%
  addReaction("TCA_cana",   "TCA_buffer", rate = "k_reflux * TCA_cana", description = "Reflux into the buffer")

mymodel <- odemodel(reactions, modelname = "ba_bayes_ode",
                    compile = FALSE, solver = "CppODE")
x <- Xs(mymodel)

observables <- eqnvec(buffer   = "s * TCA_buffer",
                      cellular = "s * (TCA_cana + TCA_cell)")
g <- Y(observables, f = x, condition = NULL,
       modelname = "ba_bayes_obs", attach.input = TRUE,
       compile = FALSE)

# Tighten model assumptions with steady state constraint
mysteadies <- steadyStates(reactions, forcings = "k_import")

innerpars <- getParameters(x, g)
trafo <- eqnvec() %>%
  define("x~x", x = innerpars) %>% # identity
  define("TCA_buffer~0") %>%
  define("x~y", x = names(mysteadies), y = mysteadies) %>% 
  insert("x~10^y", x = .currentSymbols, y = toupper(.currentSymbols)) %>% 
  branch(conditions = c("closed", "open")) %>% 
  define("k_reflux~10^K_REFLUX_OPEN", conditionMatch = "open") %>% 
  insert("S~0") # fixed structural non identifiablility

p <- P(trafo, condition = "closed", modelname = "ba_bayes_p",
       compile = FALSE)

compile(g, x, p, output = "ba_bayes", cores = 6)

prd <- g * x * p
outerpars <- getParameters(p)
cat(sprintf("Outer parameters (%d): %s\n",
            length(outerpars), paste(outerpars, collapse = ", ")))

## ----------------------------------------------------------------------------
## 2. Likelihood + weakly-informative prior on log10 scale
## ----------------------------------------------------------------------------
prior_ref <- setNames(rep(-1, length(outerpars)), outerpars)
prior_obj <- constraintL2(prior_ref, sigma = 4)        # SD = 4 on log10
lik_obj   <- normL2(data, prd)
obj       <- lik_obj + prior_obj                       # MAP target

## ----------------------------------------------------------------------------
## 3. MAP via mstrust (the prior keeps mstrust off the flat ridge)
## ----------------------------------------------------------------------------
init <- structure(runif(length(outerpars), -1, 0), names = outerpars)

cat("\n== mstrust(20 fits) for the MAP ==\n")
ms <- mstrust(obj, init, sd = 2,
              studyname = "ba_bayes_ms", cores = ncores,
              fits = 20, iterlim = 500)
map_frame <- as.parframe(ms)
print(plotValues(map_frame))

mapfit <- trust(obj, parinit = as.parvec(map_frame),
                rinit = 0.1, rmax = 10, iterlim = 100)
cat(sprintf("\nMAP OFV = %.4f\n", mapfit$value))
cat("MAP arg =\n"); print(round(mapfit$argument, 3))

hess_eig <- eigen(mapfit$hessian, only.values = TRUE)$values
cat(sprintf("Hessian eigenvalue range: %.2e .. %.2e  (cond = %.2e)\n",
            min(hess_eig), max(hess_eig), max(hess_eig) / min(hess_eig)))

## ----------------------------------------------------------------------------
## 4. SMC with MAP-Laplace cloud + fixed-MAP Langevin preconditioner
##
## The MAP-Laplace cloud seeds the level-0 particles at N(MAP, 2 * H^{-1}).
## At beta = 0 they cover the high-density region; adaptive-ESS bisection
## then anneals beta -> 1, tightening onto the full posterior. SMC handles
## the non-identifiable ridge naturally because high-temperature levels
## allow long jumps along the flat direction.
## ----------------------------------------------------------------------------
mapArg <- mapfit$argument
mapCov <- 2 * solve(mapfit$hessian)                    # -2 log L -> log p
mapL   <- chol((mapCov + t(mapCov)) / 2)

prior_obj <- constraintL2(mapArg, sigma = 1)        # SD = 1 on log10

priorSample <- function(n) {
  Z <- matrix(rnorm(n * length(mapArg)), nrow = n)
  out <- sweep(Z %*% mapL, 2L, mapArg, FUN = "+")
  colnames(out) <- names(mapArg)
  out
}

tgt <- flatTarget(likObj      = lik_obj,
                  priorObj    = prior_obj,
                  priorSample = priorSample)

cat("\n== mcmc(SMC, langevin + fixed-MAP-metric) ==\n")
t_smc <- system.time(
  chain_smc <- mcmc(target           = tgt,
                    sequenceType     = "sequential",
                    sequenceSchedule = "adaptiveEss",
                    moveType         = "langevin",
                    metric           = "fixed",
                    populationSize   = 2000L,
                    moveControl      = langevinControl(GFixed = mapfit$hessian),
                    sequenceControl  = smcControl(verbose = TRUE),
                    cores            = ncores)
)
cat(sprintf("[SMC]  elapsed = %.1fs, levels = %d, logZ = %+8.3f\n",
            t_smc[["elapsed"]], chain_smc$nLevels, chain_smc$logEvidence))

## ----------------------------------------------------------------------------
## 5. Posterior summary
## ----------------------------------------------------------------------------
postMean <- colMeans(chain_smc$samples)
postSD   <- apply(chain_smc$samples, 2L, sd)
summary_tbl <- data.frame(
  parameter     = names(mapArg),
  MAP           = round(mapArg, 4),
  posteriorMean = round(postMean, 4),
  posteriorSD   = round(postSD, 4),
  diffOverSD    = round((postMean - mapArg) / pmax(postSD, 1e-8), 2),
  row.names     = NULL)
cat("\n== MAP vs posterior summary ==\n")
print(summary_tbl)

## ----------------------------------------------------------------------------
## 6. Diagnostic plots (ggplot objects)
## ----------------------------------------------------------------------------
p_marginals <- plot(chain_smc)
p_trace     <- plotTrace(chain_smc)
p_pairs     <- plotPairs(chain_smc)

print(p_marginals)
print(p_trace)
print(p_pairs)

## ----------------------------------------------------------------------------
## Alternative sampler configurations
##
## The default above is moveType = "langevin" + metric = "fixed", which is
## the cheapest reasonable choice (one gradient per inner move, MAP Hessian
## as preconditioner stays valid across particles). Three legitimate
## alternatives, ordered by per-move compute cost:
##
## (i) Langevin + riemannFisher
##     Local-curvature MALA. The Fisher info is recomputed at every step.
##     Cost: 1 grad + 1 Hessian per step (Hessian = GN from first-order
##     sensitivities, so still cheap). Helps when the posterior curvature
##     varies strongly across the posterior; on a Gaussian-like posterior
##     it's roughly the same answer as fixed-MAP but with extra noise from
##     the per-step Fisher.
##
##         chain_alt <- mcmc(target           = tgt,
##                           sequenceType     = "sequential",
##                           sequenceSchedule = "adaptiveEss",
##                           moveType         = "langevin",
##                           metric           = "riemannFisher",
##                           populationSize   = 300L,
##                           sequenceControl  = smcControl(malaSteps = 5L,
##                                                          essThreshold = 0.7),
##                           cores            = ncores)
##
## (ii) HMC + riemannFisher (Fisher mass per particle, frozen across trajectory)
##     Each particle gets a leapfrog trajectory whose mass matrix is the
##     local Fisher info at the particle's current position. Much better
##     mixing per move than MALA (long trajectory traverses the level set).
##     Cost per move: leapfrogSteps grads + 1 Hessian factorisation.
##     Default leapfrogSteps = 10, so ~10x Langevin per move; gains come
##     from being able to drop malaSteps to 1-2 instead of 5.
##
##         chain_alt <- mcmc(target = tgt, sequenceType = "sequential",
##                           sequenceSchedule = "adaptiveEss",
##                           moveType = "hmc", metric = "riemannFisher",
##                           populationSize = 300L,
##                           moveControl = hmcControl(leapfrogSteps = 10L),
##                           sequenceControl = smcControl(malaSteps = 2L,
##                                                         essThreshold = 0.7),
##                           cores = ncores)
##
## (iii) HMC / NUTS + fixed-MAP mass (the dMod analogue of Stan's dense_e)
##     One Hessian factorisation up front, leapfrog uses MAP curvature for
##     all particles. NUTS auto-tunes trajectory length, saving the user
##     from picking leapfrogSteps. This is the most efficient long-trajectory
##     option for posteriors that are well-approximated by a Gaussian.
##
##         chain_alt <- mcmc(target = tgt, sequenceType = "sequential",
##                           moveType = "nuts", metric = "fixed",
##                           populationSize = 300L,
##                           moveControl = nutsControl(MFixed = mapfit$hessian),
##                           sequenceControl = smcControl(malaSteps = 1L,
##                                                         essThreshold = 0.7),
##                           cores = ncores)
##
## When to pick what:
##   * Posterior near-Gaussian, MAP Hessian PSD, dim < ~20 -> langevin + fixed
##     (the default in this file). Cheapest, plenty of mixing per move.
##   * Posterior curvature varies strongly -> langevin + riemannFisher.
##   * Higher dim (>30) or stronger correlations -> hmc + fixed or nuts.
##   * Diagnostics in (e.g. omega axis flat in NLME) -> nuts + fixed for
##     reliable mixing without per-particle tuning.
##   * Avoid HMC + riemannFisher unless you have a reason: per-particle
##     Fisher in tail regions can be ill-conditioned and break the
##     leapfrog integrator.
##
## ----------------------------------------------------------------------------
## Model reduction from Bayesian posterior diagnostics
##
## The BA_transport closed-only model carries a structural non-identifiability
## between (S, TCA_CELL, TCA_CANA). With a flat prior the posterior would be
## improper (infinite mass on the compensation ridge). The sigma = 1 prior
## above is a workaround. The principled fix is to *reduce the model* so the
## non-identifiability disappears. The diagnostics live entirely in the
## posterior — no profile likelihood, no separate frequentist analysis.
##
## Step 1: posterior shrinkage to find weakly-informed parameters.
##
##   The shrinkage statistic compares the posterior SD per parameter to the
##   prior SD. Values close to 0 mean the data has not moved the marginal
##   (the prior is doing all the work); values close to 1 mean the data
##   dominates.
##
##     prior_draws <- priorSample(2000L)
##     post_sd     <- apply(chain_smc$samples, 2L, sd)
##     prior_sd    <- apply(prior_draws,        2L, sd)
##     shrinkage   <- 1 - (post_sd / prior_sd)^2
##     print(round(sort(shrinkage), 3))
##
##   For BA_transport closed-only this typically yields shrinkage ~ 0.9 for
##   K_REFLUX, ~ 0.7 for K_IMPORT, and ~ 0.0-0.2 for the (S, TCA_CELL,
##   TCA_CANA, K_EXPORT_*) cluster. The latter are the parameters whose
##   marginals are essentially the prior projected onto that axis.
##
## Step 2: posterior covariance eigenanalysis to find the *direction* of
##   the non-identifiability (not just per-parameter, but the linear
##   combination of parameters that the data does not constrain).
##
##     cov_post  <- cov(chain_smc$samples)
##     eig_post  <- eigen(cov_post, symmetric = TRUE)
##     # Largest eigenvalue = widest posterior direction = least informed
##     direction_idx <- which.max(eig_post$values)
##     ridge_dir <- setNames(eig_post$vectors[, direction_idx],
##                            colnames(chain_smc$samples))
##     cat("Posterior ridge (largest-variance eigenvector):\n")
##     print(round(ridge_dir, 3))
##     cat(sprintf("Variance along ridge: %.2f   (vs smallest: %.2e)\n",
##                 eig_post$values[direction_idx], min(eig_post$values)))
##
##   For BA_transport this eigenvector points roughly along
##   (+TCA_CELL, +TCA_CANA, -S, ...) — the (S, TCA_CELL, TCA_CANA)
##   compensation flagged at BA_transport.R line 144. The eigenvalue ratio
##   max/min quantifies how anisotropic the posterior is; values > 1e3
##   indicate a strong ridge worth reducing.
##
## Step 3: 2D pair-plot confirmation.
##
##   plotPairs(chain_smc) draws the corner plot. Strongly correlated
##   pairs show up as narrow diagonal lines in the 2D marginals — those
##   are the ridge components. For BA_transport you will see tight
##   correlation between (S, TCA_CELL), (S, TCA_CANA), and (TCA_CELL,
##   TCA_CANA), confirming the same ridge identified by the eigenvector.
##
## Step 4: structural reduction implied by the ridge direction.
##
##   The ridge direction tells you *what to reduce*: in BA_transport the
##   coefficients on (S, TCA_CELL, TCA_CANA) being large-and-balanced
##   means these three together carry only one degree of identifiable
##   freedom. Two principled fixes:
##
##   4a. Steady-state substitution for state initial conditions.
##
##         mysteadies <- steadyStates(reactions, forcings = "k_import")
##         trafo <- eqnvec() %>%
##           define("x~x", x = innerpars) %>%
##           define("TCA_buffer~0") %>%
##           define("x~y", x = names(mysteadies), y = mysteadies) %>%
##           insert("x~10^y", x = .currentSymbols, y = toupper(.currentSymbols))
##
##       Two outer parameters disappear (TCA_CELL, TCA_CANA become
##       algebraic functions of the rates and k_import).
##
##   4b. Absolute-scale fixing.
##
##         trafo <- trafo %>% insert("S~0")
##
##       Drops the remaining ridge component. After 4a+4b the closed-only
##       model has 4 outer parameters (the four rate constants), all
##       identifiable, Hessian condition < 100 at the MAP.
##       `priorObj = NULL` in flatTarget() now gives a proper posterior.
##
## Step 5: validation via re-sampling.
##
##   Re-run SMC on the reduced model and verify:
##     * shrinkage > 0.5 for every remaining parameter
##     * cov_post eigenvalue spread (max/min) < ~100
##     * plotPairs shows no narrow diagonal lines
##
##   If yes, model reduction succeeded. If no, repeat steps 1-4 on the
##   reduced parameterisation — there may be a second (smaller) ridge.
##
## Step 6: comparing reduced vs full model via Bayes factor.
##
##   Both runs report logZ (chain_smc$logEvidence). The Bayes factor
##   `exp(logZ_reduced - logZ_full)` quantifies which model the data
##   prefers, accounting for the parameter-count complexity automatically.
##   A reduced model that is "right" (removes only redundant degrees of
##   freedom) typically has logZ_reduced > logZ_full because Occam's
##   razor in the marginal-likelihood integral rewards lower-dimensional
##   models with the same data fit.
##
##     bf <- exp(chain_reduced$logEvidence - chain_smc$logEvidence)
##     cat(sprintf("Bayes factor reduced/full = %.2e\n", bf))
##
## Step 7: only if 1-6 are not possible, tighten the prior.
##
##   Sometimes the non-identifiability is data-quantity-limited rather
##   than structure-limited (more data would identify it). In that case
##   no structural reduction helps; an informative prior is the honest
##   way to keep the posterior proper:
##
##     prior_obj <- constraintL2(prior_ref, sigma = 0.5)     # SD = 0.5 on log10
##
##   sigma = 0.5 means the prior contributes meaningfully (95% mass
##   within factor 10 of the reference). Document this as informative,
##   not weakly-informative, and report the posterior shrinkage so the
##   reader knows how much the prior is contributing.
## ----------------------------------------------------------------------------

