## Benchmark: sparse-grid quadrature node-generation + per-subject integrand
## evaluation on the Theophylline NLME problem (K_eta = 3, N = 12).
##
## Establishes the cost floor for one CM-step-1 evaluation in the ECM polish
## loop, before the Phase-2 batched evaluator lands. Compares against FOCEI's
## one-subject inner-solver evaluation as the "Laplace-cost-equivalent" baseline.

rm(list = ls(all.names = TRUE))
.dmod_root <- "/home/simon/Documents/Projects/dMod"
devtools::load_all(.dmod_root, quiet = TRUE)
setwd(tempdir())
unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE),
       force = TRUE)

set.seed(1)

# 1) Build the Theophylline model + data + joint (same as testing/focei_theophylline.R)
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
subjects <- sort(unique(Theoph$Subject))
N <- length(subjects)
doses <- vapply(subjects, function(s) {
  rec <- Theoph[Theoph$Subject == s, ][1, ]
  rec$Dose * rec$Wt
}, 0.0)
data_df <- data.frame(name = "y", time = Theoph$Time, value = Theoph$conc,
                      sigma = 0.5, condition = Theoph$Subject,
                      stringsAsFactors = FALSE)
dlist <- as.datalist(data_df)
reactions <- eqnlist()
reactions <- addReaction(reactions, "Ag", "",  "Ka * Ag",     "absorption")
reactions <- addReaction(reactions, "",   "Cc", "Ka * Ag / V", "appearance")
reactions <- addReaction(reactions, "Cc", "",  "Cl/V * Cc",   "elimination")
m <- odemodel(reactions, modelname = "bq_ode", compile = TRUE,
              solver = "CppODE", deriv2 = TRUE)
x <- Xs(m)
g <- Y(c(y = "Cc"), x, modelname = "bq_obs", compile = TRUE, deriv2 = TRUE)
trafo <- eqnvec(Ka = "Ka_pop * exp(eta_Ka)", V = "V_pop * exp(eta_V)",
                Cl = "Cl_pop * exp(eta_Cl)", Ag = "Ag_init", Cc = "0")
subj_table <- data.frame(
  eta_Ka = paste0("eta_Ka_", subjects),
  eta_V  = paste0("eta_V_",  subjects),
  eta_Cl = paste0("eta_Cl_", subjects),
  Ag_init = doses, row.names = subjects, stringsAsFactors = FALSE)
trafos <- branch(trafo, table = subj_table, apply = "insert")
p <- P(trafos, method = "explicit", compile = TRUE, modelname = "bq_p",
       deriv2 = TRUE)
om <- omega(eta = c("eta_Ka", "eta_V", "eta_Cl"), subjects = subjects)
joint <- normL2(dlist, g * x * p) + constraintL2(mu = 0, Omega = om)

# Starting values (close to NONMEM reference)
init <- c(Ka_pop = 1.5, V_pop = 32, Cl_pop = 3.0)
init[om$cholPars] <- c(log(0.5), log(0.1), log(0.3))
eta_init <- setNames(rep(0, N * 3),
                     c(paste0("eta_Ka_", subjects),
                       paste0("eta_V_",  subjects),
                       paste0("eta_Cl_", subjects)))
full_pars <- c(init, eta_init)

# 2) Build a FOCEI estimator and run one outer call to populate etaModes / H_i
focei_obj <- focei(joint, om, innerControl = list(rtol = 1e-7, maxit = 30))
cat("== FOCEI single outer evaluation (Laplace baseline) ==\n")
t_focei <- system.time({
  res_focei <- focei_obj(init)
})
print(t_focei)
cat(sprintf("  OFV = %.4f\n", res_focei$value))
diag <- attr(res_focei, "emDiag")
etaModes <- diag$etaStar            # N x K
cat("  etaModes (per subject):\n"); print(round(etaModes, 3))

# 3) Pull H_i for each subject from focei diag.
# emDiag stores eigenvalues (eigs) and the logdet of H_i, not H_i directly.
# Rebuild from one explicit joint call at the modes (same conv as compute_correction).
joint_full <- joint(full_pars[c(names(init), as.vector(om$subjectEtas))],
                    deriv = TRUE, conditions = subjects)
H_per_subject <- lapply(seq_len(N), function(i) {
  eta_i_names <- om$subjectEtas[i, ]
  joint_full$hessian[eta_i_names, eta_i_names, drop = FALSE]
})
# Replace etaModes from diag with what FOCEI converged to (zero on the first
# call since the cache starts at 0).
# For the benchmark we don't need converged modes; pick arbitrary nonzero modes
# to reflect realistic late-iteration state:
eta_modes_demo <- matrix(rnorm(N * 3, sd = 0.2), N, 3)

# 4) Time node generation per subject at several levels.
cat("\n== Node generation cost ==\n")
for (level in 3:6) {
  t_nodes <- system.time({
    nodes <- lapply(seq_len(N), function(i) {
      makeSubjectNodes(eta_modes_demo[i, ], H_per_subject[[i]], level)
    })
  })
  n_nodes_per_subj <- nrow(nodes[[1]]$etaNodes)
  cat(sprintf("  level=%d : %d nodes/subject, %.4f s for all %d subjects\n",
              level, n_nodes_per_subj, t_nodes["elapsed"], N))
}

# 5) Time naive per-node joint evaluation for one subject at level 5 (31 nodes).
# This is the Phase-2 batching gap: each call re-routes through the full normL2
# machinery. Phase 2's eval_condition lift reduces this overhead.
cat("\n== Naive per-node joint() evaluation ==\n")
level_bench <- 5L
qn1 <- makeSubjectNodes(eta_modes_demo[1, ], H_per_subject[[1]], level_bench)
B   <- nrow(qn1$etaNodes)

eta_names_full <- as.vector(om$subjectEtas)
chol_pars      <- om$cholPars

# Zero etas for all "other" subjects, fixed init for chol pars.
other_etas <- setNames(rep(0, length(eta_names_full) - 3L),
                       setdiff(eta_names_full, om$subjectEtas[1, ]))

t_eval <- system.time({
  per_node <- lapply(seq_len(B), function(b) {
    eta_subj <- setNames(qn1$etaNodes[b, ], om$subjectEtas[1, ])
    full <- c(init[c("Ka_pop","V_pop","Cl_pop")],
              init[chol_pars],
              eta_subj, other_etas)
    joint(full, deriv = TRUE, conditions = subjects[1])
  })
})
cat(sprintf("  level=%d B=%d for subject 1: %.4f s (%.2f ms/node)\n",
            level_bench, B, t_eval["elapsed"],
            1000 * t_eval["elapsed"] / B))
cat(sprintf("  extrapolated for all %d subjects: %.2f s\n",
            N, t_eval["elapsed"] * N))

# 6) Summary line for the plan's verification gate.
total_floor <- as.numeric(t_eval["elapsed"]) * N
cat(sprintf("\n>>> Phase-1 cost floor for K=3 level=%d: %.1f s per CM-1 eval\n",
            level_bench, total_floor))
cat(sprintf(">>> Plan target: <60s; status: %s\n",
            ifelse(total_floor < 60, "OK", "raise Phase 1b priority")))
