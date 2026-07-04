# Worked examples for symmetryDetection().
#
# Structural identifiability and the symmetries behind it, from the canonical
# scaling to an enzyme assay whose non-identifiabilities are not visible by
# inspection, plus the multi-condition and steady-state machinery. The narrative
# derivations are in vignette("symmetryDetection").

library(dMod)
setwd(tempdir())

## 1. The canonical scaling symmetry, found by all three engines -------------
#
# A reversible reaction A <-> B observed only through alpha * A leaves the
# absolute scale free: scaling A and B by a factor and dividing alpha by it
# leaves the output unchanged.

reactions <- eqnlist() |>
  addReaction("A", "B", "k1 * A") |>
  addReaction("B", "A", "k2 * B") |>
  customTotals(list(totalAB = "A+B"))

g <- eqnvec(Aobs = "alpha * A")

summary(symmetryDetection(reactions, g, method = "observability",
                          reduceCQ = TRUE, closedForm = TRUE))
summary(symmetryDetection(reactions, g, method = "liesym", reduceCQ = FALSE,
                          liesym = liesymControl(ansatz = "uni", pMax = 1L)))
summary(symmetryDetection(reactions, g, method = "scaling", reduceCQ = FALSE))


## 2. A closed-form, non-monomial direction ----------------------------------
#
# With closedForm = TRUE every non-identifiable direction is returned as an exact
# rational function of the unknowns, reconstructed over a finite field. The same
# A <-> B model has, besides the calibration scaling, a conserved-quantity
# direction whose B-component is the rational function -(A + B)/k2.

res <- symmetryDetection(reactions, g, method = "observability", equilibrate = TRUE,
                         closedForm = TRUE)
summary(res)


## 3. An enzyme assay and its hidden symmetries -------------------------------
#
# Substrate depletion under the quasi-steady-state law, read out on an unknown
# scale s. Two non-obvious non-identifiabilities: the turnover number and the
# enzyme amount enter only as the product Vmax = kcat*Etot, and the molar units
# of S, Km, Etot trade against the readout gain s. Both are scalings, so the
# scaling engine agrees with the observability engine.

f <- eqnvec(S = "-kcat*Etot*S/(Km + S)")
g.enz <- eqnvec(y = "s*S")

enz <- symmetryDetection(f, g.enz, method = "observability", closedForm = TRUE)
summary(enz)                                          # kcat*Etot and the s-units
summary(symmetryDetection(f, g.enz, method = "scaling"))   # the same two


## 4. Two rates confined to a nonlinear curve --------------------------------
#
# Gene expression observed only at the protein level: transcription ktx and
# translation ktl cannot be separated, since more mRNA poorly translated and
# less mRNA well translated give the same protein. The two rates are confined to
# the hyperbola ktx*ktl = const; any point on it yields the same trajectory
# (only the product, setting the steady state, is identifiable). The single
# returned direction (ktx:1, ktl:-ktl/ktx, m:m/ktx) is the tangent of that curve.

gene <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")   # mRNA hidden
curve <- symmetryDetection(gene, eqnvec(y = "p"), method = "observability",
                           closedForm = TRUE)

summary(curve)
## 5. Stacking conditions resolves a switch-gated rate ------------------------
#
# The rate k1 + u*k2 is gated by a switch u. A single switch value cannot
# separate k1 and k2; two values, set per condition by an event with the value
# read from the grid, identify both. A direction counts as non-identifiable only
# if it is unobservable in every condition.

fu <- eqnvec(A = "-(k1 + u*k2) * A", u = "0")
events <- addEvent(eventlist(), var = "u", time = -1, value = "var_u",
                   method = "replace")
cond.grid <- data.frame(var_u = c(0, 1), row.names = c("ctrl", "stim"))

multi <- symmetryDetection(fu, eqnvec(y = "A"), method = "observability",
                           events = events, conditions = cond.grid)

summary(multi)


## 6. Pre-equilibration: a steady-state initial condition ---------------------
#
# Before the experiment the system rests at a steady state, so the initial state
# is not free but solves f = 0. Provide it explicitly as a state-named `trafo`
# entry when it has a closed form, or let `equilibrate = TRUE` solve f = 0
# numerically (per prime, over a finite field) with its implicit-function-theorem
# parameter sensitivities. A known dose then perturbs the resting state.

fss <- eqnvec(x = "b - a*x")    # turnover; resting steady state x* = b/a
gss <- eqnvec(y = "s*x")        # scaled readout

# explicit steady state: only the resting product s*b/a is seen (rank 1 of 3)
summary(symmetryDetection(fss, gss, method = "observability",
                          trafo = eqnvec(x = "b/a"), closedForm = TRUE))

# the same verdict with no closed form supplied (the modular solver recovers b/a)
summary(symmetryDetection(fss, gss, method = "observability",
                          equilibrate = TRUE, closedForm = TRUE))

# a known dose pins the scale and makes (a, b, s) identifiable
dose <- addEvent(eventlist(), var = "x", time = 0, value = "dose",
                 method = "replace")

res <- symmetryDetection(fss, gss, method = "observability", events = dose, 
                         conditions = data.frame(dose = 2, row.names = "stim"))$identifiable

summary(res)
## 7. Free Hill and power exponents -------------------------------------------
#
# A rate with a free exponent (a Hill term, or a bare power base^n with n a
# parameter) is not rational. Under equilibrate = TRUE the term base^n is recast
# to a coordinate E = base^n with a companion L = log(base), so the exact
# finite-field analysis still applies and the exponent is assessed like any other
# parameter. The closed form differs by case. A base with a linear turnover term
# keeps base solved and E generic, so a direction that frees the exponent carries
# an explicit log(base) factor. A base without a turnover term has a
# fractional-root steady state, so the roles invert (E solved, base and log(base)
# generic) and the direction is rational.

# self-inhibition with a free Hill exponent n: n is identifiable, and the
# direction that frees it returns with its log(base) factor
hill <- eqnlist() |>
  addReaction("0",  "FB", "k_pr_FB")                     |>
  addReaction("FB", "0",  "d_FB * FB")                   |>
  addReaction("0",  "x",  "k_pr_x * K^n / (K^n + FB^n)") |>
  addReaction("x",  "0",  "d_x * x")
res <- symmetryDetection(hill, eqnvec(xobs = "scale * x"),
                         method = "observability", equilibrate = TRUE,
                         closedForm = TRUE)

summary(res)

# a bare power with no turnover term: x* = ((kpr + kin*u)/dp)^(1/q) is a
# fractional root, the inverted recast keeps q identifiable, and the direction is
# the rational scaling xi_dp = (q - 1)*dp, xi_s = s (no log)
ev <- addEvent(eventlist(), var = "u", time = -1, value = "1", method = "replace")
symmetryDetection(eqnvec(x = "kpr - dp*x^q + kin*u", u = "0"), eqnvec(y = "s*x"),
                  method = "observability", equilibrate = TRUE, forcings = "u",
                  events = ev, conditions = data.frame(var = 1, row.names = "stim"),
                  closedForm = TRUE) |> summary()


## 8. EGF/EGFR into a MEK/ERK cascade, partially observed --------------------
#
# A single-receptor signalling network: EGF binds EGFR to form the active
# complex EGF_EGFR, which phosphorylates MEK, and pMEK in turn phosphorylates
# ERK; each kinase has a constitutive dephosphatase. Four conserved moieties
# (ligand, receptor, MEK, ERK) are auto-handled. The readout sees only the
# scaled phospho-forms scale_pMEK*pMEK and scale_pERK*pERK, so the ligand, the
# receptor and the unphosphorylated kinase pools are hidden: the system is only
# partially observed. Everything is mass-action (rational right-hand sides, no
# free exponents), so every non-identifiable direction is recovered in closed
# form, including the one non-monomial direction. The SMAD model below adds a
# feedback through a Hill term with a free exponent nhill; that direction is a
# weighted scaling whose weight is the exponent itself (xi_kinh = -nhill * kinh),
# recovered in closed form via its minimal support in log coordinates.

reactions <- eqnlist() |>
  addReaction("EGF + EGFR", "EGF_EGFR", "k_bind * EGF * EGFR") |>
  addReaction("EGF_EGFR", "EGF + EGFR", "k_unbind * EGF_EGFR") |>
  addReaction("MEK", "pMEK", "k_phos_MEK * EGF_EGFR * MEK") |>
  addReaction("pMEK", "MEK", "k_dephos_MEK * pMEK") |>
  addReaction("ERK", "pERK", "k_phos_ERK * pMEK * ERK") |>
  addReaction("pERK", "ERK", "k_dephos_ERK * pERK")
reactions <- customTotals(reactions, list(
  totalEGF  = "EGF + EGF_EGFR", totalEGFR = "EGFR + EGF_EGFR",
  totalMEK  = "MEK + pMEK",     totalERK  = "ERK + pERK"))

# only the scaled phospho-forms are measured
observables <- eqnvec(pMEK_obs = "scale_pMEK * pMEK",
                      pERK_obs = "scale_pERK * pERK")

egf <- symmetryDetection(reactions, observables, method = "observability",
                         reduceCQ = FALSE, closedForm = TRUE)
summary(egf)
# rank 11 / 15: two readout-gain scalings (each scale trades against its kinase
# total and the downstream rate), one receptor-pool scaling, and one non-monomial
# direction in the EGF/EGFR binding constants, all returned as exact rational
# functions.


## 9. Multiple stimuli: analytic segments and exact gap propagation ----------
#
# The Taylor expansion starts at the earliest event time; events at later times
# split the timeline at the distinct event times into analytic segments, one local
# jet each, all stacked (like conditions, but over time). The state is propagated
# exactly across each inter-event gap by keeping the gap length a formal
# power-series variable: an identifiability direction must annihilate every
# gap-order coefficient (a polynomial that vanishes for generic timing), and the
# gap order is raised until the rank saturates (gapOrderUsed). A protocol whose
# events all sit at the earliest time stays a single segment.

# A receptor R driven by a step forcing u, the SMAD pattern below in miniature.
# With the stimulus applied at t = 0 (the only event) this is a single segment.
f <- eqnvec(R = "kpr - kdg*R + kon*u*R", u = "0")
gR <- eqnvec(y = "scale*R")
ev <- addEvent(eventlist(), var = "u", time = 0, value = "init_u", method = "replace")
cg <- data.frame(init_u = 1, row.names = "Ctrl")
one <- symmetryDetection(f, gR, method = "observability", equilibrate = TRUE,
                         events = ev, conditions = cg, forcings = "u")

summary(one)

# a later event (a washout of u at t = 60) opens a second segment whose relaxation
# is propagated exactly across the gap and stacked on top
ev2 <- ev |> addEvent(var = "u", time = 60, value = "0", method = "replace")
two <- symmetryDetection(f, gR, method = "observability", equilibrate = TRUE,
                         events = ev2, conditions = cg, forcings = "u")

summary(two)

# Exact propagation identifies a transient-channel parameter. The inhibition
# strength kinh affects only the inhibitor-relaxed state at the stimulus onset.
# The inhibitor is applied at t = -30 (the earliest event) and the state is
# propagated across the gap to the stimulus at t = 0, so kinh becomes identifiable
# through that transient channel (gapOrderUsed >= 1).
f2 <- eqnvec(x = "kpr/(1 + kinh*inh) - kdeg*x + kstim*stim", inh = "0", stim = "0")
g2 <- eqnvec(y = "s*x")
ev3 <- eventlist() |>
  addEvent(var = "inh",  time = -30, value = "1", method = "replace") |>
  addEvent(var = "stim", time = 0,   value = "1", method = "replace")
rmes30 <- symmetryDetection(f2, g2, method = "observability", equilibrate = TRUE,
                            events = ev3, forcings = c("inh", "stim"))

summary(rmes30)


## 10. A real-world signalling network at scale ------------------------------
#
# A two-receptor TGF-beta / SMAD signalling model: receptor mRNA transcription
# with Hill-type transcriptional feedback (a free exponent nhill on each feedback
# species), translation and degradation, ligand binding to each receptor, the
# active signalling complex, SMAD2/3 phosphorylation, the trimeric SMAD complex
# with SMAD4, and three feedback genes, measured through scaled phospho- and
# total-SMAD, the SMAD4 co-IP and the receptor mRNAs. The perturbations toggle
# transcription (ActD), translation (CHX) and the proteasome (MG132), or knock a
# receptor down, set per condition by events and a per-condition trafo list. The
# model has 26 dynamic states and three conserved SMAD moieties, far beyond a
# symbolic determining system. The pre-stimulus resting state is imposed with
# equilibrate = TRUE: the states stay coordinates and f = 0 enters as a tangency
# constraint solved per prime, the TGFb dose applied on top. (The explicit route --
# feeding a steadyStates() solution through the trafo -- needs a rational steady
# state; SMAD's Hill feedback makes it non-rational, so equilibrate is the route
# here.) The free Hill exponents make the transcription rates non-rational; they are
# recast to a rational coordinate (the log handled as a generic coordinate), so the
# finite-field analysis still applies and runs in a few minutes.

addRC <- function(eq, from, to, rate, ...) addReaction(eq, from, to, rate, compartment = "Cell", ...)
addRE <- function(eq, from, to, rate, ...) addReaction(eq, from, to, rate, compartment = "extraCell", ...)

reactions <- eqnlist() |>
  addRC("", "bool_ActD",  "0") |> addRC("", "bool_CHX", "0") |> addRC("", "bool_MG132", "0") |>
  addRE("", "TGFb", "0") |>
  addRC("", "R1mRNA",  "k_pr_R1mRNA * (1 + k_inh_R1mRNA_FB3 * FB3^nhill_R1) * (1 - bool_ActD)") |>
  addRC("", "R2mRNA",  "k_pr_R2mRNA / (1 + k_inh_R2mRNA_FB4 * FB4^nhill_R2) * (1 - bool_ActD)") |>
  addRC("", "FB2mRNA", "k_pr_FB2mRNA * C3^nhill_FB2mRNA / (Km_FB2mRNA^nhill_FB2mRNA + C3^nhill_FB2mRNA) * (1 - bool_ActD)") |>
  addRC("", "FB3mRNA", "k_pr_FB3mRNA * C3^nhill_FB3mRNA / (Km_FB3mRNA^nhill_FB3mRNA + C3^nhill_FB3mRNA) * (1 - bool_ActD)") |>
  addRC("", "FB4mRNA", "k_pr_FB4mRNA * C3^nhill_FB4mRNA / (Km_FB4mRNA^nhill_FB4mRNA + C3^nhill_FB4mRNA) * (1 - bool_ActD)") |>
  addRC("R1mRNA", "", "k_dg_R1mRNA * R1mRNA") |> addRC("R2mRNA", "", "k_dg_R2mRNA * R2mRNA") |>
  addRC("FB2mRNA", "", "k_dg_FB2 * FB2mRNA") |> addRC("FB3mRNA", "", "k_dg_FB3 * FB3mRNA") |>
  addRC("FB4mRNA", "", "k_dg_FB4 * FB4mRNA") |>
  addRC("", "R1",  "k_pr_R1 * R1mRNA * (1 - bool_CHX)") |>
  addRC("", "R2",  "k_pr_R2 * R2mRNA * (1 - bool_CHX)") |>
  addRC("", "FB2", "k_pr_FB2 * FB2mRNA * (1 - bool_CHX)") |>
  addRC("", "FB3", "k_pr_FB3 * FB3mRNA * (1 - bool_CHX)") |>
  addRC("", "FB4", "k_pr_FB4 * FB4mRNA * (1 - bool_CHX)") |>
  addRC("R1", "", "k_dg_R1 * R1 * (1 - bool_MG132)") |> addRC("R2", "", "k_dg_R2 * R2 * (1 - bool_MG132)") |>
  addRC("FB2", "", "k_dg_FB2 * FB2 * (1 - bool_MG132)") |> addRC("FB3", "", "k_dg_FB3 * FB3 * (1 - bool_MG132)") |>
  addRC("FB4", "", "k_dg_FB4 * FB4 * (1 - bool_MG132)") |>
  addRC("R1 + R2", "R1_R2", "k_act_R1_R2 * R1 * R2") |>
  addRC("R1_R2", "R1 + R2", "k_deact_R1_R2 * R1_R2") |>
  addRC("R1_R2", "", "k_dg_R1_R2 * R1_R2 * (1 - bool_MG132)") |>
  addRC("R2 + TGFb", "R2_TGFb", "k_act_R2_TGFb * TGFb * R2 / (km_R2 + R2 + TGFb)", rateCompartment = "Cell") |>
  addRC("R2_TGFb", "R2 + TGFb", "k_deact_R2_TGFb * R2_TGFb") |>
  addRC("R2_TGFb", "R2_TGFb_int", "k_int_R2_TGFb * R2_TGFb") |>
  addRC("R2_TGFb_int", "TGFb", "k_decay_R2_TGFb_int * R2_TGFb_int") |>
  addRC("R2_TGFb_int", "", "k_dg_R2_TGFb_int * R2_TGFb_int * (1 - bool_MG132)") |>
  addRC("R1 + TGFb", "R1_TGFb", "k_act_R1_TGFb * TGFb * R1 / (km_R1 + R1 + TGFb)", rateCompartment = "Cell") |>
  addRC("R1_TGFb", "R1 + TGFb", "k_deact_R1_TGFb * R1_TGFb") |>
  addRC("R1_TGFb", "R1_TGFb_int", "k_int_R1_TGFb * R1_TGFb") |>
  addRC("R1_TGFb_int", "TGFb", "k_decay_R1_TGFb_int * R1_TGFb_int") |>
  addRC("R1_TGFb_int", "", "k_dg_R1_TGFb_int * R1_TGFb_int * (1 - bool_MG132)") |>
  addRC("R1 + R2_TGFb", "R1_R2_TGFb", "k_act_R1_R2_TGFb * R1 * R2_TGFb") |>
  addRC("R2 + R1_TGFb", "R1_R2_TGFb", "k_act_R2_R1_TGFb * R2 * R1_TGFb") |>
  addRC("R1_R2 + TGFb", "R1_R2_TGFb", "k_act_R1_R2_TGFb_direct * TGFb * R1_R2 / (km_R1_R2 + R1_R2 + TGFb)", rateCompartment = "Cell") |>
  addRC("R1_R2_TGFb", "", "(k_dg_R1_R2_TGFb + k_dg_R1_R2_TGFb_FB1 * FB2) * R1_R2_TGFb * (1 - bool_MG132)") |>
  addRC("Smad2", "pSmad2", "(k_phospho_pS2 / (1 + k_inh_pSmad2_FB2 * FB2)) * Smad2 * R1_R2") |>
  addRC("Smad2", "pSmad2", "(k_phospho_pS2 / (1 + k_inh_pSmad2_FB1 * FB2)) * Smad2 * R1_R2_TGFb") |>
  addRC("pSmad2", "Smad2", "k_dephos_S2 * pSmad2") |>
  addRC("Smad3", "pSmad3", "(k_phospho_pS3 / (1 + k_inh_pSmad3_FB2 * FB2)) * Smad3 * R1_R2") |>
  addRC("Smad3", "pSmad3", "(k_phospho_pS3 / (1 + k_inh_pSmad3_FB2 * FB2)) * Smad3 * R1_R2_TGFb") |>
  addRC("pSmad3", "Smad3", "k_dephos_S3 * pSmad3") |>
  addRC("pSmad2 + pSmad3 + Smad4", "C3", "k_form_S4Coip * pSmad2 * pSmad3 * Smad4") |>
  addRC("C3", "Smad2 + pSmad3 + Smad4", "k_dissolve_C3_dp2 * C3") |>
  addRC("C3", "pSmad2 + Smad3 + Smad4", "k_dissolve_C3_dp3 * C3")
reactions$compartments$Cell$volume      <- "1"
reactions$compartments$extraCell$volume <- "volumeEC"
reactions <- customTotals(reactions, list(totalSMAD2 = "Smad2 + pSmad2 + C3",
                                          totalSMAD3 = "Smad3 + pSmad3 + C3",
                                          totalSMAD4 = "Smad4 + C3"))

observables <- eqnvec(
  R1_obs = "scale_R1 * R1", R2_obs = "scale_R2 * R2", 
  pSmad2_obs = "scale_pSmad2 * (pSmad2 + C3)", pSmad3_obs = "scale_pSmad3 * (pSmad3 + C3)",
  TSmad2_obs = "scale_TSmad2 * (Smad2 + pSmad2 + C3)", TSmad3_obs = "scale_TSmad3 * (Smad3 + pSmad3 + C3)",
  Smad4_CoIP_obs = "scale_CoIP * C3", TGFBR1_mRNA_obs = "scale_R1mRNA * R1mRNA",
  TGFBR2_mRNA_obs = "scale_R2mRNA * R2mRNA",
  TGFb_obs = "scale_TGFb * TGFb")

# switches at t = -30 (pre-stimulus regime), TGFb dose at t = 0 (the stimulus)
events <- eventlist() |>
  addEvent(var = "TGFb",       time = 0,   value = "init_TGFb",      method = "replace") |>
  addEvent(var = "bool_CHX",   time = -30, value = "var_bool_CHX",   method = "replace") |>
  addEvent(var = "bool_MG132", time = -30, value = "var_bool_MG132", method = "replace") |>
  addEvent(var = "bool_ActD",  time = -30, value = "var_bool_ActD",  method = "replace")

# one condition per perturbation, as a per-condition covariate grid: each row bakes
# the ActD/CHX/MG132 switch values and the TGFb dose (the event-value placeholders),
# and the receptor knockdowns rename a synthesis rate. The perturbation name is the
# condition (row) name, so it is dropped as a column before branching.
cond.grid <- data.frame(Pertubation = c("Ctrl", "ActD", "CHX", "MG132", "R1Knd", "R2Knd"),
                        init_TGFb = 1, stringsAsFactors = FALSE)
cond.grid$var_bool_ActD  <- ifelse(cond.grid$Pertubation == "ActD",  1, 0)
cond.grid$var_bool_CHX   <- ifelse(cond.grid$Pertubation == "CHX",   1, 0)
cond.grid$var_bool_MG132 <- ifelse(cond.grid$Pertubation == "MG132", 1, 0)
cond.grid$k_pr_R1mRNA <- ifelse(cond.grid$Pertubation == "R1Knd", "k_pr_R1mRNA_R1Knd", "k_pr_R1mRNA")
cond.grid$k_pr_R2mRNA <- ifelse(cond.grid$Pertubation == "R2Knd", "k_pr_R2mRNA_R2Knd", "k_pr_R2mRNA")
rownames(cond.grid) <- cond.grid$Pertubation
cond.grid$Pertubation <- NULL

# The conditions enter through a per-condition trafo list built with branch(): it
# broadcasts a base trafo over the grid and define()s each row's columns into that
# condition's copy (a numeric bake for the switches/dose, a rename for a knockdown)

cond.trafo <- eqnvec() |> 
  define("x~x", x = getParameters(reactions, events)) |> 
  branch(table = cond.grid, apply = "insert")

# switches and ligand are forcings. equilibrate imposes the resting state implicitly;
# the TGFb dose is applied on top, and the perturbations enter through the trafo list.
# In the equilibrate path each sample point's steady state is solved per condition, so
# cores.GLp (per-point parallelism) is the knob that matters here, not cores.conditions;
# raise it towards the number of physical cores to speed the run up.
outsmad <- symmetryDetection(
  reactions, observables, method = "observability",
  events = events, trafo = cond.trafo,
  forcings = c("bool_ActD","bool_CHX","bool_MG132","TGFb"),
  equilibrate = TRUE, reduceCQ = TRUE, closedForm = TRUE,
  cores.GLp = 8)

# 9 of the 10 non-identifiabilities close in closed form (a few minutes):
#   - the mRNA/protein synthesis scalings (each k_pr_X paired with its readout scale),
#   - the pSmad-feedback inhibition group,
#   - both receptor Hill feedbacks, each a parameter-weighted scaling recovered exactly
#     over Q(nhill) by the toric peel: xi_kinh = -nhill * kinh (the inhibition strength
#     against the feedback protein's synthesis rate; the Hill exponents are identifiable),
#   - the ligand initial value init_TGFb.
# The 10th direction stays open (reported by support only): the SMAD-pool confound
# couples the conserved totals and observation scales with the complex-formation rate
# AND the feedback-mRNA Hill exponents nhill_FB* -- a scaling of the pool that also
# shifts those exponents is transcendental (exp-sensitivity ~ base^n * log base), so it
# has no rational closed form. Its sample bank cannot fill (perturbing the coupled
# parameters leaves the recast manifold), and the reconstruction bails to support-only.
summary(outsmad)
