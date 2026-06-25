# Worked examples for symmetryDetection().
#
# Structural identifiability and the symmetries behind it, from the canonical
# scaling to an enzyme assay whose non-identifiabilities are not visible by
# inspection, plus the multi-condition and steady-state machinery. The narrative
# derivations are in vignette("symmetryDetection").

library(dMod)


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

symmetryDetection(reactions, g, method = "observability", reduceCQ = TRUE, closedForm = TRUE)
symmetryDetection(reactions, g, method = "liesym", ansatz = "uni", pMax = 1,
                  reduceCQ = FALSE)
symmetryDetection(reactions, g, method = "scaling", reduceCQ = FALSE)


## 2. A closed-form, non-monomial direction ----------------------------------
#
# With closedForm = TRUE every non-identifiable direction is returned as an exact
# rational function of the unknowns, reconstructed over a finite field. The same
# A <-> B model has, besides the calibration scaling, a conserved-quantity
# direction whose B-component is the rational function -(A + B)/k2.

res <- symmetryDetection(reactions, g, method = "observability",
                         closedForm = TRUE)
lapply(res$nonIdentifiable, function(d) d$vector)


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
lapply(enz$nonIdentifiable, function(d) d$vector)   # kcat*Etot and the s-units
symmetryDetection(f, g.enz, method = "scaling")     # the same two, instantly


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
curve$nonIdentifiable[[1]]$vector


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
c(identifiable = multi$identifiable, conditions = multi$conditions)


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
symmetryDetection(fss, gss, method = "observability",
                  trafo = eqnvec(x = "b/a"), closedForm = TRUE)

# the same verdict with no closed form supplied (the modular solver recovers b/a)
symmetryDetection(fss, gss, method = "observability",
                  equilibrate = TRUE, closedForm = TRUE)

# a known dose pins the scale and makes (a, b, s) identifiable
dose <- addEvent(eventlist(), var = "x", time = 0, value = "dose",
                 method = "replace")
symmetryDetection(fss, gss, method = "observability", events = dose,
                  conditions = data.frame(dose = 2, row.names = "stim"))$identifiable


## 7. A real-world signalling network at scale -------------------------------
#
# A two-receptor TGF-beta / SMAD signalling model: receptor mRNA transcription
# with Hill-type transcriptional feedback (a free exponent nhill on each feedback
# species), translation and degradation, ligand binding to each receptor, the
# active signalling complex, SMAD2/3 phosphorylation, the trimeric SMAD complex
# with SMAD4, and three feedback genes, measured through scaled phospho- and
# total-SMAD, the SMAD4 co-IP and the receptor mRNAs. The perturbations toggle
# transcription (ActD), translation (CHX) and the proteasome (MG132), or knock a
# receptor down, set per condition by events and the condition grid. The model
# has 26 dynamic states and three conserved SMAD moieties, far beyond a symbolic
# determining system. With equilibrate = TRUE the pre-stimulus resting state is
# imposed exactly: forcings held at 0, no events, the steady state solved per
# prime and the TGFb dose then applied as the t0 event. The free Hill exponents
# make the transcription rates non-rational; under equilibrate they are recast to
# a rational coordinate (the log handled as a generic coordinate), so the exact
# finite-field analysis still applies and runs in under a minute.

addRC <- function(eq, from, to, rate, ...) addReaction(eq, from, to, rate, compartment = "Cell", ...)
addRE <- function(eq, from, to, rate, ...) addReaction(eq, from, to, rate, compartment = "extraCell", ...)

reactions <- eqnlist() |>
  addRC("", "bool_ActD",  "0") |> addRC("", "bool_CHX", "0") |> addRC("", "bool_MG132", "0") |>
  addRE("", "TGFb", "0") |>
  addRC("", "R1mRNA",  "k_pr_R1mRNA * (1 + k_inh_R1mRNA_FB3 * FB3^nhill_R1) * (1 - bool_ActD)") |>
  addRC("", "R2mRNA",  "k_pr_R2mRNA / (1 + k_inh_R2mRNA_FB4 * FB4^nhill_R2) * (1 - bool_ActD)") |>
  addRC("", "FB2mRNA", "k_pr_FB2mRNA * C3 * (1 - bool_ActD)") |>
  addRC("", "FB3mRNA", "k_pr_FB3mRNA * C3 * (1 - bool_ActD)") |>
  addRC("", "FB4mRNA", "k_pr_FB4mRNA * C3 * (1 - bool_ActD)") |>
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

# log10 readouts equal scale * h in identifiability content, supplied rationally
observables <- eqnvec(
  pSmad2_obs = "scale_pSmad2 * (pSmad2 + C3)", pSmad3_obs = "scale_pSmad3 * (pSmad3 + C3)",
  TSmad2_obs = "scale_TSmad2 * (Smad2 + pSmad2 + C3)", TSmad3_obs = "scale_TSmad3 * (Smad3 + pSmad3 + C3)",
  Smad4_CoIP_obs = "scale_CoIP * C3", TGFBR1_mRNA_obs = "scale_R1mRNA * R1mRNA",
  TGFBR2_mRNA_obs = "scale_R2mRNA * R2mRNA")

# switches at t = -30 (pre-stimulus regime), TGFb dose at t = 0 (the stimulus)
events <- eventlist() |>
  addEvent(var = "TGFb",       time = 0,   value = "init_TGFb",      method = "replace") |>
  addEvent(var = "bool_CHX",   time = -30, value = "var_bool_CHX",   method = "replace") |>
  addEvent(var = "bool_MG132", time = -30, value = "var_bool_MG132", method = "replace") |>
  addEvent(var = "bool_ActD",  time = -30, value = "var_bool_ActD",  method = "replace")

# one condition per perturbation; receptor knockdowns reparametrise a synthesis rate
cond.grid <- data.frame(Pertubation = c("Ctrl", "ActD", "CHX", "MG132", "R1Knd", "R2Knd"),
                        init_TGFb = 1, stringsAsFactors = FALSE)
cond.grid$var_bool_ActD  <- ifelse(cond.grid$Pertubation == "ActD",  1, 0)
cond.grid$var_bool_CHX   <- ifelse(cond.grid$Pertubation == "CHX",   1, 0)
cond.grid$var_bool_MG132 <- ifelse(cond.grid$Pertubation == "MG132", 1, 0)
cond.grid$k_pr_R1mRNA <- ifelse(cond.grid$Pertubation == "R1Knd", "k_pr_R1mRNA_R1Knd", "k_pr_R1mRNA")
cond.grid$k_pr_R2mRNA <- ifelse(cond.grid$Pertubation == "R2Knd", "k_pr_R2mRNA_R2Knd", "k_pr_R2mRNA")
rownames(cond.grid) <- cond.grid$Pertubation

# switches and ligand are forcings; equilibrate solves the resting state with
# them at 0 (the ligand-receptor complexes fall out as zero automatically) and
# then applies the TGFb dose as the t0 event; cores parallelises the conditions
smad <- symmetryDetection(
  reactions, observables, method = "observability",
  events = events, conditions = cond.grid,
  forcings = c("bool_ActD", "bool_CHX", "bool_MG132", "TGFb"),
  t0 = 0, equilibrate = TRUE, reduceCQ = TRUE, closedForm = TRUE, cores = 6)
# 
c(identifiable = smad$identifiable, rank = smad$rank, dim = smad$dim,
  conditions = smad$conditions)
lapply(smad$nonIdentifiable, function(d) d$support)   # the involved unknowns
# rank 56 / 65: 9 structural non-identifiabilities, all confounded rate/scale
# groups (feedback strengths with their target's production rate, receptor
# synthesis with its readout scale, the SMAD totals with the observation scales);
# both Hill exponents nhill_R1, nhill_R2 are identifiable
