# Open problem: 2 SMAD directions never get a closed form

> **Update / superseded.** The "gauge problem" hypothesis below turned out to be only
> partly right. A log-coordinate decoupling gauge removes the normalization artifact
> (9 → 7 coupled variables) but does **not** close the directions: the residual 7-variable
> coupling is intrinsic — it is the steady-state log-sensitivity of the feedback species
> through its loop, created by `equilibrate` **eliminating** the feedback states. The real
> fix is to keep the states (solve implicitly), so the directions stay low-degree. See
> **`notes/symdetection_implicit_steadystate_plan.md`** for the full diagnosis and the
> implicit-determining-system plan. The original notes below are kept for context.

## Symptom

In `inst/examples/symmetryDetection.R` example #10 (SMAD network, 26 states, 6
conditions, with the added `R1_obs`/`R2_obs` observables), `closedForm = TRUE`
leaves **2 directions** as `type = "general, closed form not found"`:

- `k_inh_R1mRNA_FB3, k_pr_FB3`
- `k_inh_R2mRNA_FB4, k_pr_FB4`

Goal: these two should also get a closed-form vector.

Relevant code (`R/symmetryDetection.R`): `.sym_interpolate_direction`,
`.sym_general_rational_entry`, `.sym_canon_gauge`, the relevance probe (`relProbe`,
~line 1269), and `reconstControl()` caps.

## What is already validated / ruled out (do NOT re-test)

- **Not the variable count.** A synthetic harness (feed a known rational into
  `.sym_general_rational_entry` via `.sym_eval_modq`) reconstructs 9-variable
  rationals, even a degree-9 product, in ~3 s. The Ben-Or-Tiwari path is O(terms),
  independent of variable count.
- **Not the degree cap alone.** `generalDegNum = 12` does not close them.
- **Not pivot-shift over-reporting.** Instrumented `relProbe`: only **1/71** leaves
  is structurally pivot-unstable at `probeRetries = 40`. The 9-leaf coupling is
  genuine, not a probe artifact.
- **Dense fit works but is unusable.** `relevanceCap >= 9` routes to the dense fit,
  which runs 30+ min. Not a path.

## The entry structure

The failing entry couples **9 leaves**:
`k_act_R1_R2, k_deact_R1_R2, k_dg_R1mRNA, k_dg_R2, k_inh_R1mRNA_FB3, k_pr_FB3,
k_pr_R1mRNA, k_pr_R2mRNA, nhill_R1`.

Physically genuine: only 3 (`k_inh_R1mRNA_FB3, k_pr_FB3, nhill_R1`). The other 6
enter through the **feedback loop**: FB3 inhibits R1mRNA, but
R1mRNA -> R1 -> R1_R2 complex -> pSmad -> C3 -> FB3mRNA -> FB3, so FB3 affects its
own production.

## Leading hypothesis: it is a gauge problem

The direction is physically a "weighted scaling" (`k_inh ~ lambda^(-nhill)`,
`k_pr_FB3 ~ lambda`, exactly like r2's `xi_dp = (q-1) dp`). But the free-column /
canonical gauge produces a genuinely 9-variable, high-degree representative of it.
The simple form lives in a different gauge the algorithm does not select.
`.sym_canon_gauge` only reduces the z-support modulo the scalings; it does not
reduce the per-entry leaf coupling.

## Ideas to brainstorm

1. **Generalized / weighted-scaling detection**: directions of the form
   `xi_i = c_i * z_i` with `c_i` a constant (parameter, not just integer). Then
   `xi_i / z_i = c_i` is constant (the cheap constEntry path) and the 9-variable
   blow-up disappears. Trick: detect via constancy of
   `(xi_i / xi_anchor) * (z_anchor / z_i)` across points.
2. **Better gauge selection**: pick the minimal-leaf-coupling representative, i.e.
   reduce the direction modulo the *other* directions over the rational function
   field (not just modulo the integer scaling lattice as canon_gauge does now).
3. **Verification-guarded leaf reduction**: hold candidate gauge leaves fixed during
   reconstruction; accept only if `.sym_verify_direction` passes at a fresh point
   (sound: a wrong reduction is rejected). Caveat: the free-column section genuinely
   depends on those leaves, so this needs the right gauge first.

## Constraints when working on this

- SMAD runs take ~150-300 s. **Always wrap test runs in a `timeout` watchdog**, and
  iterate on the **fast synthetic harness** (known rationals into
  `.sym_general_rational_entry`), not the full model.
- Do NOT launch several heavy SMAD runs in parallel; that oversubscribed the machine
  and hung for hours.

## Status of the branch

- Diagnostic feature is built: on a failed closed form, `summary()` prints the
  reason and which `reconstControl` cap was hit; `print()` notes the count.
- `symmetryDetection` S3 class with lightweight `print` + detailed `summary`.
- Inverted recast (no-bare-handle power laws) done earlier this session.
- 172 tests green. Nothing committed.
