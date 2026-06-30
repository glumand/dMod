# Plan: implicit steady-state determining system for closed-form symmetry directions

Status: **IMPLEMENTED, TESTED ON SMAD, THEN REVERTED (2026-06-30, user decision).** The
full joint path worked on small models and ran on SMAD (192s, closed 12/14) but the 2 target
Hill directions still did NOT close: they are `-nhill`-WEIGHTED scalings (the free Hill
exponent is the weight), not integer scalings, so the peel misses them and the rational
reconstruction re-entangles with the upstream chain (the relevance probe re-solves the resting
steady state, so state values still depend on the chain rates — keeping states as coordinates
did not decouple the reconstruction). Code reverted to the explicit/eliminated version; SMAD's
blocker has "some other cause" to brainstorm. The implementation below is recoverable from the
2026-06-30 session transcript if revisited. See `memory/project_symdetection_implicit_ss.md`
for the failure analysis and the open threads.

### Progress

- Kernel: `symObsNullMulti`/`symObsNullChain` take `nzExtra` (decouples `nz` from
  `zSlots.size()`); icSeed-backed state log-columns. Backward-compatible, recompiled.
- Python `solveSteadyStateModular`: returns `dfRows`/`dfStateCols`/`dfParamCols`
  (resting Jacobian `[Jx|Jt]`) and `xrest` (pre-t0-event resting values).
- Python `compileObservabilityTapeMulti`: `jointSteadyState` flag appends real states
  to `znames` as extra log z-columns; emits `nzExtra`/`stateZNames`.
- R `symmetryDetection(implicitSteadyState=TRUE)` → `.observability_analytic_multi`
  joint branch: per-condition icSeed with value `x*` (post-dose) and log dual; `df_rest`
  rows mapped to mixed coords (`.sym_joint_df_mat`, params additive, state cols `·x*_rest`
  via modular multiply — naive `*` overflows R doubles) and stacked via
  `.sym_stack_df_rows`; peel uses `zval=1` on the state log-columns.
- **Validated:** `f=b-a*x, g=s*x`, equilibrate, multiply-dose at t0. Eliminated →
  rank 2/3, scaling `{b:1, s:-1}`. Joint → rank 3/4, scaling `{b:1, s:-1, x:1}`: the
  parameter projection matches the eliminated path, plus the state weight. The
  no-stimulus model is the known-degenerate case (constant trajectory); use a stimulated
  model for joint validation.

### Remaining

1. Recast/Hill for SMAD: keep base FB as a state log-column, handle E/L; `df_rest` in
   recast coords. The Python `dfRows` is built from `solveStates` (post-recast) — verify
   the column set and the log multiplication for the recast base.
2. Reporting polish: optional split of the state weights into a `stateVector` field.
3. A unit test for the joint path; then the SMAD acceptance run.

## TL;DR

The two SMAD directions that never get a closed form
(`{k_inh_R1mRNA_FB3, k_pr_FB3}`, `{k_inh_R2mRNA_FB4, k_pr_FB4}`) are **not** a gauge
artifact. The root cause is that `equilibrate = TRUE` **eliminates the feedback
species** (solves `f(x,θ)=0` for the states `x*` and substitutes). For a Hill
feedback in a loop, `x*(θ)` is a high-degree determinant ratio, so the direction's
parameter entries become irreducibly high-degree rationals (7–9 variables) that the
interpolation cannot and should not reconstruct.

The fix is to **stop eliminating the states**: set up the determining system in the
**joint `(x, θ)` space** with the steady-state condition kept as a constraint, and
take the nullspace there. The entries are then low-degree, the parameter components
come out compact and state-free for these scaling symmetries (the old output format,
for free), and the transcendental Hill term is handled by the existing recast.

## Problem recap

`inst/examples/symmetryDetection.R` example #10 (SMAD, 26 states, 6 conditions,
`equilibrate = TRUE`, `closedForm = TRUE`). Rank 58/65. Seven scaling directions close;
the two Hill-feedback directions are reported `type = "general, closed form not found"`
with reason "an entry couples 9 variables and the sparse fit hit its caps".

These are weighted scalings whose weight is the Hill exponent: the invariant is
`k_inh_R1mRNA_FB3 · FB3^nhill_R1`. The physically clean generator is
`ξ_{k_inh} = -nhill·k_inh`, `ξ_{k_pr_FB3} = k_pr_FB3`, `ξ_{FB3} = FB3` (and 0 else) —
but in the eliminated parameter-only space it blows up.

## What was tried and ruled out (do NOT redo)

All at the SMAD scale, watchdog'd, serial (~150 s each):

1. **Log-coordinate decoupling gauge** (eta = v/z, decouple, back-substitute
   xi = eta·z). Implemented and sound (verification-guarded, all tests green). It
   removes the *normalization* artifact: **9 → 7 variables**, but does **not** close —
   7 intrinsic variables remain. Sound but ineffective for SMAD.
2. **Scaling reduction in log space** (`reduceScal` modulo the integer scaling span,
   recovering weights by a centred lift). Stays at **7 variables** — only the
   representative changes (anchor moves to `k_pr_FB3mRNA`). The 7 are not a scaling
   admixture.
3. **Raised caps** (`generalDegNum=10, generalDegDen=8, termCap=200`). Still fails —
   the entry is genuinely high-degree (total degree > 18 or > 200 terms).
4. **Relevant-set dump** (decisive): the failing entry's relevant variables are
   `{k_act_R1_R2, k_deact_R1_R2, k_dg_R1mRNA, k_dg_R2, k_pr_R1mRNA, k_pr_R2mRNA,
   nhill_R1}` — i.e. the rates that set FB3's steady state through its feedback loop.
   **`E_FB3` / `L_FB3` (the recast coords) are NOT relevant.** So the chain coupling
   is *not* a back-substitution artifact; the augmented direction genuinely is a
   high-degree rational of the chain rates.
5. **Recast inversion check** (`solveSteadyStateModular`, lines ~1677-1688): FB3 has a
   linear turnover term (`k_dg_FB3·FB3`), so it is in the **non-inverted** branch →
   **FB3 (base) is solved/eliminated**, `E = FB3^nhill` is held generic. The
   elimination of FB3 is the source of the coupling. The "inverted" switch cannot be
   repurposed (it requires the base to have *no* turnover term and solves `E` instead),
   and keeping only FB3 generic merely relocates the elimination one step upstream
   (to FB3mRNA → C3 → …). The coupling is the whole loop.

Conclusion: the elimination of the feedback-loop states is the cause; no gauge in
parameter space and no stronger interpolation removes it. The compact form lives in
the un-eliminated `(x, θ)` space.

## The approach: implicit / joint determining system

### Math

The parameter non-identifiability direction is the projection of a generator tangent
to the steady-state manifold that preserves the output jet. Let `x0` be the (free)
initial-condition states and `θ` the parameters. The symmetry `ξ = (ξ_x, ξ_θ)` must
satisfy:

- **output-preserving:** `Obs(x0,θ) · ξ = 0`, where `Obs = ∂(output jet)/∂(x0, θ)` is
  the standard observability/identifiability matrix with states as free ICs (this is
  exactly what the **non-equilibrate** observability mode already builds).
- **tangent to the steady-state manifold:** `df(x0,θ) · ξ = 0`, where
  `df = ∂f/∂(x0, θ)` is the (low-degree) RHS Jacobian. This encodes
  `ξ_x = -(∂f/∂x)^{-1}(∂f/∂θ) ξ_θ` *without ever forming the inverse* — the inverse is
  exactly the high-degree determinant ratio we must avoid.

So `ξ ∈ ker([Obs ; df])` in the joint `(x0, θ)` coordinate space. The entries are
low-degree (bounded by the model RHS degree, ~2–3 for mass action), because nothing is
eliminated.

### Sampling must be ON the resting manifold (corrected 2026-06-30)

The earlier draft claimed `df·ξ = 0` holds at *generic* `(x0,θ)` so we could sample
off-manifold cheaply. **That is wrong** and was disproven at the standalone harness
(`scratchpad/implicit_math_check.py`, `implicit_onmanifold.py`, `implicit_loop_fast.py`):

- For the clean generator, `df_rest·ξ = 1·f_rest` **row-by-row** — it equals the RHS
  itself, so it vanishes **only on** `{f_rest=0}`, not generically. A dynamical
  symmetry maps the manifold to itself, hence ξ is *tangent* to `{f=0}` **on** the
  manifold, but not off it. Off-manifold `ker([Obs;df])` is a different, high-degree
  space that does **not** contain the clean direction.
- **Therefore samples must lie on the resting steady-state manifold**:
  `solveSteadyStateModular` (solve `f_rest=0` per (point,prime)) stays essential. We do
  **not** avoid solving `f=0`; we avoid forming `x*(θ)` *symbolically* / folding its IFT
  duals into the parameter columns.

### Resting field ≠ observed field (new requirement)

The constraint and the observation use **different** vector fields:

- `df` (constraint rows, tangency to the IC manifold) uses **`f_rest`** = the RHS with
  all forcings/stimuli **off** (the resting steady state the IC is drawn from).
- `Obs` (output jet) uses **`f_stim`** = the RHS **with** forcings/events on (the
  observed transient). For the clean scaling, `Obs_stim·ξ ≡ 0` holds *identically*
  (off-manifold too), because the scaling is a Lie symmetry of the stimulated field and
  the observed species are invariant.

So the determining system is `ξ ∈ ker([ Obs(f_stim) ; df(f_rest) ])` **evaluated at
resting-manifold points** `x = x*_rest(θ)`.

### Per-condition steady states ⇒ states need LOG z-columns (corrected 2026-06-30, user)

Each experimental condition has its **own** resting steady state `x*_c` (conditions can
differ in kinetic substitutions, not only in forcings). Consequences for the joint mode:

- A single **additive** shared state z-column does **not** work: a scaling's additive
  nullspace entry on a state would be `w·x*_c` — different per condition — but the kernel
  returns one shared nullspace vector over shared z-columns. Contradiction.
- Fix: make the **state z-columns LOG coordinates**. Seed each condition's state via
  `icSeed` with value `x*_c` and dual `∂state_i/∂z_i = x*_c` (the value itself) on its own
  state column, `0` on the parameter columns. The log dual absorbs `x*_c`, so a scaling
  `ξ_state = w·state` has the **condition-independent** weight `w` on that column — shared
  across conditions, exactly what `.sym_peel_scalings` needs. (Parameters stay additive,
  dual `1`, because their values are already shared.) This is why the earlier log-gauge
  work pointed the right way; it was just applied post-hoc on the eliminated parameters
  instead of being built into the determining system here.
- Mixed-coordinate bookkeeping: the kernel's Obs rows then use `∂jet/∂p` (param cols,
  additive) and `∂jet/∂state · x*_c` (state cols, log). `df_rest,c` rows likewise:
  `[ Jt (param cols, additive) | Jx · x*_c (state cols, log) ]`, stacked per condition.
  Peel uses `zval = param value` for param cols and `zval = 1` for state log-cols.
- `nz = nParams + nStates`. State columns are **icSeed-backed, not leaf-backed**, so
  `zSlots` keeps only the `nParams` leaf slots and the kernel must take `nz` decoupled
  from `zSlots.size()` (small C++ change: pass `nz` explicitly; the `dualCol` loop runs
  over `zSlots.size()`, the extra state columns get their duals solely from `icSeed`).

### The directions are integer-weight scalings over joint coordinates

Confirmed at the loop harness (`implicit_loop_fast.py`, genuine Hill feedback loop): the
failing Hill directions are **scalings**, and in the joint `(x,θ)` coordinates they have
**integer weights** — e.g. states `+1`, `k_prA +1`, `k_stim +1`, `k_inh −nhill`. The
existing **`.sym_peel_scalings`** recovers exactly this kind of vector (it searches for an
integer-weight vector in the modular nullspace, and the clean integer-weight vector is in
it). **No rational interpolation, no gauge fight, no `df`-stacked general reconstruction is
needed for the SMAD cases** — keep states as dual leaves, add the `df_rest` rows, and let
the scaling peel run over the enlarged coordinate set. The eliminated path fails precisely
because `x*_rest(θ)` is algebraic-irrational (e.g. a square root for `nhill=2`), so folding
the state out turns the clean integer state-weight into a non-rational parameter expression.

This also answers "can we turn the implicit result into an explicit (pure-parameter) one":
for these directions **no rational pure-parameter form exists** — the joint form is the
only closed form.

### Output behaviour (with closedForm = TRUE)

- The direction is reconstructed symbolically in `(x, θ)`. Entries are low-degree
  rationals.
- For the scaling-type symmetries (all SMAD cases), the **parameter components are
  state-free and compact** (`-nhill·k_inh`, `k_pr_FB3`). Splitting param vs state
  components is trivial (we know which coordinate is which), so the **old pure-parameter
  output falls out for free** — no projection, no `x*(θ)` substitution.
- The **state components** (`ξ_FB3 = FB3`) are extra, genuinely useful info; report them
  optionally (e.g. a `stateVector` field) or drop them.
- Only if a *general* (non-scaling) symmetry has a parameter component that genuinely
  references a state would the pure-parameter form need `x*(θ)` substitution (high
  degree) — but in that case no compact pure-parameter form exists anyway, and the
  joint form is the honest closed answer (better than today's support-only).

The transcendental Hill term is handled exactly as today: the recast `E = FB3^nhill`,
`L = log(FB3)` keeps `f` rational; `E`, `FB3`, `L` are algebraically independent with
the differential relations `E' = nhill·E·FB3'/FB3`, `L' = FB3'/FB3`. In the implicit
path **FB3 stays a coordinate** (it is not solved out).

## Implementation sketch (files and what changes)

Investigate/confirm each before coding; the exact kernel wiring is the main unknown.

1. **`inst/code/symmetryDetectionVersion2.py` — `compileObservabilityTapeMulti` /
   `solveSteadyStateModular`.** Add an implicit constraint mode that:
   - keeps the states as free leaves (do not solve `f=0`, do not seed `icSeed` with
     `x*` + IFT duals),
   - emits the steady-state constraint Jacobian rows `df = ∂f/∂(x, θ)` to be stacked
     with the observability rows (or emits `f` so the kernel can form them),
   - keeps the recast (`_apply_power_recast`) so `f`/`df` stay rational; **do not**
     eliminate the recast base.
   The non-equilibrate path already treats states as free ICs — reuse that machinery and
   add the constraint rows, rather than writing a new observability builder.

2. **`src/symmetry_kernel.cpp` — `symObsNull` / `symObsNullMulti` / `symObsNullChain`.**
   Likely needs to accept extra constraint rows and RREF `[Obs ; df]` together (the
   nullspace then spans the joint `(x, θ)` columns). Alternatively stack on the
   R/Python side and pass a combined matrix. Decide which is cleaner; the kernel return
   shape (`R`, `pivots`, `rank`, `dim`) stays the same, just with more columns
   (`dim = n_states + n_params`).

3. **`R/symmetryDetection.R` — reconstruction + reporting.** The reconstruction
   (`.sym_interpolate_direction`, sparse/general entry fits) is already generic over the
   leaf set, so it handles joint `(x, θ)` columns with no change beyond `znames` now
   including states. Add reporting that splits state-free parameter components (old
   format) from state components (new optional field). The whole **log-coordinate gauge
   becomes unnecessary** — plan to revert it (see working-tree state below).

4. **Conserved moieties.** In the implicit path `∂f/∂x` is *allowed* to be singular
   (we never invert it), so the moiety rank-deficiency that forces "solve also for
   totals/params" in the explicit path is a non-issue — the conservation laws live in
   `ker(df)` naturally. Confirm the totals handling composes (it may simplify).

## Open questions — status after the standalone harness (2026-06-30)

1. **Exact determining system.** RESOLVED for the cases of interest:
   `ker([ Obs(f_stim) ; df(f_rest) ])` **evaluated on the resting manifold** is the right
   object. Spurious pure-state directions (`ξ_θ=0`) can appear (and an unstimulated toy
   showed an inflated nullity because a steady-state IC gives a constant trajectory);
   filter the reported set to directions with **nonzero parameter support** (and, for
   conserved moieties, the conservation tangents live in `ker(df)` naturally — see §4).
   The real SMAD model is stimulated, so the observability is non-degenerate.
2. **Generic vs on-manifold sampling.** RESOLVED: **must be on-manifold.** Generic
   off-manifold sampling gives the wrong (high-degree) nullspace; `df·ξ=f` vanishes only
   on `{f=0}`. Keep `solveSteadyStateModular` for sample placement.
3. **Mechanism.** RESOLVED for SMAD: the directions are **integer-weight scalings over the
   joint `(x,θ)` coordinates** → recovered by the existing **`.sym_peel_scalings`**, no
   rational reconstruction. Keep the general `df`-stack + rational path only as a fallback
   for hypothetical non-scaling joint directions.
4. **Relation to the existing non-equilibrate mode.** Still worth checking at the real
   harness: the non-equilibrate path already treats states as free-IC leaves; the new
   pieces are (a) seed those state leaves to the **resting** `x*` value (not a generic IC)
   so samples sit on `{f_rest=0}`, (b) keep each state's **own** dual column instead of
   folding IFT duals into params, (c) add the `df_rest` constraint rows (stack on the R
   side via `.sym_rref_modp`; no C++ change needed — the kernel can't accept extra rows but
   `nz` is fully data-driven so extra state *columns* are fine), (d) run `.sym_peel_scalings`
   over the enlarged coordinate set, (e) report state-weight vs param-weight split.

## Current working-tree state (uncommitted) — decide fate first

`R/symmetryDetection.R` currently contains the log-coordinate-gauge implementation
(superseded by this plan):
- new `.sym_logcoord_gauge` and `.sym_logcoord_backsub`;
- `zvals` threading through `.sym_interpolate_direction`, `.sym_sparse_entry`,
  `.sym_general_rational_entry`, the `nv` wrapper, `relProbe` (now carries `zvals`),
  and `.sym_verify_direction` (3rd arg);
- a `logCoords` flag in `reconstructOne` and the stage-3 wiring after the canon block.
`tests/testthat/test-symmetryDetection.R` has a new passing unit test
("log-coordinate gauge reconstructs a parameter-weighted scaling").

**Recommendation: revert the log-gauge changes** (`git checkout -- R/symmetryDetection.R`
and remove the new test) before starting the implicit work — the implicit path makes the
gauge unnecessary and the dormant code adds confusion. Keep the synthetic-harness *test
pattern* (the hand-built `fakeKcall`) as a model for testing the implicit reconstruction.
The full suite is currently green; confirm again after reverting.

## Verification plan

1. **Small fast model first.** Build a single Hill feedback in a short loop (so the
   explicit path either closes trivially or is fast) and use it to settle the open
   questions (determining system, generic sampling) at the fast harness. Note: earlier
   small models with a feedback Hill loop either closed under the explicit path or were
   slow (>400 s) — pick the model so the *implicit* path is exercised and fast; iterate
   there, not on SMAD.
2. **SMAD watchdog run (acceptance).** The two directions must close with compact
   state-free parameter components (`ξ_{k_inh} = -nhill·k_inh`, etc.), rank stays 58/65.
   Runner (scratchpad is session-local; reproduce):
   ```r
   suppressMessages(devtools::load_all(".", quiet=TRUE))
   ex <- readLines("inst/examples/symmetryDetection.R")
   eval(parse(text = paste(ex[266:343], collapse = "\n")))   # model + observables + events + cond.grid
   smad <- symmetryDetection(reactions, observables, method="observability",
     events=events, conditions=cond.grid,
     forcings=c("bool_ActD","bool_CHX","bool_MG132","TGFb"),
     t0=0, equilibrate=TRUE, reduceCQ=TRUE, closedForm=TRUE, cores=6)
   ```
   Always wrap in `timeout`/`R.utils::withTimeout`; **never run several SMAD jobs in
   parallel** (oversubscribes the machine). The Bash tool's own cap is 2 min — run such
   jobs with `run_in_background` and a generous inner `timeout`.
3. `devtools::test(filter = "symmetryDetection")` stays green.

## Key code references

- `solveSteadyStateModular`, recast + inverted logic: `inst/code/symmetryDetectionVersion2.py`
  ~1155–1500 (recast detection `_detect_power_atoms` ~1494, application `_apply_power_recast`
  ~1513, inversion decision ~1677–1698).
- `compileObservabilityTapeMulti`: same file ~1575+, `znames`/`zSlots` emit ~1814–1826.
- Observability kernels: `src/symmetry_kernel.cpp` — `symObsNull` ~601, `symObsNullMulti`
  ~664, `symObsNullChain` ~807, `symSolveMod` ~962. Return list:
  `(ok, R[rank×nz], pivots[0-indexed], rank, dim=nz)`.
- R reconstruction: `.sym_interpolate_direction` (~1560 pre-edit numbering),
  `.sym_general_rational_entry`, `.sym_sparse_entry`, `.sym_null_residues` (~623),
  `.sym_verify_in_nullspace` (~870, gauge-independent, sound).
