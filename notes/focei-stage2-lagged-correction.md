---
output:
  html_document: default
  pdf_document: default
title: "FOCEI Stage 2: Lagged log|H|-Correction for the Outer Gradient"
---

# FOCEI Stage 2: Lagged $\log|H|$-Correction

Design document. Folge-Session zum FOCEI-Stage-1, der bereits in dMod gemerged ist. Ziel: den in Stage 1 weggelassenen Korrekturterm $\partial \log|H_i|/\partial \theta$ einbauen, ohne pro Outer-Iteration die volle Korrektur teuer neu zu rechnen, mit BDF-artigem Lag-Mechanismus.

---

## 1. Ausgangslage aus Stage 1

Stage 1 nutzt das Envelope-Theorem für den Outer-Gradienten:
$$
\frac{\mathrm{d}\,\mathrm{OFV}}{\mathrm{d}\theta}
\;\approx\; \frac{\partial J}{\partial \theta}\bigg|_{\eta_i^*}
\quad\text{(Envelope-Pfad, exakt)}
$$
und lässt den $\log|H_i|$-Korrekturanteil weg. Auf Theophylline (`testing/focei_theophylline.R`) führt das zu Punktschätzern innerhalb 3-6% der NONMEM-FOCEI-Referenz, mit leicht aufgeblähten $\Omega$-Diagonalen (z.B. $\mathrm{sd}(\eta_{Ka}) = 0.63$ vs. Referenz $0.5$).

Der vollständige Outer-Gradient ist
$$
\frac{\mathrm{d}\,\mathrm{OFV}}{\mathrm{d}\theta}
= \frac{\partial J}{\partial \theta}\bigg|_{\eta^*}
+ \underbrace{\frac{\partial \log|H_i|}{\partial \theta}\bigg|_{\eta^* \text{ fix}}
   + \frac{\partial \log|H_i|}{\partial \eta}\bigg|_{\eta^*} \cdot \frac{\mathrm{d}\eta_i^*}{\mathrm{d}\theta}}_{=:\,c(\theta)\,\text{ Korrekturterm}}
$$
Die Korrektur $c(\theta)$ braucht im Allgemeinen 2.-Ordnung-Sensitivitäten ($\partial^2 f / \partial\eta\,\partial\theta$ und $\partial^2 f / \partial\eta^2$), keine 3. Ordnung. Siehe Diskussion im Brainstorming-Doc.

## 2. Zwei Pfade zur Berechnung von $c(\theta)$

| Pfad | Methode | Kosten pro Outer-Schritt | Implementierungsaufwand |
|---|---|---|---|
| **A: FD über $\theta$** | $c(\theta) \approx [\log\|H_i\|(\theta + h e_k) - \log\|H_i\|(\theta - h e_k)] / (2h)$ pro Komponente $k$. Inner wird mit Warmstart neu gelöst. | $\dim(\theta)$ extra Inner-Cycles (mit Warmstart 1-2 Newton-Steps each) | klein, 1-2 Tage |
| **B: 2.-Ordnung-Sens via `_s2`** | Analytisch via $\mathrm{tr}(H_i^{-1} \partial H_i/\partial \theta)$ mit $\partial^2 f$ aus integriertem Sensitivitätssystem | eine erweiterte ODE-Integration, Faktor $\sim \dim(p)+1$ teurer | groß, 1-2 Wochen (CppODE _s2 Integration, Xs/Y/Pexpl Chain Rule) |

Pfad A reicht vermutlich. Stage 2 startet mit A.

## 3. Lagging-Strategie (BDF-Analogie)

### 3.1 Beobachtung

Die Korrektur $c(\theta)$ ändert sich von Outer-Iteration zu Outer-Iteration glatt. Bei BDF wird die Jacobi-Matrix der ODE aus dem gleichen Grund eingefroren und nur sporadisch neu berechnet. Die Newton-Iteration konvergiert trotzdem gegen die exakte implizite Lösung, weil das Newton-Residuum am Ende sauber zu null gefahren wird.

Übertragen auf FOCEI: $c(\theta)$ wird gecacht und nur bei bestimmten Triggern neu gerechnet. Trust konvergiert trotzdem gegen den exakten FOCEI-Fixpunkt, weil

- der OFV-Wert pro Iteration **immer** frisch ist (Inner-Solve + Eigendekomposition von $H_i$),
- der Korrekturterm vor der Konvergenz-Deklaration **immer** einmal frisch geholt wird.

### 3.2 Korrektheitsargument

**OFV-Pfad ist nie gelaggt.** $\log|H_i|$ als Wert in der OFV wird pro Outer-Iteration frisch aus der Eigendekomposition gerechnet, das passiert ohnehin. Trust-Region akzeptiert/ablehnt Schritte basierend auf dem realen OFV-Vergleich (actual reduction). Diese Schicht ist also nicht von Lagging betroffen.

**Gradienten-Pfad ist gelaggt.** Trust nutzt den Gradienten nur für die Wahl der nächsten Suchrichtung (predicted reduction). Wenn der Gradient veraltet ist, wird die Predicted-vs-Actual-Ratio schlecht, Trust schrumpft den Radius, und der Refresh-Trigger feuert. Trust-Region liefert die Feedback-Schleife frei Haus.

**Konvergenz-Test mit Refresh.** Vor der finalen Konvergenz-Deklaration wird $c(\theta)$ einmal explizit gerechnet. Damit ist der Gradient am deklarierten Fixpunkt exakt, und der Punkt ist identisch zur nicht-gelaggten FOCEI-Schätzung.

### 3.3 Refresh-Trigger

```
refresh_correction <- function(state) {
  triggers <- c(
    state$step_rejected,                                    # rho < 0
    state$rho < 0.25,                                       # poor predicted reduction
    norm(state$theta - state$cache_anchor) > tau * norm(state$cache_anchor),
    state$iter %% M == 0,                                   # periodic safety net
    state$step_norm < tol_pre_convergence                   # near convergence
  )
  any(triggers)
}
```

Defaults:
- `tau = 0.1` (relative Anker-Distanz)
- `M = 5` (periodischer Refresh)
- `tol_pre_convergence` = $10\times$ der Trust-Konvergenz-Toleranz (refresh feuert, bevor Trust "konvergiert" sagt)

Nach `refresh_correction == TRUE`: $c(\theta)$ neu rechnen via Pfad A, Cache-Anker auf aktuelles $\theta$ setzen, Cache aktualisieren.

### 3.4 Erwartete Einsparung

Auf Theophylline-Größe ohne Lag: $\dim(\theta) = 6$ FD-Perturbationen pro Outer-Schritt $\times$ 16 Outer = 192 zusätzliche Inner-Cycles (warm).

Mit Lag und obigen Defaults: typisch 3-5 Refreshes über 16 Outer = 18-30 zusätzliche Inner-Cycles. Ersparnis Faktor $\sim 6$-$10$.

Auf größeren Problemen ($\dim(\theta) = 20+$, $N = 100+$ Subjekte) wird der Faktor relativ noch größer, weil die FD-Schicht polynomiell mit $\dim(\theta)$ skaliert während Refreshes seltener werden.

## 4. Implementierungsskizze

### 4.1 Erweiterung des `focei()`-Closures

```r
focei <- function(joint, omegaSpec,
                  innerControl = list(),
                  correction = c("none", "fd_lagged", "fd_eager"),
                  correctionControl = list(h = 1e-4, tau = 0.1, M = 5L),
                  cores = 1L) {
  ...
  # Cache state
  cache_correction <- NULL          # numeric vector of length dim(theta)
  cache_anchor     <- NULL          # theta at which correction was computed
  prev_step_rho    <- NA_real_      # last accepted step's rho
  iter_since_refresh <- 0L
  ...
}
```

### 4.2 Berechnung der Korrektur (Pfad A)

```r
compute_correction <- function(outer_pars, fixed) {
  # 1. Anchor: compute log|H_i| at the unperturbed point (already in cache)
  base_logdetH <- sum(diag_logdet)   # from current evaluation
  
  # 2. FD: for each theta_k, perturb +/- h, run inner solves with warmstart,
  #    accumulate log|H_i| over subjects, take central FD
  correction <- numeric(length(outer_pars))
  for (k in seq_along(outer_pars)) {
    theta_plus  <- outer_pars; theta_plus[k]  <- theta_plus[k]  + h_k
    theta_minus <- outer_pars; theta_minus[k] <- theta_minus[k] - h_k
    
    # Reuse warmstart cache for inner; this is THE key efficiency point
    logdet_plus  <- inner_solve_all(theta_plus,  fixed)$sum_logdetH
    logdet_minus <- inner_solve_all(theta_minus, fixed)$sum_logdetH
    
    correction[k] <- (logdet_plus - logdet_minus) / (2 * h_k)
  }
  setNames(correction, names(outer_pars))
}
```

Schrittweite $h_k$ adaptiv: $h_k = \sqrt{\epsilon_{\text{mach}}} \cdot \max(|\theta_k|, 1)$.

### 4.3 Refresh-Logik im myfn-Body

```r
myfn <- function(..., fixed = NULL, deriv = TRUE, env = NULL) {
  ... (existing inner solves and OFV aggregation, unchanged) ...
  
  # Decide whether to refresh
  do_refresh <- correction == "fd_eager" ||
                is.null(cache_correction) ||
                norm_diff(outer_input, cache_anchor) > tau ||
                iter_since_refresh >= M ||
                last_step_was_rejected
  
  if (deriv && correction != "none" && do_refresh) {
    cache_correction   <<- compute_correction(outer_input, fixed)
    cache_anchor       <<- outer_input
    iter_since_refresh <<- 0L
  } else {
    iter_since_refresh <<- iter_since_refresh + 1L
  }
  
  # Apply correction to gradient
  grad_outer[outer_active] <- grad_outer[outer_active] +
                                cache_correction[outer_active]
  ...
}
```

### 4.4 Trust-Region-Hook für `step_rejected`

Trust gibt aktuell keinen direkten Zugriff auf $\rho$ pro Schritt aus dem objfn-Callback. Zwei Optionen:

- **Option (a)**: focei-Wrapper detektiert Rejection heuristisch über sukzessive Outer-Aufrufe mit gleichem $\theta$ (Trust ruft den objfn nochmal mit dem reduzierten Schritt). Etwas fragil.
- **Option (b)**: kleine Erweiterung von `R/trust.R` um einen optionalen Hook `on_step(rho, accepted)` der bei jeder Schritt-Entscheidung gerufen wird. Sauberer.

Empfohlen: Option (b), zusätzliche Zeile in `R/trust.R` ist trivial und nützlich auch für andere Profile-Likelihood-Erweiterungen.

## 5. Checkliste für die Implementierung

- [ ] `correction` Argument zu `focei()` hinzufügen, Default `"fd_lagged"`.
- [ ] Cache-Variablen im Closure: `cache_correction`, `cache_anchor`, `iter_since_refresh`, `last_rho`.
- [ ] `compute_correction()` als interne Hilfsfunktion, ruft `inner_solve_one` mit perturbiertem $\theta$ und Warmstart aus aktuellem Cache.
- [ ] Refresh-Logik im myfn-Body, Hooks am Anfang und Ende.
- [ ] Optional: `R/trust.R` um `on_step`-Callback erweitern für sauberen Rejection-Trigger.
- [ ] Tests:
   - [ ] Korrektur am konvergierten Punkt von Stage 1: $\|c(\hat\theta_1)\|$ messen, sollte klein aber nicht null sein.
   - [ ] FOCEI mit Korrektur auf Theophylline laufen lassen, Schätzer sollten sich Richtung NONMEM-Referenz bewegen, $\Omega$-Diagonalen sollten sinken.
   - [ ] Lag vs. Eager: identische Endschätzer (innerhalb numerischer Toleranz), aber Lag braucht weniger Inner-Cycles.
   - [ ] Refresh-Counter: prüfen dass auf Theophylline 3-5 Refreshes total auftreten.

## 6. Was bewusst nicht in Stage 2 gehört

- **2.-Ordnung-Sensitivitäten via `_s2`**: das ist Pfad B, separates Folgeprojekt. Erfordert CppODE-Aufbohrung. Erst sinnvoll wenn $\dim(\theta)$ groß genug ist dass FD-Cost dominant wird.
- **Outer-Hessian via 3.-Ordnung-Sens**: erst in Stage 3 relevant, wenn überhaupt. Trust kommt mit Gauss-Newton-Schur-Komplement gut zurecht.
- **Standardfehler-Berechnung**: separater PR.
- **PEtab-Anbindung**: separater PR.
- **SAEM/Laplace-EM**: andere Architektur, nicht in dieser Linie.

## 7. Erwartetes Ergebnis nach Stage 2

Auf Theophylline:
- Punktschätzer für $(K_a, V, Cl)$ innerhalb $\sim 1\%$ der NONMEM-FOCEI-Referenz.
- $\mathrm{sd}(\eta)$-Werte deutlich näher (z.B. $\mathrm{sd}(\eta_{Ka}) \approx 0.5$ statt $0.63$).
- Outer-Konvergenz weiterhin in $\sim 20$ Iterationen.
- Total Inner-Cycles ca. doppelt so viel wie Stage 1 (Refreshes hinzu), nicht $\dim(\theta)$-fach.

Damit erreicht dMod auf der Genauigkeitsachse Parität mit nlmixr2/NONMEM, gewinnt aber durch die analytische Envelope-Komponente plus FD-on-correction-only an Effizienz: NONMEM macht FD über die *gesamte* OFV (envelope wird verschenkt), wir machen FD nur über den kleinen Korrekturterm.

## 8. Referenzen für die Folge-Session

- Brainstorming-Doc: `notes/focei-brainstorm.md` (Theorie, Architekturentscheidungen, Verfahrensvergleich)
- Stage-1-Implementierung: `R/nlme.R`, `R/objClass.R`, Tests in `tests/testthat/test-focei.R`
- Benchmark: `testing/focei_theophylline.R`
- BDF-Lag-Referenz für die Analogie: Hairer & Wanner, "Solving Ordinary Differential Equations II", Kapitel IV.8 (Implementation of implicit Runge-Kutta methods).

---

## 9. Path-B Addendum (analytical correction landed)

Sobald `deriv2` für die gesamte `g*x*p`-Pipeline verfügbar wurde
(`Pexpl/Pequil` als 3D-Tensor, `Xs.CppODE/Y` als 4D-Tensor, plus
Compositions-Chain-Rule), kann der Korrekturterm rein analytisch ausgewertet
werden -- der FD-Pfad in §4.2 entfällt damit. Die in dMod gemergte Stage-2
Implementierung (`R/nlme.R`, signature `focei(joint, omegaSpec, model,
data, correction, correctionControl, ...)`) verwendet konsequent die
**Gauss-Newton-Form**

$$
H_i^{GN} = G_i^\top \Sigma^{-1} G_i + \Omega^{-1}, \qquad G_i = \partial g/\partial \eta\big|_{\eta_i^*},
$$

und damit nur Ableitungen 2. Ordnung. Die explizite Komponente

$$
\partial H_i^{GN}/\partial \theta_k = (\partial^2 g/\partial\eta\partial\theta_k)^\top \Sigma^{-1} G_i
+ G_i^\top \Sigma^{-1} (\partial^2 g/\partial\eta\partial\theta_k)
$$

wird einmal pro Refresh aus einem `model(times, full_pars, deriv = TRUE,
deriv2 = TRUE)`-Call zusammengesetzt; für die Cholesky-Parameter, die in den
Predictions nicht auftauchen, kommt der Beitrag $\partial \Omega^{-1}/\partial
\mathrm{chol\_par}$ über zentrale FD auf `omegaSpec$build_L` (vernachlässigbar
billig).

Der implizite Pfad $\partial \log|H_i|/\partial \eta_l \cdot \mathrm{d}\eta_l^*/\mathrm{d}\theta_k$
nutzt den Cross-Block der GN-Joint-Hessian (`joint_at_modes$hessian[eta, theta]`)
für $\mathrm{d}\eta_l^*/\mathrm{d}\theta_k$ und dieselbe `deriv2`-Slice für
$\partial \log|H_i|/\partial \eta_l = \mathrm{tr}(H_i^{-1} \partial H_i/\partial \eta_l)$.

Wichtige Konvention: dMod's `normL2`/`constraintL2` führen Werte ohne den
$1/2$-Vorfaktor (`(res/sigma)^2`, nicht $\tfrac{1}{2}(res/sigma)^2$), daher
ist die Hessian aus der Joint **2x** der konventionellen Gaussian-NLL-Hessian.
`compute_correction()` mirrors this with an explicit `conv2 = 2` factor on
both $\partial H_i/\partial \theta$ and $\partial H_i/\partial \eta$, so dass
$H_i^{-1}$ aus dem Inner-Solver und $\partial H_i$ aus der Korrektur in
derselben Konvention leben.

Lagging und Refresh-Trigger sind unverändert wie in §3.3 spezifiziert (cold
start, periodisch via `M`, Anker-Distanz via `tau`, Reject via `on_step`-Hook
auf `R/trust.R`, plus implizite Konvergenz-Korrektheit weil OFV stets frisch
ist). Der Hook ist ein optionales `on_step = NULL`-Argument auf `trust()` und
hat keinen Einfluss auf bestehende `mstrust()`/`profile()`-Aufrufer, die ihn
nicht setzen.

Stage-2.5-Restriktion: das Fehlermodell darf $\theta$ enthalten, aber **kein
$\eta$** (kein Interaktionsterm im klassischen FOCEI-Sinn). Diese Erweiterung
braucht `errmodel`-2nd-Order und eine sorgfältigere Kontraktion in
`nll_ALOQ`; bewusst ein Folge-PR.

Verifikation auf einem K=1, N=4-Toy-Modell (`tests/testthat/test-focei.R`):
analytische Korrektur stimmt mit zentraler FD über die volle OFV innerhalb
~10 % überein, Lagged- und Eager-Modus liefern beim Cold-Start identische
Gradienten, und der Lagging-Cache hält erwartungsgemäß über kleine
Parameterschritte und refresht bei großen.
