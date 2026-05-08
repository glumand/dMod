---
output:
  html_document: default
  pdf_document: default
---
# FOCEI in dMod: Brainstorming

Living document. Wir bauen das inkrementell aus, während wir die Designentscheidungen abklopfen. Noch kein Plan, noch kein Code.

---

## 1. Ziel

FOCEI (First-Order Conditional Estimation with Interaction) für nichtlineare gemischte Effekte in dMod implementieren. Skopus aktuell:

- Nur FOCEI, keine Stufenleiter über FO/FOCE.
- Allgemeine $\Omega$ mit voller Kovarianzstruktur, Cholesky-parametrisiert.
- Bestehende dMod-Likelihood- und Fehlermodellinfrastruktur (`normL2`, `err`) wiederverwenden.
- Outer-Optimierung mit `trust()`.
- FOCEI-Marginal-Likelihood-Approximation symbolisch ableiten, Modell-/Sensitivitätspipeline bleibt wie sie ist.
- PEtab-Anbindung später.

---

## 2. Theorie

### 2.1 Modell

Für $N$ Subjekte $i = 1, \dots, N$ mit Beobachtungen $y_{ij}$ zu Zeiten $t_{ij}$:
$$
y_{ij} = f(t_{ij}, \varphi_i) + \varepsilon_{ij}, \qquad \varphi_i = h(\theta, \eta_i)
$$
mit
- $\theta$: feste Effekte (Populationsparameter),
- $\eta_i \sim \mathcal{N}(0, \Omega)$: zufällige Effekte (individuelle Abweichungen),
- $\varepsilon_{ij} \sim \mathcal{N}\!\bigl(0, \Sigma_{ij}\bigr)$: Residualfehler. $\Sigma_{ij}$ ist ein **beliebiger benutzerdefinierter Ausdruck**, kann von $f_{ij}$, $\theta$ und $\eta_i$ abhängen. Wir machen keine Annahme über die Form. Wird über das `err`-Argument in `normL2` symbolisch kompiliert.

Die individuelle Parametrisierung $\varphi_i = h(\theta, \eta_i)$ ist ebenfalls beliebig (z.B. $\varphi_i = \theta \exp(\eta_i)$). Wird in dMod heute schon durch `Pexpl` abgebildet.

### 2.2 Marginale Likelihood

Pro Subjekt:
$$
L_i(\theta, \Omega) = \int p(y_i \mid \eta_i, \theta)\, p(\eta_i \mid \Omega)\, \mathrm{d}\eta_i.
$$
Für nichtlineares $f$ analytisch nicht lösbar. Verschiedene Approximationen:

| Methode | Idee |
|---|---|
| FO     | Linearisierung um $\eta_i = 0$ |
| FOCE   | Innen $\eta_i^* = \arg\max p(y_i\mid\eta_i)\,p(\eta_i)$ pro Subjekt, dann Linearisierung um $\eta_i^*$ |
| FOCEI  | FOCE plus Interaktionsterm (wenn $\Sigma$ von $\eta$ abhängt, z.B. proportionaler Fehler) |
| Laplace | 2. Ordnung, volle Hessische am Mode |
| SAEM   | Stochastische Approximation des EM-Algorithmus, Monte-Carlo-basiert |

### 2.3 FOCEI Zielfunktion

Pro Subjekt, ausgewertet bei $\eta_i^*$:
$$
-2\log L_i \;\approx\; \underbrace{\sum_j \!\left[\frac{(y_{ij} - f_{ij})^2}{\Sigma_{ij}} + \log \Sigma_{ij}\right]}_{\text{Residualbeitrag}}
\;+\; \underbrace{\eta_i^{*\top} \Omega^{-1} \eta_i^* + \log|\Omega|}_{\text{Prior auf } \eta}
\;+\; \log|H_i|
$$
mit Informationsmatrix
$$
H_i \;=\; G_i^\top \Sigma_i^{-1} G_i \;+\; \Omega^{-1} \;+\; \text{(Interaktion via } \partial \Sigma / \partial \eta\text{)}, \qquad
G_i = \frac{\partial f}{\partial \eta}\bigg|_{\eta_i^*}.
$$

Außen wird über $(\theta, \Omega, \Sigma)$ minimiert, innen pro Subjekt $\eta_i^*$ gesucht. Bilevel.

### 2.4 Prozeduraler Ablauf, und Abgrenzung gegen EM

FOCEI ist **geschachtelt**, nicht alternierend. Pro Outer-Schritt wird das Inner für jedes Subjekt voll konvergiert, *bevor* das Outer überhaupt einen Schritt macht.

```
Initialisiere (theta, Omega, Sigma) und eta_i = 0 fuer alle i.

repeat   // Outer-Schleife (trust)
  // Outer fragt OFV an aktueller (theta, Omega, Sigma) an
  for i = 1..N:                                      // Inner pro Subjekt
    eta_i* <- argmin_eta J(theta, eta, Omega, Sigma; data_i)
              via inneres trust(), warmstart aus letzter Outer-Iter
    H_i    <- Hessian von J an eta_i* (faellt aus trust raus)
  end
  OFV <- sum_i [ J_i(eta_i*)  +  log|H_i| ]          // FOCEI-Aggregation
  dOFV <- analytischer Gradient nach (theta, Omega, Sigma)
  // Outer verarbeitet (OFV, dOFV) und macht einen Schritt
until Outer konvergiert
```

Der entscheidende Punkt: das Inner wird in jeder Outer-Iteration *vollständig* gelöst. Es gibt keinen "burn-in dann alternierend"-Modus und keine partielle Inner-Konvergenz.

**Warum das oft fälschlich als alternierend wahrgenommen wird.** In der Praxis startet die Inner-Suche pro Subjekt mit dem $\eta_i^*$ aus der vorherigen Outer-Iteration (Warmstart). Wenn das Outer nahe am Optimum ist, ist das Inner schon nach 1 bis 2 Newton-Schritten wieder konvergent. Das *Verhalten* sieht aus wie "Outer-Schritt, Inner-Update, Outer-Schritt, Inner-Update", aber formal ist es strikt geschachtelt: das Outer wartet auf vollständige Inner-Konvergenz.

**Abgrenzung zu echten EM-Verfahren:**

| Verfahren | Inner-Anteil | Outer-Anteil | Ablauf |
|---|---|---|---|
| **FOCE / FOCEI** | Mode-Suche $\eta_i^* = \arg\max p(y_i\mid\eta_i)p(\eta_i)$ pro Subjekt | minimiert approx. marginal $-2\log L$ über $(\theta,\Omega,\Sigma)$ | **geschachtelt**, Inner stets voll konvergent |
| Laplace | wie FOCE, mit voller Hessischer | wie FOCE | geschachtelt |
| **SAEM** | Monte-Carlo-Sampling aus Posterior $p(\eta_i\mid y_i,\theta)$ | M-step über $(\theta,\Omega,\Sigma)$ | **echt alternierend**, stochastische Approximation, Burn-in dann Annealing |
| EM (analytisch) | E-step: Posterior-Erwartung | M-step: Maximierung | alternierend, für NLME aber nicht tractable |

**Was es bei FOCEI also nicht gibt:**

- Keinen Burn-in mit nachgeschalteter Annealing-Phase. Das ist ein SAEM-Konzept (Monolix), nötig wegen MC-Rauschen.
- Keine zwei Phasen mit unterschiedlichem Algorithmusverhalten.
- Keine Stochastik im Optimierungspfad (deterministische Bilevel-Optimierung).

**Was FOCEI in der Praxis trotzdem als "Phasen-Verhalten" sehen kann:**

- Anfangs (weit weg vom Optimum) braucht das Inner mehr Newton-Schritte pro Subjekt, später (nahe Optimum mit Warmstart) sehr wenige. Das ist Effizienz-Effekt, keine Algorithmus-Phase.
- Eine *Inner-Toleranz*, die mit der Outer-Gradientennorm skaliert (anfangs locker, später strenger), kann die Effizienz weiter verbessern. Das ist eine Implementationsoption, kein konzeptioneller Phasenwechsel.

In dieser Implementierung also: ein einziger geschachtelter Loop, Warmstart pro Subjekt, optional adaptive Inner-Toleranz. Kein EM, kein Burn-in.

### 2.5 Verwandte Verfahren und Architekturentscheidung

Sobald wir die Frage "warum nicht alternierend?" konsequent zu Ende denken, landen wir bei verschiedenen verwandten Methoden, die die gleiche Laplace-Approximation nutzen, aber unterschiedlich daraus eine Schätzgleichung bauen. Lohnt sich ein Überblick, weil die *Mathematik unterhalb* fast identisch ist und wir sie nur einmal bauen müssen.

#### Verfahren im Vergleich

Sei $J_i(\theta, \eta_i, \Omega, \Sigma)$ die joint negative log density (Residual + Prior auf $\eta$, ohne $\log|H|$), $\eta_i^*$ ihr Mode, $H_i = \partial^2 J_i / \partial \eta^2$ am Mode.

| Verfahren | Outer-Zielfunktion | Optimierungsstrategie | Ω-Update |
|---|---|---|---|
| **FOCEI** | $\sum_i [J_i(\eta_i^*) + \log|H_i|]$ (Laplace-Approx von $-2\log L_{\text{marginal}}$) | nested, Outer mit voller Chain Rule durch $\eta^*$ | als Teil des Outers numerisch |
| **Laplace 2. Ordnung** | wie FOCEI plus Korrekturterme aus echter Hessischer | nested | wie FOCEI |
| **Laplace-EM** | ELBO $\sum_i \bigl[-\tfrac{1}{2} E_{q_i}[J_i] + \tfrac{1}{2}\log|H_i^{-1}|\bigr]$ mit $q_i = \mathcal{N}(\eta_i^*, H_i^{-1})$ | alternierend (E-step: $\eta_i^*$, $H_i$; M-step: $\theta, \Omega, \Sigma$) | **closed form**: $\Omega^{(t+1)} = \tfrac{1}{N}\sum_i (\eta_i^* {\eta_i^*}^\top + H_i^{-1})$ |
| **ITS** | $\sum_i J_i$ allein (kein $\log|H|$) | alternierend, Block-Koordinatenabstieg | Stichprobenkovarianz von $\eta_i^*$ |
| **SAEM** | stochastische ELBO mit MC-Posterior | alternierend, mit Annealing | closed form aus MC-Samples |

Beobachtung: alle Verfahren teilen denselben Inner-Mode $\eta_i^*$ und die gleiche Hessische $H_i$, sie unterscheiden sich nur in (a) der skalaren Outer-Zielfunktion und (b) der Strategie zum Optimieren. Die teure Mathematik (ODE-Integration, Sensitivitäten, Cholesky-Ableitungen) ist gemeinsam.

#### Wer macht das heute mit exakten Ableitungen?

| Software | Verfahren | Sensitivitäten |
|---|---|---|
| NONMEM | FO, FOCE, FOCEI, Laplace 2. Ordnung | überwiegend FD |
| Monolix | SAEM | überwiegend FD |
| nlmixr2 | FOCEI, SAEM | gemischt, AD nur in Teilen |
| TMB / glmmTMB | Laplace direkt (nicht EM) | volle AD, *aber* GLMM/lineare Modelle, kein ODE-RHS |
| Pumas (Julia) | FOCEI, Laplace | volle AD inkl. ODE |
| **dMod-Plan** | FOCEI primär, Laplace-EM optional | **volle exakte Ableitungen via CppODE / funCpp** |

Die echte Lücke im R-Ökosystem: **NLME mit ODE-Modellen und exakten Ableitungen durch die komplette Pipeline (Modell, Sensitivitäten, FOCEI-Aggregation, Cholesky-Diff)**. Genau das, was dMods Infrastruktur ohnehin schon hergibt.

#### Warum FOCEI und nicht Laplace-EM als primäres Ziel?

**Pro Laplace-EM:**

- $\Omega$-Update in geschlossener Form ist mathematisch elegant und reduziert die Outer-Dimension um $K(K+1)/2$ Parameter.
- Konzeptionell sauber alternierend, einfacher Pseudocode.
- Variations-Lower-Bound monoton wachsend (zumindest die ELBO, nicht $\log L$).

**Contra Laplace-EM:**

- ELBO und Laplace-Marginal-Likelihood sind **unterschiedliche Funktionen**. Konvergenzpunkt ist nicht der FOCEI-Punkt. Das wirft Validierungsprobleme gegen NONMEM-Referenzen auf.
- EM-Monotonie greift nur auf der ELBO, nicht auf $\log L$. Der psychologische Vorteil "EM konvergiert immer" trügt bei nichtexakter Posterior.
- M-step für $\theta$ und $\Sigma$ bleibt eine nichtlineare Optimierung mit ODE-Auswertungen, also nicht billiger als der FOCEI-Outer.
- Pharmakometrie-Welt kennt FOCEI, Reviewer und Regulierungsbehörden auch. Laplace-EM ist exotisch und braucht zusätzliche Begründung.

**Pro FOCEI:**

- Standardziel, validierbar gegen NONMEM/nlmixr2/Monolix-Referenzen.
- Lingua Franca der Pharmakometrie.
- Mathematik ist sauberer in einer einzigen Outer-Schleife mit `trust`.

**Contra FOCEI:**

- Cross-Term $\partial \log|H|/\partial \eta \cdot \mathrm{d}\eta^*/\mathrm{d}\theta$ in der Outer-Ableitung ist nicht trivial.
- Kein closed-form $\Omega$-Update, $\Omega$ wird im Outer zusammen mit $\theta, \Sigma$ optimiert.

Auf das Wesentliche reduziert: FOCEI ist die **standardvalidierbare Wahl**, Laplace-EM die **strukturell elegantere**.

#### Architektur: ein Math-Layer, mehrere Outer-Strategien

Die saubere Konsequenz: bauen wir die gemeinsame Mathematik einmal sauber auf, sind FOCEI und Laplace-EM (und auch ITS als Burn-in-Initialisierung) verschiedene Outer-Strategien auf demselben Substrat.

```
nlme_components(joint, etaSpec, OmegaSpec)
   ->  list(
         inner_solve(theta, Omega, Sigma)    # liefert eta_i*, H_i pro Subjekt
         joint_at_modes(...)                 # J_i(eta_i*) und Ableitungen
         Hi_logdet(...)                      # log|H_i| und Ableitungen
         Omega_chol_algebra(...)             # Inv, det, dInv/dchol etc.
       )

focei(components)        # nested, Outer = trust auf J + log|H|
laplace_em(components)   # alternierend, closed-form Omega update + Outer auf Rest
its(components)          # alternierend ohne log|H|, fuer Burn-in
```

Der gemeinsame Math-Layer (`nlme_components`) macht den teuren Teil: ODE-Integration, Sensitivitäten, Cholesky-Algebra, alles mit exakten Ableitungen durch die existierende dMod/CppODE-Pipeline. Die Outer-Strategien obendrauf sind dünne Wrapper, die nur die jeweiligen Schätzgleichungen orchestrieren.

#### Empfehlung

1. **FOCEI als primäres Implementierungsziel**, Validierung gegen Benchmark-Modelle aus NONMEM und nlmixr2.
2. **Math-Layer so faktoriert, dass Laplace-EM später als Outer-Variante reinkommt**. Kosten gering, sobald die FOCEI-Maschinerie steht. ITS als optionales Initialisierungsverfahren ebenfalls aus diesem Substrat.
3. **Die echte Story** dieser Implementierung gegenüber dem Stand der Technik ist *nicht* FOCEI vs. Laplace-EM, sondern **NLME mit exakten Ableitungen durch alle Schichten**. Das gilt für FOCEI und Laplace-EM gleichermaßen, ist aber die wirklich neue Lücke im R-Ökosystem.

---

## 3. Existierende dMod-Bausteine

| Baustein | Status |
|---|---|
| ODE-Lösung mit Sensitivitäten | vorhanden via CppODE / deSolve |
| Bedingungsspezifische Parameter | vorhanden via `parlist`, PEtab-Multi-Condition |
| Komponierbare Pipeline `g*x*p` | vorhanden, Kettenregel intern |
| `normL2` Residualbeitrag inkl. Gradient/Hessian | vorhanden |
| `trust()` Optimizer mit Bound Constraints | vorhanden |
| Parallelisierung über Bedingungen via `mclapply` | vorhanden |
| Symbolische Jacobians, AD via CppODE | vorhanden (`derivMode = "dual"/"symbolic"`) |
| Profile Likelihoods | vorhanden |

Die mentale Umstellung: in dMod ist heute "Bedingung" eine Datengruppe. In FOCEI wird **Subjekt = Bedingung**, plus eine populationsweite Kopplung über $\Omega$, die im aktuellen `normL2`-Pfad nicht vorgesehen ist.

---

## 4. Designidee: joint likelihood plus FOCEI-Wrapper

### 4.1 Joint negative log posterior

Der Nutzer schreibt das volle penalized-likelihood-Objektiv für alle Subjekte zusammen:
$$
J(\theta, \eta_{1..N}, \Omega, \Sigma)
\;=\; \sum_i \!\left[\, \sum_j \frac{(y_{ij} - f_{ij})^2}{\Sigma_{ij}} + \log\Sigma_{ij}
\;+\; \eta_i^\top \Omega^{-1} \eta_i + \log|\Omega| \,\right]
$$

In dMod-Komposition (Variante A, der Nutzer baut die joint selbst):
```r
joint <- normL2(data, g * x * p) + constraintL2(mu = etas, Omega = OmegaSpec)
```

`normL2(...)` liefert den Residualanteil mit beliebigem `err`-Modell. Der Prior-Beitrag $\sum_i \eta_i^\top \Omega^{-1} \eta_i + N\log|\Omega|$ kommt aus einer Erweiterung von `constraintL2`: das heutige `constraintL2(mu, sigma)` ist skalar mit diagonaler $\sigma$. Wir geben ihm einen optionalen `Omega`-Pfad, der bei MVN-Cholesky-Spec auf $\eta^\top\Omega^{-1}\eta + \log|\Omega|$ inkl. Ableitungen nach $\eta$ und nach den Cholesky-Einträgen umschaltet. Skelett (objfn-Vertrag, `dP`-Chain-Rule) bleibt identisch.

### 4.2 FOCEI-Wrapper

Eine Funktion, die das Bilevel-Splitting macht:
```r
focei_obj <- focei(
  joint,                     # objfn auf (theta, eta_1..N, Omega, Sigma)
  etaSpec,                   # welche Parameter sind eta, pro Subjekt
  OmegaSpec,                 # welche Parameter parametrisieren Omega (Cholesky)
  innerControl = list(...)   # rtol, maxit, warmstart cache
)
```
`focei_obj` ist wieder ein `objfn`, aber jetzt nur über $(\theta, \Omega, \Sigma)$. Bei jedem Aufruf:

1. Für jedes Subjekt $i$: `trust()` auf `joint` bei fixiertem $(\theta, \Omega, \Sigma)$, frei in $\eta_i$. Liefert $\eta_i^*$ und die innere Hessische $H_i$ zurück.
2. FOCEI-Beitrag: $\mathrm{OFV}_i = J_i(\theta, \eta_i^*, \Omega, \Sigma) + \log|H_i|$.
3. Summe über $i$.
4. Gradient/Hessian nach $(\theta, \Omega, \Sigma)$ aus der vorhandenen `joint`-Maschinerie plus dem $\log|H_i|$-Term.

### 4.3 Drei wichtige Beobachtungen

**(a) Inner und Outer sauber separierbar via Envelope Theorem.**
Bei $\eta_i^*$ ist $\partial J / \partial \eta = 0$. Damit verschwindet $\mathrm{d}\eta_i^*/\mathrm{d}\theta$ aus dem Gradienten von $J_i(\theta, \eta_i^*, \dots)$ nach $\theta$. Wir brauchen nur die direkten partiellen Ableitungen, die `joint` ohnehin liefert. Für die Outer-Likelihood ohne $\log|H_i|$ ist der Outer-Gradient damit *trivial*: `joint` an $(\theta, \eta_i^*, \Omega, \Sigma)$ auswerten und die Komponenten zu $(\theta, \Omega, \Sigma)$ ablesen.

**(b) `trust()` liefert $H_i$ frei Haus.**
Damit ist $\log|H_i|$ in der OFV ohne Extraarbeit auswertbar. Schwieriger: die *Ableitung* von $\log|H_i|$ nach $(\theta, \Omega, \Sigma)$. Drei Optionen:

- $\log|H_i|$ in der Outer-Ableitung schlicht ignorieren. Outer-Schritte werden mit dem dominanten Teil gewählt. Konvergenz langsamer, aber pragmatisch (analog nlmixr2).
- Numerisch differenzieren in der Outer-Schleife. $\dim(\theta) + \dim(\Omega) + \dim(\Sigma)$ extra Inner-Solves pro Outer-Schritt, mit Warmstart bezahlbar.
- Vollanalytisch: bräuchte $\partial^3 f / \partial\eta^2 \partial\theta$. Nicht praktikabel.

Siehe Abschnitt 5 (funCpp-Codegen) für die saubere Variante: kompilierter symbolischer OFV-Knoten mit Gradient durch Chain Rule, der $\log|H_i|$ konsistent mitdifferenziert.

**(c) $H_i$ in Gauss-Newton-Form genügt.**
$$
H_i \;\approx\; G_i^\top \Sigma_i^{-1} G_i + \Omega^{-1} + \text{(Interaktion)}, \qquad G_i = \partial f / \partial \eta_i.
$$
Erste-Ableitungs-Information, die das `prdframe`-`deriv`-Attribut bereits liefert. Wir brauchen `_s2` / `_sdcv`-Files (zweite Ableitungen aus CppODE) **nicht** für FOCEI selbst, nur falls wir später echtes Laplace machen.

### 4.4 Was der Wrapper konkret wissen muss

- **$\eta$**: pro Subjekt $K$-Vektor. Mapping aus `parlist`-Struktur, da Subjekte = Conditions.
- **$\Omega$**: $K \times K$, symmetric PD, $K(K+1)/2$ freie Parameter. Cholesky $\Omega = LL^\top$, $L$ untere Dreiecks, $L_{kk} = \exp(\omega_k)$, $L_{kl}$ frei für $k > l$.
- **$\Sigma$**: im `err`-Modell, das `normL2` schon kennt. In der Outer-Sicht von $\eta$ und $\Omega$ trennbar via Namens-Partition.
- **$\theta$**: der Rest.

Die Aufteilung ist primär eine Namens-Partition über den Gesamtparametervektor, ähnlich wie heute Bedingungs-Parameter zugewiesen werden. Mögliche Hilfsfunktion: `nlmePartition()`.

---

## 5. Codegen via funCpp und 3. Ableitung

### 5.1 Idee

Die OFV-Aggregation pro Subjekt
$$
\mathrm{OFV}_i = \sum_j \frac{r_{ij}^2}{\Sigma_{ij}} + \log\Sigma_{ij}
\;+\; \eta_i^{*\top} \Omega^{-1} \eta_i^*
\;+\; \log|\Omega|
\;+\; \log|H_i(G_i, \Sigma, \Omega)|
$$
ist ein geschlossener symbolischer Ausdruck in den Inputs

- $f_{ij}$, $G_i = \partial f / \partial \eta$ (numerischer ODE/Sensi-Output zur Auswertungszeit),
- $\Omega$-Cholesky (Outer-Param),
- $\Sigma$-Params (Outer-Param, beliebige Form),
- $\eta_i^*$ (aus Inner),
- $y_{ij}$ (Daten).

Genau dafür ist `funCpp` gemacht. Gleiches Pattern wie heute `Y` (Observable) und `Pexpl` (Parametertrafo): symbolischer Ausdruck wird mit Ableitungen kompiliert und sitzt als Knoten in der dMod-Pipeline. Der FOCEI-Wrapper wird dann nur ein weiterer kompilierter Endknoten.

### 5.2 Welche Ableitungsordnung?

Aufschlüsselung des Outer-Gradienten $\mathrm{d}\mathrm{OFV}/\mathrm{d}\theta$ via Chain Rule durch die OFV-Inputs:

- Direkter Pfad: braucht $\partial f/\partial\theta$ (1. Ordnung) und $\partial G/\partial\theta = \partial^2 f / \partial\eta\partial\theta$ (gemischte 2. Ordnung).
- Impliziter Pfad via $\eta_i^*(\theta)$: $\mathrm{d}\eta_i^*/\mathrm{d}\theta = -H_i^{-1} \partial^2 J/\partial\eta\partial\theta$, ebenfalls 2. Ordnung gemischt. Envelope eliminiert den Beitrag im $J$-Anteil, NICHT im $\log|H|$-Anteil.

**Outer-Gradient: maximal 2. Ordnung gemischt.** Das liefert CppODE heute schon via `_s2`.

Outer-Hessian:
- $\mathrm{d}^2\mathrm{OFV}/\mathrm{d}\theta^2$ enthält $\partial^3 f / \partial\eta\partial\theta^2$, also echte 3. Ordnung.

**Outer-Hessian exakt: braucht 3. Ordnung.** Trust-Region kommt mit Gauss-Newton-Approx der Outer-Hessischen üblicherweise gut zurecht (insb. bei exaktem Gradient).

### 5.3 Empfehlung

Stufe 1: exakter Gradient, Gauss-Newton-approximierter Outer-Hessian. Reicht funCpp mit bestehender 2. Ordnung. Kein Codegen-Aufbohrung.

Stufe 2 (optional, später): falls Konvergenz zickt, funCpp auf 3. Ordnung erweitern (sowohl symbolisch als auch in der Dual-AD-Variante, da Forward-Mode-AD kompositionell auf höhere Ordnungen geht). Eigenständiges Subprojekt im CppODE-Repo.

---

## 6. Verbleibende offene Designfragen

1. **Inner-Toleranz adaptiv?** Je näher der Outer am Optimum, desto enger muss das Inner konvergiert sein, sonst rauscht der Outer-Gradient. Standardrezept (Toleranzfaktor x Outer-Gradientennorm) oder konfigurierbar?

2. **Warmstart-Speicher**: Closure-Variable im FOCEI-Wrapper hält $\eta_i^*$ pro Subjekt (bestätigt). Frage: was bei Outer-Reset (Multistart, neuer Datensatz)? Reset-Hook im Wrapper?

3. **Parallelisierung über Subjekte**: `mclapply` über Inner-Solves ist trivial, aber die Reihenfolge der Side-Effekte (Warmstart-Cache schreiben) muss ggf. ein deterministischer Reduce werden.

4. **API-Konvention für `etaSpec` / `OmegaSpec`**: konkretes Format. Idee: `etaSpec = list(eta_kel = c("subj1", "subj2", ...), eta_V = c(...))`, oder über die existierende `parlist`-Struktur abgeleitet. Klärt sich erst in Phase 1 der Implementation.
