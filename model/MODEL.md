# `cmfate` — what the model is, how it works, and how the data tuned it

A methods and provenance document for the Tier-1 cardiomyocyte cell-cycle fate
network. Companion to [RESULTS.md](RESULTS.md) (what came out) and
[TODO.md](TODO.md) (what is still open). Every number below is either read from
the spec files or produced by code in this directory.

---

## The short answer: which data tuned this?

**Not a combination of Baniol and our knockout data.** The three datasets in play
have three different jobs, and only one of them ever touches a parameter.

| Dataset | What it contributed | Fitted parameters it sets |
|---|---|---|
| **Murganti 2022** (hiPSC-CM, live imaging) | The **only** calibration data. Four fate fractions from Fig 2A, inverted through the fate partition. | **4** |
| **Baniol 2021** (285-cell Smart-seq2) | **Structure** — which nodes exist, which edges are real, which gene is the right anchor, and two measured maturation couplings. Plus held-out validation targets. | **0** |
| **Our E2f7/E2f8 KO** (30,030 cells) | **Nothing.** It is the motivation at the front and a prediction target at the back. | **0** |

The model has never been fitted to our own data, and currently *cannot* be scored
against most of it — only **29 of 63** model node genes are present in the app's
2,181-gene curated panel. That is a limitation, but it is also why the in-silico
knockout in [RESULTS.md §7.5](RESULTS.md) is a genuine prediction rather than a
restatement: nothing about our cells was available to tune against.

Verifiable directly — the model package contains no reference to the knockout
export at all:

```bash
grep -rn "ko_export\|app_data\|tbl_cellcycle" cmcycle/ tests/    # no matches
```

---

## 1. What the model is

`cmfate` is a **normalized-Hill logic ODE** (Netflux formalism; Kraeutler, Soltis
& Saucerman 2010) over **55 nodes** and **77 reactions**. It takes eight inputs
and returns the four cardiomyocyte cell-cycle outcomes.

```
inputs (8)          Clonidine, BetaAR, Nrg1, IGF1, MechLoad,
                    ROSenv, Maturation, InVitro
        ↓
signalling → E2F layer → cycle & cytokinesis machinery
        ↓
gates (3)           SPhase, MitoticEntry, Abscission
        ↓
fates (4)           Quiescent, Division, Binucleation, Polyploidization
```

The whole network is **three tracked text files** — one row per species, one row
per reaction, and a manifest — rather than a binary workbook, because a curated
network is reviewed one edge at a time and every reaction carries its own
`evidence` column.

| File | Contents |
|---|---|
| [`cmfate_species.csv`](cmcycle/data/cmfate_species.csv) | 55 nodes: `Yinit, Ymax, tau, module, genes, data_orient` |
| [`cmfate_reactions.csv`](cmcycle/data/cmfate_reactions.csv) | 77 rules: `Rule, Weight, n, EC50, evidence` |
| [`cmfate_model.toml`](cmcycle/data/cmfate_model.toml) | fate/gate declarations, `[calibration]`, six named contexts |

Six contexts are defined purely as input vectors — `hipsc_cm` (M = 0.12),
`mncm_invitro` (0.50), `mouse_p1_invivo` (0.50), `mouse_p0_invivo` (0.30),
`mouse_p7_invivo` (0.55), `adult` (0.95).

**What Tier 1 deliberately cannot do.** Absolute phase durations, duration
*distributions*, ploidy and cell counting, cumulative-EdU versus instantaneous
Ki-67, or FUCCI trace shape. `tau` here is a relaxation constant, **not a phase
duration** — matching three measured durations with three `tau` values would be a
definition, not a prediction. Those questions are what the specified-but-unbuilt
Tier-2 mechanistic ODE exists for.

---

## 2. How the engine works

A state vector goes in and a derivative comes out:

```
dY/dt = (drive · Ymax − Y) / tau
```

Four stages, two of which are composition rules that are worth understanding
because **both of this model's structural bugs lived in them**.

**Per reaction — AND.** Each reactant is passed through `act()` (or `inhib() = w −
act()` for a negated reactant) and the results are combined by AND. The
`(w, n, EC50)` triple belongs to the **reaction**, not the reactant, so one curve
applies to every reactant in that rule.

**Per node — weighted OR.** Every reaction pointing at a node is combined by
weighted OR. Because OR can only *add* drive, **no route can be vetoed by
another**. An inhibitor must therefore be an AND term *inside* a reaction, never a
separate reaction pointing at the same node.

That single property produced the same class of bug four separate times. The worst
instance: an OR'd `Midbody` route bypassed the obligatory Ect2/RhoA arm entirely,
which made it structurally *impossible* for Ect2 to be rate-limiting — the model
could not express its own central claim. Making RhoA obligatory widened the arm's
span from 1.6× to 42× between hiPSC and mouse P1. A passing test suite would have
hidden all four.

**The activation function.** `act()` is the normalized Hill curve, reparameterized
so `(n, EC50)` set shape and half-max directly. It carries a guard: `EC50^n < 0.5`,
enforced at load. This is not a correction to the published formalism — Netflux's
own defaults (n = 1.4, EC50 = 0.50) sit safely inside the ceiling of 2^(−1/n). A
*high-threshold gate* is what walks into it, and this model has three; at
EC50 = 0.70, n = 1.4 the intermediate β goes negative and activities are pushed
above 1.

**Perturbations enter through `Ymax`.** Knockdown is `Ymax → 0`, overexpression
holds `Y = 1`, graded knockdown scales `Ymax`. Output-only nodes never enter the
integrator — their steady state is algebraic, which matters because the fate nodes
carry `tau` up to 35 h.

### The fate layer

The four fates are a **complementary product over the three gates**, so they sum
to exactly 1 with no normalization step and no engine change:

```
!SPhase                              => Quiescent
SPhase & MitoticEntry & Abscission   => Division
SPhase & MitoticEntry & !Abscission  => Binucleation
SPhase & !MitoticEntry               => Polyploidization
```

Writing `g()` for the shared activation function, at steady state:

```
Q = 1 − g(S)
D = g(S)·g(M)·g(A)
B = g(S)·g(M)·(1 − g(A))          =>  Q + D + B + P = 1  identically
P = g(S)·(1 − g(M))
```

Three structural invariants make the identity hold, and the linter enforces all
three: exactly one reaction per fate node, `Weight = 1` on all four, and identical
`(n, EC50)` across all four. Verified to **1e−6** in every context and at every
dose. The identity is exact only at steady state; during a transient the sum
deviates by the differing `tau`, and that deviation is informative — so
normalization belongs in a reporting helper and never in the model.

---

## 3. How the data tuned it

### 3.1 The complete parameter budget

77 reactions each carry `(Weight, n, EC50)`. Here is where every one of those
numbers comes from:

| Category | Count | Source |
|---|---|---|
| **Fitted by solver** | **4** | Murganti Fig 2A, at `hipsc_cm` only |
| Input weights | 8 | *Are* the stimulus — these define the six contexts, not free parameters |
| `(n, EC50)` shape params | 77 reactions → **6 distinct classes** | Netflux defaults, plus 3 deliberate gates |
| Internal weights | 69 | Hand-set from literature and evidence, on a coarse grid |

Two facts about that table matter more than the counts.

**`(n, EC50)` are class-assigned, not per-reaction fitted.** All 77 reactions draw
from just six parameter classes:

| n | EC50 | Reactions | What it is |
|---|---|---|---|
| 1.4 | 0.35 | 47 | amplifying default |
| 1.4 | 0.50 | 21 | Netflux neutral default |
| 1.4 | 0.40 | 4 | the fate layer (invariant: all four identical) |
| 3.0 | 0.50 | 2 | restriction-point switch |
| 3.0 | 0.45 | 2 | restriction-point switch |
| 2.0 | 0.090 | 1 | the abscission switch (`r080`) |

**Weights sit on a coarse grid.** Of 19 distinct weight values, 16 fall exactly on
a 0.05 grid. That is the signature of hand-setting from evidence, not numerical
optimization — an optimizer would produce arbitrary decimals. The three off-grid
values are each explained in the `evidence` column (`0.434` is `0.62 × 0.70`, the
CycD arm scaled to position the restriction-point fold inside the context range).

### 3.2 The four fitted numbers

The `[calibration]` block is the whole of it:

```toml
[calibration]
context = "hipsc_cm"
targets = { SPhase = 0.0992, MitoticEntry = 0.7698, Abscission = 0.5028 }
knobs   = { SPhase = "r063", MitoticEntry = "r069", Abscission = "r080" }
source  = "Murganti 2022 Fig 2A, n=570 cells from 20 movies, inverted through the fate partition"
```

Three gate weights, plus the *position* of the abscission switch (`r080`'s
`n = 2.0, EC50 = 0.090`), which was chosen by hand rather than solved. That is four
quantities.

**Why the inversion is closed-form.** Competence is deliberately *not* gated on
prevalence, so the three gates are independent by construction. Murganti's four
measured fate fractions (Q 90.35%, D 5.09%, B 3.16%, P 1.40%) therefore invert
into three **separate one-dimensional solves**, not a joint fit:

```
g(SPhase)       = 1 − Q         = 0.0965
g(MitoticEntry) = (D+B)/g(S)    = 0.8549
g(Abscission)   = D/(D+B)       = 0.6170
```

**A subtlety worth documenting, because it looks like a bug.** Those three
`g`-values are *drive* values. The `targets` in the manifest — 0.0992, 0.7698,
0.5028 — are the corresponding **node activities**, i.e. the `Y` at which the fate
layer's shared curve (n = 1.4, EC50 = 0.40) returns the `g` above. Both sets are
correct; they live in different spaces. Confirmed numerically:

| Gate | node activity (manifest target) | → `g(x)` | Fig 2A inversion |
|---|---|---|---|
| SPhase | 0.0992 | 0.0965 | 0.0965 ✓ |
| MitoticEntry | 0.7698 | 0.8549 | 0.8549 ✓ |
| Abscission | 0.5028 | 0.6170 | 0.6170 ✓ |

Anyone comparing RESULTS.md §7.1 (which quotes the `g`-values) against the
manifest (which quotes activities) will think the numbers have drifted. They have
not.

**The solve itself** is one bisection per gate, 32 iterations, pure standard
library — `spec.calibrate()`. It hits all three targets to 1e−9 and takes ~70 s,
cached beside the spec against a SHA-256 fingerprint of the three spec files, so
editing any of them invalidates the cache automatically. It also fails loudly and
usefully:

> `SPhase: even w=1 on r063 only reaches 0.0412 < 0.0992. The upstream drive is too weak — raise an activator, not this knob.`

Result: the model reproduces the fit context **exactly** — Q 0.9035, D 0.0509,
B 0.0316, P 0.0140 against Murganti's 0.9035 / 0.0509 / 0.0316 / 0.0140.

**Say plainly what that means: the fate layer is fitted exactly and predicts
nothing.** All predictive content is downstream of those four numbers.

### 3.3 What Baniol contributed — structure, not parameters

The 285-cell re-analysis constrains the network's **shape**. Four distinct ways,
none of which sets a fitted value:

**Two measured maturation couplings.** `M = mean z(FAO) − mean z(glycolysis)` per
cell, over 12 genes each. Within cycling ventricular cells alone (n = 89, P0 and
P7 pooled so stage cannot be doing the work), Ect2 correlates with M at
**r = −0.563** and E2f6 at **+0.396**. These are the model's two load-bearing
couplings, and making them *measured functions* removes its largest degree of
freedom.

**Node identity — the E2F family must split.** E2f7 is flat across maturation
(t = −0.20) while E2f8 rises (t = +2.28), so they cannot be one node. E2f6 is the
only member expressed in noncycling cells (0.38 in 100% of noncycling P7 vCM,
t = +4.12) and E2f1~E2f6 is only +0.18 n.s., so E2f6 belongs on its own
maturation-driven branch rather than the activator programme.

**Edges tested before being wired.** Spearman on candidate edges *first*:
E2f1~E2f7 **+0.49** and E2f1~E2f8 **+0.46** supported the delayed-feedback arm;
E2f8~Ccna2 **−0.56** confirmed E2F8 is G1/S-restricted. One negative result
mattered as much — Gmnn~Cdt1 correlates **positively** (+0.29) because the
oscillation is post-translational, so **scRNA-seq cannot validate the FUCCI
reporter layer at all**.

**Five gene anchors the data overturned.** Detection rates were checked before
anything was wired: `Wee1 → Pkmyt1`, `Mapk14 → Mapk12` (p38γ, not p38α),
`Cdh1 → Fzr1` (*Cdh1* is E-cadherin, detected in 4→0% of cells — wiring mitotic
exit to it would have hit an absent adhesion gene), `Ccnd1 → Ccnd2`,
`Adra1a → Adra1b`.

This provenance is carried per row. Of 55 species, **40 name their genes** and 34
declare a `data_orient` (33 `direct`, 1 `inverse` for Rb, 6 explicitly `excluded`
where the transcript is not a valid proxy).

### 3.4 Literature-fixed, and adjusted-during-development

**Literature-fixed.** Netflux defaults for shape; Liu 2019 (β-adrenergic/PKA
blocks CM cytokinesis) for the adrenergic arm; Puente 2014 for the oxidative-stress
DDR arm; Gérard & Goldbeter 2009 for cell-cycle topology.

**Adjusted during development — explicitly not independent evidence.** Recorded
so the count stays honest: the E2F repressor-pool weights (an OR pool of three
saturates, so "any repressor suffices" became "the pool is always on"); the
abscission arm's parameter class (attenuating → amplifying, since four hops at
EC50 = 0.50 delivered E2F activity of 0.4 as 0.02); the CycD arm scaled 0.70×; and
the removal of two structures that defeated the model's own claim — an OR'd floor
on Ect2, and the OR'd Midbody bypass described in §2.

### 3.5 The constraint budget

There are roughly **20–25 effective independent constraints** available across
both papers, which puts the defensible ceiling at about **13 fitted parameters**.
The model uses **4**. The rest of its numbers must stay literature-fixed or
evidence-derived, and that is the discipline the `evidence` column exists to
enforce.

---

## 4. What that buys, and what it costs

| | Status |
|---|---|
| Fit context (`hipsc_cm`) reproduced | exact — all four fractions to 1e−5 |
| **Clonidine triad — 3 contexts, 3 outcomes** | **3/3 dominant fates correct**, held out; only `Maturation` and `InVitro` change |
| **Mouse in-vivo cycling level** | **17.70% observed vs 17.17% ceiling = 1.03×, nothing fitted** |
| Entry-fold magnitudes | 1.02 / 2.13 / 1.83 vs observed 2.44 / 2.12 / 1.52; mean error 26% |
| hiPSC drug response | **under-predicted ~2×** — a structural ceiling, not a missing edge |
| Maturation slope of entry | **13.1× modelled vs 3.2× measured** — open, pinned by a test that fails by design |

The mechanism behind the triad is a clean 2×2: **maturation** closes the
abscission arm (`E2Fact → Ect2 → RhoA → Midbody → Abscission`), and **culture**
closes the mitotic-entry brake (`ROS → DDR → Ccng1/Pkmyt1 → MitoticEntry`,
collapsing mitotic entry from 0.531 in vivo to 0.251 in vitro). Which of the two
shuts first is the whole 2×2.

Two residuals are recorded rather than tuned away. The hiPSC miss is provably
structural — at low maturation everything the drug could relieve is already near
maximum (Rb 0.14, p21 0.34, p27 0.02), and gain is capped at the reaction weight,
so no re-wiring *within this formalism* can give a larger response at low
maturation than at high. Two candidate fixes were tested and **rejected**:
`Autophagy ⊣ p21` helps hiPSC (1.02 → 1.49×) but blows mNCM out to 6.70×, taking
mean error from 26% to **121%**; gating it on `!Maturation` gives 55%. Both are
worse than doing nothing.

---

## 5. Where the E2f7/E2f8 knockout data actually sits

**Today: input to nothing, output target for one prediction.**

It motivates the project — cycling rises in the KO (cardiomyocytes, P7: 31.6% vs
25.6% WT) but nothing in that measurement says which *fate* the extra flux
reached, and the app's `sig_ploidy = sig_prolif − sig_cytokinesis` heuristic is
structurally blind to the difference (the two modules fall in lockstep, Δ = +0.005,
p = 0.93).

The model's corresponding prediction is that in-silico E2f7/E2f8 double knockdown
raises Ect2 at **both** ages (+0.448 at P0, +0.258 at P7) but converts that into
division only at P0 (+0.107 vs +0.005, a **22×** difference), with real epistasis
— the double is not the sum of the singles (+0.18 on Ect2 at P0), which matters
because our data *is* a double knockout.

**What blocks scoring it.** Only **29 of 63** model node genes are in the curated
panel; `E2f2`–`E2f6`, `Rb1`, `Ccng1`, `Chek1`, `Wee1`, `Pkmyt1`, `Rhoa`, `Nisch`,
`Adrb1` and `Mapk12` are all missing. Widening the panel is the single change that
would let our own cells validate the network rather than just its E2F layer.

Two things the KO data *has* already settled, both in [RESULTS.md §8](RESULTS.md):
the knockout is **functionally confirmed** (24 of 26 canonical E2F targets up at
P7, median log2FC +0.91 against a housekeeping baseline of +0.19) even though
E2f7/E2f8 mRNA is *higher* in KO than WT; and the E2f8 phase-scoring confound is
real in principle but negligible in practice (≤ 0.5 pp on the KO−WT cycling gap).

---

## 6. Reproducing and checking this

```bash
cd model

# the calibration, printed
python3 -c "from cmcycle import spec; spec.load_calibrated(verbose=True)"

# the fit context reproduced exactly
python3 -c "
from cmcycle import spec, model
net = spec.load_calibrated()
print(net.ss(inputs=model.contexts(net)['hipsc_cm']))"

# the held-out test
python3 -c "
from cmcycle import spec, model
print(model.clonidine_triad(spec.load_calibrated()))"

python3 -m cmcycle.preflight     # six closed-form cross-paper checks, no solver
python3 -m cmcycle.baniol        # the re-analysis (needs CARDIAC_RNASEQ_ROOT)
python3 -m pytest                # includes two tests that pin known failures
```

Two tests are worth knowing by name.
`test_the_mouse_in_vivo_cycling_level_is_reproduced_unfitted` pins the 1.03×
validation. `test_the_maturation_slope_of_entry_is_too_steep` asserts that the
slope error **still exists** — so it will fail, deliberately, on the day someone
fixes it.

### A note on reading the older write-up

[RESULTS.md](RESULTS.md) §1–§4 (the re-analysis) is verified clean against
`figures/results.json`. **§6 onward has drifted** from the regenerated
`figures/model_results.json` — the G1/S switch fix and the screen rescoring
updated the JSON and the figures but not the prose. Most consequential: §7.5
reports the knockout's pro-division effect as "57× larger at P0", where the
regenerated result is **22×** (+0.1073 vs +0.0049) and `fig10` prints the current
numbers on its own face. Treat the JSON and the figures as authoritative.
