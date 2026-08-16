# `tier2/` — the mechanistic cardiomyocyte cell-cycle fate ODE

Tier 2 of the two-tier model described in [../MODEL.md](../MODEL.md) §6 and specified in
[../TODO.md](../TODO.md) item 3. Julia, built on the published generic cell-cycle model
in [`Cell_Cycle_Model`](https://github.com/Jwomack7512-bio/Cell_Cycle_Model) (MIT).

| | Tier 1 (`../cmcycle`) | Tier 2 (here) |
|---|---|---|
| Formalism | normalized-Hill logic (Netflux) | mass-action + Michaelis–Menten ODE |
| Nodes | 55 nodes, 77 reactions | 63 species, 218 parameters |
| Time | `tau` is a relaxation constant | real time |
| Fates | product of three steady-state gate activities | emergent from one cell's trajectory |
| Species | 55 nodes | 63 inherited + 9 Tier-2 = 72 |

Tier 1 cannot express absolute phase durations, duration *distributions*, ploidy and cell
counting, cumulative-EdU versus instantaneous Ki-67, or FUCCI trace shape. Those five are
what Tier 2 exists for.

## Why this base and not Gérard & Goldbeter 2009

`../TODO.md` item 3 specifies extending gg2009. The published Julia model is a better
starting point: it supplies everything item 3 asks to inherit (Rb–E2F, four cyclin–CDK
waves, APC/C–Cdh1/Cdc20/Emi1, the Wee1/Cdc25 gate, a checkpoint arm) **plus** the single
highest-value thing item 3 asks to *add* — `CDT1`, `Geminin`, `Geminin_CDT1` as FUCCI
observables, which gg2009 lacks — along with `LMNA/LMNAp` for envelope breakdown, `PTTG1`
for anaphase, and `p21`.

Added in Phase 2: the E2F sub-family split (`E2F6`, `E2F7`, `E2F8`) and the cytokinesis
arm (`Ect2`, `RhoA`, `Centralspindlin`, `AurKB`, `Anillin`, `Midbody`), and the maturation
axis `M`, and the oxidative-stress/DDR arm (`Ccng1`). `p27` is deliberately deferred —
see step 5. The test suite asserts the absent ones are absent, so a half-landed module
cannot go unnoticed.

## Status: Phase 2 complete

**441 tests pass** (52 Phase 0, 139 Phase 1, 75 steps 1–2, 72 step 3, 60 step 4, 43 step 5). The source repo has no test suite, so these
are the first the inherited model has had.

`src/inherited/` is a **byte-identical** copy of the published `model_files/`
(`state.jl`, `parameters.jl`, `diff_eqns.jl`). Do not edit it — extensions belong in
their own files so the diff against the published model stays legible to a reviewer.

Phase 1 added **no biology and no parameters**; it only measures the inherited model.
Phase 2 adds modules that are **inert by default** — the extended model reduces
bit-exactly to the published one until an enable parameter is set.

```bash
cd model/tier2
julia +1.11 --project=. -e 'using Pkg; Pkg.instantiate()'   # first time, ~4 min
julia +1.11 --project=. -e 'using Pkg; Pkg.test()'
```

### Pinned environment

`Manifest.toml` is copied from the published model's own resolved environment, so the
whole solver stack matches what its figures were produced with. This is not fussiness:
**`OrdinaryDiffEqCore` 2.x breaks `Rosenbrock23`'s LU factorisation when the Jacobian is
a `ComponentMatrix`**, and `AutoTsit5(Rosenbrock23())` — the integrator every published
script uses — switches into `Rosenbrock23` on this stiff system. Resolving freely picks
up 2.x and the model does not run at all.

Instantiate from the Manifest. Do not `Pkg.update()` without re-running the golden test.

## What Phase 0 turned up: four different `α` values

`α` is a global time-scale multiplying the entire RHS, so the cycle period scales exactly
as `1/α`. The source repo carries **four** values, and the one in `parameters.jl` is not
the one the published figures were made with:

| `α` | period | provenance |
|---|---|---|
| **1.447** | **28.11 h** | hardcoded in `publication_fig_12_*_test.jl` and `*_varyingABE.jl`; reproduces the committed `experimental_vs_simulation_summary.csv` (DMSO 28.1 h) ✓ |
| 1.60 | 25.43 h | `publication_fig_12_*_varyingVolo.jl` |
| 1.725 | 23.58 h | implied by the default target in `Helpers/estimate_alpha_for_target_dt.jl` ✓ |
| 2.3 | 17.69 h | the default in `parameters.jl` — **matches no published number** |

Note that `publication_fig_12_experimental_vs_simulation.jl` (no `_test` suffix) calls
plain `params()`, i.e. `α = 2.3`, yet the committed CSV it writes reports 28.1 h. The CSV
was produced by the `_test` variant. That is a provenance hazard worth resolving before
submission — the timing numbers in the manuscript depend on which script was run.

Two things were ruled out along the way:

- **Not a measurement artefact.** The threshold-crossing method
  (`Helpers/estimate_alpha_for_target_dt.jl`) and the peak-detection method (the figure
  scripts) are independent implementations and agree to ~2 ms on the baseline
  oscillation. Both are ported here (`doubling_time`, `peak_period`) and the test suite
  asserts they agree.
- **Not the `9ab1b1d` Km rewiring.** Extracting the pre-commit `model_files/` and running
  them gives an identical 17.685 h period, so that commit's parameter-set change
  (`Km_CDC25Cp_CDH1` removed, four `Km_*` added) is numerically inert.

`CmTier2.PUBLISHED_ALPHA` is therefore **1.447**, the value the figures actually used,
with `PARAMFILE_ALPHA = 2.3` kept as a named constant. All four are pinned by test, so
whenever they get reconciled the test fails and forces the decision to be recorded.

## Inherited equation defects, pinned rather than silently fixed

Found by reading the equations rather than the variable names. All three are pinned by
`test/runtests.jl` so Phase 2 cannot fix them by accident:

1. **`ks_CDT1_E2F` and `ks_Geminin_E2F` are constants**, not multiplied by `E2F`
   (`diff_eqns.jl:425,432`), despite the names and the comments ("Synthesis of CDT1 via
   E2F"). CDT1 and geminin currently oscillate purely through degradation. Since FUCCI
   trace shape and the E2F split are both Tier-2 targets, Phase 2 must connect them —
   and that changes FUCCI dynamics, which is exactly why this is a follow-up model
   rather than a patch to the manuscript in flight.
2. **`d.p53p` is missing its `*p53` factor** (`diff_eqns.jl:382`), which `d.p53`'s
   matching loss term carries, so p53 ↔ p53p is not mass-balanced. Harmless in the
   published figures because the DDR arm is dormant (`kf_ATMp = 0`, asserted by test),
   but Phase 2 turns that arm on for `ROSenv`/`Ccng1`.
3. **`E2F` has exactly six downstream consumers** — `CCNE`, `CDC25A`, `CCNA`, `PLK1`,
   `PTTG1`, `EMI1`. The Phase 2 split must rewire each and nothing else; a test asserts
   the count so a newly added consumer cannot be missed.

Also not carried over: `model_files/diff_eqn_plk1p_test.jl`, a 442-line shadow copy of
`diff_eqns.jl` with different `Km` wiring. `include_all_scripts()` in the source repo
glob-includes every `.jl` in a folder, so a stale second `modelDiffEq!` can silently win.

## Phase 1: what the inherited model measures

Landmarks are root-found with `ContinuousCallback` rather than scanned off the saved
grid, because durations are the deliverable, not a diagnostic. Verified not to perturb
the trajectory: the difference against a callback-free solve converges to zero with
solver tolerance (2×10⁻³ at default → 5.9×10⁻⁵ at 10⁻⁶ → 3.3×10⁻⁹ at 10⁻¹⁰), which a
real perturbation would not do.

Cycles open at the **restriction point** (`ppRB`), not at CycE activity. The obvious
CycE-keyed rule is leaky and measurably so: under `con_ABE = 10` µM, `ppRB` peaks at
0.013 — never committing — while CycE-CDK2 still reaches 0.52 and crosses any sensible
marker, because `d.CCNE` carries a `ks_CCNE_pRBE2F` term plus basal synthesis. A
CycE-keyed classifier calls that arrested cell "cycling", which is the one error a fate
classifier must not make.

At the published α = 1.447, over a 400 h settled window:

| quantity | model | published |
|---|---|---|
| cycle period | 28.11 h | 28.1 h (committed CSV) ✓ |
| G1 / S / G2 / M | 17.33 / 8.64 / 1.93 / 0.22 h | — |
| S/G2/M, cyclin-defined | 18.95 h | — |
| S/G2/M, FUCCI-defined | 12.80 h | 16.38 / 17.29 / 24.50 h (Murganti Fig 2E); Baniol 15.1 ± 4.0 ✓ |
| mitosis (NEB → exit) | 1.74 h | ~4.2 h for CM (Tier 1 preflight) |
| cells after 14 cycles | 2¹⁴ = 16384 | — |

Two S/G2/M durations are carried deliberately. Murganti Fig 2E times a FUCCI trace —
mAG (geminin) appearing, to division — so scoring a cyclin-defined duration against it
would be a category error, and not a small one: the two differ by ~6 h here.

> **Note.** Phase 1 originally read the FUCCI duration as 2.06 h and recorded that as a
> known failure. Phase 2 showed the failure was in the ruler, not the model — the cutoff
> was a plotting default. The table above uses the calibrated cutoff. See
> [the retraction](#retracted-the-fucci-defect-was-a-plotting-default).

### Fate classification, and what Phase 1 deliberately will not do

Three of Tier 1's four fates are decidable from the inherited model:

| fate | criterion |
|---|---|
| `:Quiescent` | restriction point never passed |
| `:Polyploidization` | S entry, no nuclear-envelope breakdown (re-replication) |
| `:MitoticCompletion` | S entry → NEB → anaphase → mitotic exit |

`:MitoticCompletion` is Tier 1's Division **+** Binucleation. Splitting it needs the
abscission decision, which needs the Ect2/RhoA/Centralspindlin/AurKB/Anillin/Midbody arm
— absent from the inherited model, added in Phase 2. It is deliberately not guessed at:
nothing in the current equations represents the furrow, so there is no honest way to
call it. `bookkeep(...; assume_abscission = ...)` makes the placeholder explicit at
every call site rather than hiding it in a default.

Note also that for a single deterministic trajectory the fate fractions are 0 or 1 —
whichever the limit cycle settles into. The graded fractions Murganti Fig 2A reports are
a population quantity and need the Phase 3 ensemble, since `../TODO.md` item 4 notes
fate fractions are tail statistics a mean-field treatment cannot reach.

## Layout

```
Project.toml        name/uuid/[compat] — the source repo's env declared none
Manifest.toml       pinned to the published model's resolved stack (see above)
src/
  CmTier2.jl        module: solve_baseline, solve_drug, doubling_time, peak_period
  inherited/        byte-identical copy of the published model_files/ — do not edit
  observables.jl    FUCCI states and fractions, total pools, phase durations
  events.jl         EventThresholds, EventLog, root-found cell-cycle landmarks
  fates.jl          Cycle, classification, nuclei/ploidy/cell bookkeeping
  tier2_model.jl    extended state/params/RHS: E2F split, cytokinesis arm, maturation
  contexts.jl       the six named contexts, read from Tier 1's manifest
test/runtests.jl    441 tests: the Phase 0, 1 and 2 gates
scripts/            (Phase 3+)
```

## Phase 2: steps 1–2

### How the extension preserves the published model

`tier2DiffEq!` **calls** the inherited `modelDiffEq!` and never copies it. `diff_eqns.jl`
destructures `u` and `p` positionally, so an extended `ComponentVector` with its new
components appended *at the end* is destructured to exactly the same 63 states and 218
parameters. Tier 2 is then additive corrections, each carrying a factor that is zero at
default parameters:

```julia
d.CCNE += ks_CCNE_E2F * E2F * (rep - 1)        # zero when rep == 1
d.CDT1 += ks_CDT1_E2F * w   * (E2F/E2F_ref - 1) # zero when w == 0
```

So the reduction is **structural, not a tolerance**. Verified bit-exact —
`max |d_tier2 − d_published| = 0.000e+00` over 100 random states, checked on random
states rather than a trajectory because a trajectory only visits the limit cycle and
would miss a correction that is nonzero elsewhere. Copying the 442-line RHS and editing
it would have recreated exactly the `diff_eqn_plk1p_test.jl` hazard flagged in Phase 0.

Six enable parameters gate everything, all zero by default:
`w_CDT1_E2F`, `w_Geminin_E2F`, `kd_CDT1_CDK2`, `ks_E2F7_E2F`, `ks_E2F8_E2F`, `ks_E2F6`.

### RETRACTED: the FUCCI "defect" was a plotting default

Phase 1 reported that the inherited FUCCI layer had no G1/S state. **That is withdrawn.**
It was an artefact of the 0.05 cutoff, which is the source repo's `plot_fucci_backgrounds`
*plotting* default — not a calibrated value.

The cutoff is a property of the reporter and the microscope, not of the model. It can be
calibrated with no external data at all, because the model measures phase durations two
independent ways: from cyclin peaks (`phase_times`, which never looks at Cdt1 or geminin)
and from the FUCCI channels (`fucci_fractions`, which never looks at a cyclin). Requiring
them to agree gives a sharp optimum:

| cutoff | Cdt1⁺ | cyclin G1 | geminin⁺ | cyclin S+G2+M | double-neg | error |
|---|---|---|---|---|---|---|
| 0.050 | 0.517 | 0.613 | 0.127 | 0.387 | 0.356 | 0.712 |
| 0.025 | 0.581 | 0.613 | 0.314 | 0.387 | 0.105 | 0.211 |
| **0.020** | **0.609** | **0.613** | **0.391** | **0.387** | **0.000** | **0.0006** |
| 0.015 | 0.813 | 0.617 | 0.187 | 0.387 | 0.000 | 0.400 |

Two measurements sharing no equations agree to under one percentage point, and the
double-negative population — which a real asynchronous FUCCI culture does not have —
vanishes exactly there. That is a genuine internal validation of the inherited licensing
layer, and it is the reason this model is a better Tier-2 base than one without FUCCI
observables at all.

`FUCCI_THRESHOLD` is now 0.02 and counts as **one declared fitted parameter**;
`PUBLISHED_FUCCI_THRESHOLD = 0.05` is kept and pinned so nobody restores it and
re-derives the retracted conclusion.

### REJECTED: step 1, coupling Cdt1/geminin synthesis to E2F

`ks_CDT1_E2F` and `ks_Geminin_E2F` are constants despite their names — a real naming
defect. Connecting them to E2F was the planned step 1. It was implemented, measured and
**rejected**: it makes FUCCI structure worse. Each variant scored at *its own* best
cutoff, so the comparison is fair:

| variant | error |
|---|---|
| **published (off)** | **0.0006** |
| gem 1.0 | 0.4663 |
| cdt1 1.0 | 0.0973 |
| gem + cdt1 1.0 | 0.1524 |
| gem + cdt1 0.5 | 0.0402 |

Mechanism of the failure: geminin's peak moves from +2.0 h (mitosis, correct for mAG) to
+16.5 h (late G1) and its amplitude falls 0.133 → 0.043, because E2F is a sharp late-G1
spike while geminin needs to accumulate across S/G2/M. The constants happen to produce
the right dynamics; the defect is in the naming, not the behaviour.

`kd_CDT1_CDK2` was also tried — clearing Cdt1 at S onset via CDK2, which is real biology
(CDK2 phosphorylation licenses SCF-Skp2; Cdt1's only route here is SCF, which peaks at
mitosis rather than G1/S). Too blunt: `CCNE_CDK2 + CCNA_CDK2` never falls below ~0.07, so
it acts as near-constant degradation and erases Cdt1. Doing it properly needs an S-phase
marker the inherited model does not have — there is no DNA replication variable. Left
wired and defaulted off.

Both negative results are kept as tests rather than deleted, following this project's
practice of committing rejected fixes.

### Step 2: the E2F sub-family split

`E2F6`, `E2F7`, `E2F8` added as states; the inherited `E2F` is the activator pool.
E2F7/E2F8 are induced by E2F and repress it back — the canonical delayed negative
feedback — through a single saturating denominator rather than a product, since they
compete for the same promoters. E2F8 additionally carries CycA-driven clearance, which
is what makes it G1/S-restricted rather than merely short-lived (Baniol: `E2f8~Ccna2`
= −0.56). E2F6 sits on its own branch, maturation-driven from step 4.

Baniol's correlations are scored as **signs, not magnitudes**, deliberately. Those are
Spearman correlations across single cells at Smart-seq2 depth, where dropout attenuates
|r| substantially; a deterministic trajectory has no measurement error and will always
give larger |r|. Comparing magnitudes without an attenuation correction would be scoring
the noise model rather than the biology. All three signs are reproduced.

**The result worth having.** Repression strength lengthens the period monotonically, so
in-silico knockdown *shortens* it — the direction the lab's own data shows (cycling
cardiomyocytes at P7: KO 31.6 % vs WT 25.6 %). Nothing was fitted to the KO data:

| genotype | period | vs WT |
|---|---|---|
| WT (both repressors) | 39.33 h | 1.00× |
| E2f7 KD | 29.42 h | 0.75× |
| E2f8 KD | 36.96 h | 0.94× |
| **E2f7/E2f8 double KD** | **28.11 h** | **0.71×** |

The double lands exactly on the published period, as it must. And the epistasis is real:
singles give −9.91 h and −2.37 h (sum −12.28) against a double of −11.22 h, i.e.
sub-additive. Tier 1 also found genuine epistasis on Ect2, which matters because the
lab's data *is* a double knockout.

A tension to carry into Phase 3: strong enough repression to reach Baniol's correlation
magnitudes roughly doubles the cycle period. Repression strength is therefore a Phase 3
calibration target, not something to pick by eye now — the enable parameters stay at zero.

### Step 3: the cytokinesis arm

Six new species, all absent from the inherited model:

```
                 Ect2 ──┐
AurKB ── Centralspindlin┴─→ RhoA ──┐
                                   ├─→ Midbody ──→ abscission
                        Anillin ───┘
```

**RhoA is obligatory by construction, not by parameter choice.** `d.Midbody` has exactly
one production term and RhoA is a factor in it, so no setting of any parameter can make a
midbody without RhoA — asserted by a test that greps the source and counts the term.
Tier 1 hit the opposite arrangement four separate times; the worst OR'd a second
`(Centralspindlin & AurKB)` route onto the same node, which bypassed the Ect2/RhoA arm
entirely and made it *structurally impossible* for Ect2 to be rate-limiting. Removing it
widened the arm's hiPSC-to-P1 span from 1.6× to 42×. Here Centralspindlin and AurKB act
**upstream, through RhoA**, never in parallel with it.

Timing comes from CDK1: Ect2's GEF activity and centralspindlin bundling are both blocked
by CDK1 phosphorylation, so the arm is held off until MPF is destroyed. That is what makes
the midbody a late-mitotic structure rather than something accumulating through G2.

Ect2 is E2F-driven and E2F8-repressed — Tier 1's `r071` is `E2Fact & !Maturation & !E2F8`,
with the `!Maturation` arm arriving in step 4. This is where step 2 earns its keep: the
repressor Ect2 needs already exists.

#### Ect2 is rate-limiting, and only for division

Tier 1's central claim, now mechanistic. Knocking Ect2 down converts **every** division to
binucleation while leaving S-phase entry and mitotic entry completely untouched — a model
where Ect2 knockdown also reduced entry would be describing general toxicity, not a
cytokinesis-specific block:

| condition | S-entries | mitoses | Division | Binucleation | midbody max |
|---|---|---|---|---|---|
| WT (arm on) | 7 | 7 | **6** | 0 | 0.0868 |
| Ect2 × 0.5 | 7 | 7 | 0 | **6** | 0.0472 |
| Ect2 × 0.1 | 7 | 7 | 0 | **6** | 0.0101 |
| Ect2 × 0 | 7 | 7 | 0 | **6** | 0.0000 |
| RhoA KD | 7 | 7 | 0 | **6** | 0.0000 |
| Anillin KD | 7 | 7 | 0 | **6** | 0.0000 |
| AurKB KD | 7 | 7 | 0 | **6** | 0.0000 |

The midbody scales smoothly with Ect2 (0.0868 → 0.0472 → 0.0293 → 0.0101 → 0), so the
*fate* is binary for one deterministic cell but the *mechanism* is graded. That is what
will produce graded fate fractions once Phase 3 runs a heterogeneous population — a single
trajectory can only ever return one fate.

#### Caveats to carry forward

- **The abscission threshold does a lot of work.** The midbody peaks at 0.0868 against a
  0.05 cutoff — a 1.7× margin — so the division/binucleation call is sensitive to it. This
  is the same exposure Tier 1 records for its `r080` switch (`n = 2.0, EC50 = 0.090`),
  which `../MODEL.md` calls "the model's largest single sensitivity". It needs the same
  treatment Tier 1 gave it: profile the threshold and show the conclusions are invariant
  across a range. Phase 3.
- **No back-coupling into the core oscillator.** The arm reads the cycle but does not drive
  it, so a cell that fails abscission keeps the same molecular oscillation and simply ends
  up binucleate. That keeps the reduction property trivially exact, but it means ploidy
  does not feed back on the cycle — and real polyploid cardiomyocytes cycle differently.
  A Phase 4 question, not a permanent claim.
- Phase 1's `:MitoticCompletion` label is preserved when the arm is off. Reporting
  `:Binucleation` merely because no midbody formed would be a guess dressed as a result.

### Step 4: the maturation axis

`M` is a parameter, not a state — a developmental coordinate that is effectively constant
over one cell cycle. The six named contexts are **read from Tier 1's manifest**
(`../cmcycle/data/cmfate_model.toml`) rather than duplicated, so the two tiers cannot
drift on what "P7 mouse" means.

#### The couplings cost one parameter, not two

Baniol measured `Ect2 ~ M` at **−0.563** and `E2f6 ~ M` at **+0.396** within cycling
ventricular cardiomyocytes (n = 89, P0 and P7 pooled so stage cannot be doing the work).
But `M` is `mean z(FAO) − mean z(glycolysis)`, z-scored *within* one 285-cell dataset, and
`../MODEL.md` is explicit that it "has no absolute cross-system scale".

A correlation between two z-scored quantities is a regression slope in z-units. So what
the measurement actually fixes is the **sign** of each coupling and their **ratio**,
−0.563 / +0.396 = **−1.422**, which is scale-free and therefore transfers even though
neither slope does. It does *not* fix absolute strength.

The honest parameterisation is therefore one shared `maturation_gain` with the ratio
welded in, not two independent slopes — one fitted parameter for the whole axis instead of
two. That is the difference between using the measurement and merely citing it. Ect2 is
suppressed through a saturating denominator so it cannot go negative at `adult`'s M = 0.95.

#### Maturation closes the abscission arm

Tier 1's mechanism, reproduced. Ect2 and the midbody fall monotonically with M and the
fate switches inside the observed range:

| context | M | Division | Binucleation | Ect2 | midbody | period |
|---|---|---|---|---|---|---|
| `hipsc_cm` | 0.12 | **6** | 0 | 0.434 | 0.0742 | 39.33 |
| `mouse_p0_invivo` | 0.30 | **6** | 0 | 0.396 | 0.0609 | 39.33 |
| `mouse_p1_invivo` | 0.50 | **6** | 0 | 0.362 | 0.0508 | 39.32 |
| `mouse_p7_invivo` | 0.55 | 0 | **6** | 0.354 | 0.0487 | 39.34 |
| `adult` | 0.95 | 0 | **6** | 0.302 | 0.0369 | 39.33 |

The period does not move with M (spread 0.05 %, and the two contexts sharing M = 0.50 give
bit-identical periods) — the arm reads the cycle without driving it, as designed.

**`mncm_invitro` and `mouse_p1_invivo` are currently indistinguishable**: both sit at
M = 0.50 and differ only in `InVitro`, `ROSenv` and the adrenergic inputs, which belong to
the signalling layer in step 5. Tier 1 separates them (Polyploidization vs Binucleation)
through the ROS→DDR arm, so Tier 2 cannot reproduce the full clonidine triad until step 5
lands. Stated rather than glossed.

#### Two knobs, deliberately not one preset

`MATURATION_ON` carries only the Ect2 coupling — the **fate** knob. `E2F6_EXIT_ON` carries
E2F6, which Tier 1 calls the "cell-cycle exit enforcer" and which acts on the **period**:
at M = 0.55 the cycle runs 39.3 h at `ks_E2F6 = 0`, 48.1 h at 0.01, 57.7 h at 0.02, 85.6 h
at 0.05. They calibrate against different data — fate fractions versus cycling fractions —
so folding them together would tie two independent Phase 3 targets to one number. Measured
with both on at once at their first-guess values: everything binucleates and the period
stretches to 100 h, because the knobs fight.

#### A real bug the step surfaced

Turning E2F6 on drove the solver to `MaxIters` at t = 1518 after 2×10⁶ steps. Cause: total
geminin spends 12–19 % of each cycle within ±10 % of the FUCCI cutoff, and a
`ContinuousCallback` on a level a signal lingers near re-triggers and restarts the step
every time — **995,356 logged crossings** in one run. A refractory guard cleaned the log
but not the solver cost; the fix is a Schmitt-trigger deadband on every detector, which
removes the chatter at source. Recorded times are unchanged (the trigger still fires at
the true level; only re-arming uses the deadband) and the callbacks now cost almost
nothing — 7,602 steps against 7,045 with no events at all. Pinned by regression test,
because this would have been far more painful to find inside a 10⁵-run Phase 3 ensemble.

### Step 5: oxidative stress → DDR → Ccng1 → the mitotic-entry brake

The other half of Tier 1's 2×2. Step 4 built maturation closing the abscission arm; this
is culture closing the mitotic-entry brake.

**Turning the arm on is not a refit.** The inherited model *ships* it dormant:
`kf_ATMp = 0`, with the damage input commented out in `d.ATM`, and `ks_p21 = 1e-4`, which
leaves p21 identically zero in every published figure. Those are off-switches, not
estimates — the source paper's own Figure 6 activates the arm by hand the same way.

ROS enters as an input, not a state: it is a property of the environment, constant on the
timescale of a cycle, exactly like `M`.

#### The held-out contrast works

`mouse_p1_invivo` and `mncm_invitro` sit at the **same M = 0.50** and differ only in
culture, so the pair isolates the ROS arm. Tier 1 fitted neither and Tier 2 has nothing
fitted to either:

| context | M | ROSenv | InVitro | predicted | Tier 1 observed |
|---|---|---|---|---|---|
| `hipsc_cm` | 0.12 | 0.25 | 1.0 | Polyploidization | Division ❌ |
| `mouse_p1_invivo` | 0.50 | 0.20 | 0.0 | **Binucleation** | Binucleation ✓ |
| `mncm_invitro` | 0.50 | 0.25 | 1.0 | **Polyploidization** | Polyploidization ✓ |
| `adult` | 0.95 | 1.00 | 0.0 | Polyploidization | — |

And the *mechanism* is right, not just the label: in culture, NEB never fires while
S-phase entry continues. Polyploidization here is genuinely S-without-mitosis, not simply
fewer cycles. The brake acts through MPF — raising it drives `CCNB_CDK1` down, `LMNAp`
follows, and mitotic entry stops exactly when `LMNAp` can no longer reach its threshold.

#### KNOWN MISS: hiPSC-CM

Predicted polyploid, observed dividing. Recorded, not tuned away.

`hipsc_cm` and `mncm_invitro` carry **identical** oxidative input (ROSenv 0.25,
InVitro 1.0), so the DDR brake hits them equally and only `M` distinguishes them — and `M`
acts on the abscission arm, not on mitotic entry. Tier 1 gets hipsc right because hipsc
*is* its calibration context: its MitoticEntry gate is fitted to 0.7698 there, and
`../MODEL.md` says plainly that "the fate layer is fitted exactly and predicts nothing" at
that context. So this is a genuine prediction failure, and it says the model lacks
whatever keeps the DDR response weak in immature cardiomyocytes.

**Not fixed by gating the DDR arm on `!Maturation`** — that is one free parameter fitted
to one outcome, and Tier 1 tested the directly analogous move (gating the clonidine
response on `!Maturation`) and rejected it: mean fold error went 26 % → 55 %, worse than
doing nothing.

#### Two things declared rather than buried

- **A parameter degeneracy.** Only the product `ks_Ccng1_p53 × kf_CCNB_Ccng1 / kd_Ccng1`
  sets the brake strength — Ccng1 is at quasi-steady state, so halving synthesis and
  doubling effect is the same model. Verified by test. Phase 3 must fit **one** effective
  parameter here, not three.
- **`p27` is deferred.** It would cost three species and ~8 parameters (mirroring p21's
  complexes), none of which any available measurement constrains, and it is not required
  for the triad. Tier 1 can afford it because a logic node is free; in a mass-action model
  it is not. Adding eleven unconstrained quantities for no measurable gain is exactly what
  the constraint budget exists to prevent.

#### The p53 mass-balance defect, now fixed

`d.p53p` gained p53p at a rate independent of available p53, because the `*p53` factor
that `d.p53`'s matching loss term carries is missing. Harmless in the published figures —
`Chk2p` is identically zero there, which is why nobody saw it — but not harmless with the
arm live. `fix_p53_massbalance` defaults **off** so the reduction property stays
*unconditional* (bit-exact for any state, not merely reachable ones) and ships with
`DDR_ON`, which is exactly when the defect can bite.

#### A second context-loading bug caught here

`context_params` initially read an unlisted input as 0. But Tier 1's manifest declares
`default_on = [... "ROSenv" ...]`, so an unlisted default-on input sits at **1.0**.
`adult` names no `ROSenv`, so the adult heart was getting *no* oxidative stress — and
produced adult cycling faster than P0 (39.3 h against 52.2 h), which is backwards. Puente
2014, the source of this arm, is precisely about postnatal ROS rising. Fixed and pinned.

## Next: Phase 3 — calibration and the population layer

Phase 2 is complete. Remaining from its original list, deferred with reasons:

3. **The cytokinesis arm** — Ect2 → RhoA → Centralspindlin/AurKB/Anillin → Midbody, with
   RhoA obligatory (Tier 1 learned that the hard way four times: an OR'd bypass made it
   structurally impossible for Ect2 to be rate-limiting, which is the model's central
   claim). This is what splits `:MitoticCompletion` into Division and Binucleation.
4. **The maturation axis `M`**, with Ect2 and E2f6 coupled to it by the measured
   correlations (r = −0.563 and +0.396) rather than by free parameters.
5. **`Ccng1`, `p27`**, and turning the DDR arm on — which requires fixing the `d.p53p`
   mass-balance defect first. Note `p21` is currently inert (`ks_p21 = 1e-4`, and p21 is
   identically zero in the default configuration), so it needs activating before it can
   act as a brake.

**Gate:** with the new modules disabled, every state must match the Phase 0 golden
trajectory exactly — enforced by test, not by inspection. That is what lets the paper
claim the cardiomyocyte model *contains* the published one.
