# `tier2/` — the mechanistic cardiomyocyte cell-cycle fate ODE

Tier 2 of the two-tier model described in [../MODEL.md](../MODEL.md) §6 and specified in
[../TODO.md](../TODO.md) item 3. Julia, built on the published generic cell-cycle model
in [`Cell_Cycle_Model`](https://github.com/Jwomack7512-bio/Cell_Cycle_Model) (MIT).

| | Tier 1 (`../cmcycle`) | Tier 2 (here) |
|---|---|---|
| Formalism | normalized-Hill logic (Netflux) | mass-action + Michaelis–Menten ODE |
| Nodes | 55 nodes, 77 reactions | 63 species, 218 parameters |
| Time | `tau` is a relaxation constant | real time |
| Fates | product of three steady-state gate activities | *(Phase 1)* emergent from one cell's trajectory |

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

Still to add (Phase 2): the E2F sub-family split, the Ect2/RhoA/Centralspindlin/AurKB/
Anillin/Midbody cytokinesis arm, the maturation axis `M`, `Ccng1`, `p27`, and
nuclei-and-ploidy bookkeeping. The test suite asserts all of these are currently absent,
so a half-landed module cannot go unnoticed.

## Status: Phase 1 complete

**182 tests pass** (52 Phase 0, 130 Phase 1). The source repo has no test suite, so these
are the first the inherited model has had.

`src/inherited/` is a **byte-identical** copy of the published `model_files/`
(`state.jl`, `parameters.jl`, `diff_eqns.jl`). Do not edit it — extensions belong in
their own files so the diff against the published model stays legible to a reviewer.

Phase 1 adds **no biology and no parameters**; it only measures the inherited model.

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
| **S/G2/M, FUCCI-defined** | **2.06 h** | **16.38 / 17.29 / 24.50 h** (Murganti Fig 2E) ✗ |
| mitosis (NEB → exit) | 1.74 h | ~4.2 h for CM (Tier 1 preflight) |
| cells after 14 cycles | 2¹⁴ = 16384 | — |

Two S/G2/M durations are carried deliberately. Murganti Fig 2E times a FUCCI trace —
mAG (geminin) appearing, to division — so scoring a cyclin-defined duration against it
would be a category error, and not a small one: the two differ by ~17 h here.

### The Phase 1 gate found a real limitation

The cyclin-defined duration (18.95 h) is squarely in the published range for a model
with nothing fitted to it. The **FUCCI-defined one is not**, and the reason is specific:

- Total geminin exceeds the published 0.05 cutoff for only **12.6 %** of the cycle; the
  mAG reporter marks S/G2/M, which should be ~40 %.
- Cdt1 and geminin **never overlap**, so the G1/S double-positive state has frequency
  exactly **zero**. Murganti Fig 1C reports 1.6 % and Baniol Fig 1D reports 19.1 % at P0
  — and the double-positive is precisely the population Baniol's Suppl 1G correction
  operates on, which is what gave Tier 1 its 1.03× unfitted validation.
- 36 % of the cycle reads double-negative, far more than a real FUCCI trace.

So the inherited model **has** the FUCCI observables but does not reproduce FUCCI phase
structure. The likely root cause is the Phase 0 defect: `ks_CDT1_E2F` and
`ks_Geminin_E2F` are constants, so licensing is decoupled from the cycle and geminin
peaks at mitosis instead of accumulating through S and G2.

This is recorded as **deliberately-failing tests** (`KNOWN FAILURE: the inherited FUCCI
layer has no G1/S state`), following Tier 1's
`test_the_maturation_slope_of_entry_is_too_steep`. They fail the day Phase 2 fixes this,
which is the point — the fix has to be deliberate and recorded, not silent.

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
test/runtests.jl    182 tests: the Phase 0 and Phase 1 gates
scripts/            (Phase 3+)
```

## Next: Phase 2 — the cardiomyocyte modules

In rough dependency order:

1. **Connect Cdt1/geminin to E2F**, fixing the FUCCI phase structure Phase 1 measured as
   broken. Highest value: it makes both papers' primary readout simulatable, and it is
   the prerequisite for scoring the ten `fucci_fraction` targets.
2. **The E2F sub-family split** — E2Fact / E2F6 / E2F7 / E2F8. Six consumer lines to
   rewire, asserted by test so none is missed.
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
