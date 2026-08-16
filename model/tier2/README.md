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

## Status: Phase 0 complete

`src/inherited/` is a **byte-identical** copy of the published `model_files/`
(`state.jl`, `parameters.jl`, `diff_eqns.jl`). Do not edit it — extensions belong in
their own files so the diff against the published model stays legible to a reviewer.

52 tests pass. The source repo has no test suite, so these are the first the inherited
model has had.

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

## Layout

```
Project.toml        name/uuid/[compat] — the source repo's env declared none
Manifest.toml       pinned to the published model's resolved stack (see above)
src/
  CmTier2.jl        module: solve_baseline, solve_drug, doubling_time, peak_period
  inherited/        byte-identical copy of the published model_files/ — do not edit
test/runtests.jl    52 tests: the Phase 0 gate
scripts/            (Phase 1+)
```

## Next: Phase 1

Events for restriction-point passage, nuclear-envelope breakdown, anaphase and mitotic
exit, off states that already exist (`ppRB`, `LMNAp`, `PTTG1`, `CCNB_CDK1`); nuclei /
ploidy / cell counters; FUCCI state read-out. Gate: simulated S/G2/M durations against
the five `duration` rows in [../cmcycle/data/cmcycle_targets.csv](../cmcycle/data/cmcycle_targets.csv)
— a result Tier 1 structurally cannot produce.
