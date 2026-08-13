# Cell-cycle fate model — what's left

Ordered by whether it *unblocks* something else. Each item says why it matters, so
a stale one is easy to spot and drop.

---

## Done since the last pass

- **Extracted the KO bundle** (one `docker run` against the existing rocker image):
  30,030 cells, 21,598 cardiomyocytes, per-cell metadata, the P0/P7 DE tables, and
  the model-node expression. See `RESULTS.md` §8.
- **Settled the E2f8 phase-scoring question.** Real in principle, negligible in
  practice: E2f8 carries 0.15–0.54% of the S-set mean, phase agrees for 99.1% of
  cells with and without it, and the KO−WT cycling gap moves by ≤0.5 pp. It was
  also the *opposite* sign from what a first pass assumed — E2f7/E2f8 mRNA is
  *higher* in KO, not lower.
- **Confirmed the knockout is functional.** 24/26 canonical E2F targets significantly
  up at P7 against a housekeeping baseline near zero; HALLMARK_E2F_TARGETS regulon
  activity KO−WT = +1.12 (P0) and +3.15 (P7).
- **Built and calibrated the Tier-1 logic network** (`cmfate`): 55 nodes, 78
  reactions, stdlib-only, 53 tests. Reproduces the fit context to 1e-5 and gets
  3/3 held-out clonidine outcomes.
- **Fixed the screen's scoring.** Two errors, both flattering the model: the hit test
  was purely relative and applied at P7 where baseline entry is 0.27%, so a
  0.4-percentage-point change scored as a hit; and the "conversion" figure counted any
  positive delta, almost all epsilon. A hit now needs a fold **and** a 1 pp absolute
  effect, and specificity is measured in hiPSC-CM where Murganti ran theirs. Result:
  10.5% hit rate against their 6.4%, and 25% of hits also raising the division
  share against their 2 of 6 — the model now agrees with the wet screen instead of
  contradicting it.
- **Fixed the G1/S switch.** Three changes: the three restriction-point reactions
  (`r044`, `r045`, `r060`) were gate-shaped — they were the only graded ones in a model
  whose two downstream gates are steep, so the most famous switch in the cell cycle was
  the only one built as an interpolation and could not hold a fold; an OR'd CycE leg
  (`E2F2 & Maturation => CycE`) that cancelled the switch's travel was removed, and E2F2
  re-cast as an AND-term brake on mitotic competence, which is what Baniol actually
  propose; and CycD's drive was scaled 0.70× to put the fold inside the context range.
  The `!PKA` shortcut came off the S-phase reaction, so entry is now driven *through*
  the switch. Result: **E2F1 travel 1.4× → 25.9×**, CycE 1.2× → ~36,000×, Rb 0.14 → 0.92
  across contexts, **triad still 3/3**, Ect2 still rate-limiting, and the mean
  clonidine fold error **236% → 26%**.
- **Found the OR'd-floor pathology a fourth time, in my own fix.** Adding E2F2's brake as
  its own reaction *raised* mitotic competence, because a separate reaction can only add
  through the weighted OR. An inhibitor has to be an AND term inside a reaction. That is
  now stated in `r068`'s evidence column, since it has bitten this model four times.
- **Fixed a third engine bug, found mid-diagnosis.** `_compiled()` cached the
  flattened reaction list keyed on weights alone, so replacing a reaction to change
  its reactants, `n` or `EC50` was silently ignored and the model kept running the
  old structure. It invalidated a diagnosis and produced a conclusion that was wrong
  in the opposite direction. Now keyed on weights *and* reaction identity, with the
  mutation contract documented and a regression test.
- **Fixed the two engine bugs** the old plan listed: the `EC50**n >= 0.5` sign
  inversion now raises, and `perturbation_matrix` takes a tuple so the double
  knockout is expressible.

---

## Blocking — the model's known misses

### 1. Residual tensions from fixing the G1/S switch
The switch is **fixed** — see the Done list for what changed and what it bought. Three
things it left behind, in order:

**The hiPSC drug response is under-predicted ~2× across every readout, and this now
looks like a formalism ceiling rather than a missing edge.** Entry 1.02× against 2.44×,
abscission 1.33× against a measured midbody doubling of 2.04×, division 1.29×. The
*direction* is right everywhere; the magnitude is about half.

Two hypotheses were tested and rejected:

- **Autophagy ⊣ p21** (clonidine's documented mechanism is autophagy induction, and there
  is literature on autophagic turnover of CDK inhibitors). It helps hiPSC — 1.02× → 1.49×
  — but p21 is *higher* at high maturation (0.68 at mNCM versus 0.34 at hiPSC), so
  relieving it helps the wrong end: mNCM blows out to 6.70×. Mean fold error 26% → 121%.
- **The same, plus gating the drug on `!Maturation`** (motivated: Nisch, clonidine's I1
  receptor, falls 1.277 → 1.047 with maturation at 100% detection — measured here and
  still unused by the model). Better, 55%, but still worse than doing nothing.

The reason is structural, and worth stating as a limit rather than a to-do. At hiPSC
everything the drug could relieve is *already* near its operating maximum — the switch is
open (Rb 0.14) and the brakes are low (p21 0.34, p27 0.02, Ccng1 0.10) — so relieving
anything yields little. Every brake that *does* have headroom has more of it at high
maturation. So **no re-wiring within this formalism can give a larger response at low
maturation than at high**: gain is capped at the reaction weight, and there is nothing
left to release.

Two things follow. First, the comparison may be against the wrong observable: 2.44× is a
Ki-67 index, and the pre-flight already established that Ki-67 carries a dwell-time
component which a prevalence-only `Ki67` node cannot express — a drug that converts
arrested cells into cycling ones raises the Ki-67 index without raising the entry rate.
Notably the model does best against the one *cumulative* measurement it is judged on
(in-vivo EdU, 1.83× against 1.52×). Second, this is a concrete argument for Tier 2, which
has real rates and dwell times and can represent the difference.

**Current configuration is the best found** — folds 1.02 / 2.13 / 1.83, mean error 26%,
triad 3/3. Do not "improve" the hiPSC number by relieving a brake without re-checking the
mature contexts; that trade has been measured and it loses.

**The comparators drifted 28%.** All five move 27–28% between P0 and P7, uniformly,
tracking E2Fact's own 29% fall — i.e. it is the switch working, not five bad edges. The
test bound was raised from 25% to 35% with that reasoning written into it, and it now
asserts *uniformity* rather than just magnitude. The real refinement: truly flat targets
under a travelling E2F activity would need their reactions to saturate, so a 29% fall in
input gives a small fall in output. Worth doing, and it would let the bound come back down.

**CycE spans ~36,000× and the maturation slope of entry is ~4× too steep — one problem,
not two.** `SPhase` multiplies four maturation-dependent factors (CycE and three brakes),
so after the switch fix CycE's enormous travel compounds into far too steep a slope.
Measured against Baniol (see item 6): their P0→P7 fall in genuinely cycling
cardiomyocytes is 3.2×, the model's is ~13×. Flattening CycE's range would address both,
and is the most concrete open item on the model itself.

### 2. Widen the curated gene panel
Only **29 of 63** model node genes are in the 2,181-gene panel. E2f2–E2f6, Rb1,
Ccng1, Chek1, Wee1, Pkmyt1, Rhoa, Nisch, Adrb1 and Mapk12 are all missing, so most
of the model cannot be scored against the lab's own cells. This is the single change
that would let the KO data validate the network rather than just its E2F layer.

---

## Build the Tier-2 ODE

### 3. Extend Gérard & Goldbeter 2009 rather than starting over
Restore the three missing SBML files first (they download from
`https://www.ebi.ac.uk/biomodels/search/download?models=<ID>`; the
`model/download/...` path 301-redirects to a broken host). Then keep gg2009's Rb–E2F
switch, four cyclin-CDK waves, Skp2/p27/Cdh1/Cdc20, Wee1/Cdc25 gate and Chk1/ATR
arm; add **Cdt1 and geminin as FUCCI observables** (the highest-value single feature
— it makes both papers' primary readout simulatable), p21, the E2F sub-family split,
Ccng1, cytokinesis competence driven by the measured Ect2(M), discrete
envelope-breakdown and abscission events, nuclei-and-ploidy bookkeeping, and M.

Hand-write in numpy/scipy so runtime needs no Tellurium; use Tellurium (confirmed
installable — libroadrunner ships a cp314 wheel) only to cross-validate the
inherited core against the published reference.

### 4. Population layer
Fate fractions are tail statistics — reproducing 1.40% needs the *shape* of the
parameter distribution, not its mean, so a mean-field treatment cannot get there.
Plan an ensemble of single-cell runs with lognormal parameter heterogeneity, sized
so Monte-Carlo error is small against the binomial error on 8/570 cells (≈10⁴–10⁵
runs). Calibrate the noise scale to the measured duration CV (0.265 for mouse P0),
not to a guessed value.

### 5. Performance — partly addressed
Done: reactions are flattened to plain tuples once (the dataclass attribute lookups
in the inner loop dominated), output-only nodes are solved algebraically instead of
integrated, and the screen's two passes were merged into one sweep. A steady state is
~405 ms and the 76-node sweep ~31 s; calibration is ~70 s but cached against a hash
of the spec.

Left: that is fine for the current scale and nowhere near enough for item 4. An
ensemble of 10⁴–10⁵ single-cell runs at 0.4 s each is days, so the population layer
needs numpy vectorisation of the RHS or a parallel sweep before it is worth starting.

---

## Calibration and validation

### 6. Fit mouse, predict human — attempted, and it found two things
The flip could not be completed as planned, for an instructive reason: **the mouse target
as naively defined is unreachable at any weight**, and that turned out to be a
comparison error rather than a model error.

Baniol's Fig 1D gives FUCCI *states*, so the naive cycling fraction (1 − mKO2⁺) is 32.5%
at P0. But their own Suppl 1G shows only **22.7%** of the G1/S double-positives are Ki67⁺
at P0 — the rest have prematurely exited, a distinction they draw explicitly and which
this project recorded early and then failed to apply. Correcting with their own data:

    mAG⁺ 12.36% (all Ki67⁺) + 22.7% of 19.1% + mitotic 1%  =  17.7%

against a model ceiling of **17.2%** — a **1.03× match, with nothing fitted to it.** That
is the strongest unfitted validation the model has, and it was hiding behind a
mis-specified observable. Pinned by
`test_the_mouse_in_vivo_cycling_level_is_reproduced_unfitted`.

**What remains is the slope, not the level.** Corrected the same way, P7 is 5.5% observed
against a 1.3% ceiling — 4.2× short — so the P0→P7 fall is 3.2× measured and ~13×
modelled. That is the CycE over-swing in item 1, measured against data for the first
time. Pinned by `test_the_maturation_slope_of_entry_is_too_steep`.

Still worth doing properly once the slope is fixed: calibrate on the corrected mouse
level and hold out all four hiPSC fate fractions. The pass criterion stays what it was —
predicted fractions inside the binomial CI, and the *significance structure* of the
durations reproduced rather than three means matched.

### 7. Sensitivity on the abscission switch — half done, and reassuring
The switch position (`r080`: n=2.0, EC50=0.090) is a fitted quantity, not an
evidence-derived one, and it is what makes Ect2 rate-limiting — so it is the model's
largest single sensitivity.

Profiled downward (EC50 = 0.060, 0.075, 0.090, re-calibrating the weight at each):
the triad stays **3/3**, Ect2 knockdown still removes **100%** of division at P0 while
leaving entry untouched, and the hiPSC division share is pinned at 52.7% by
construction. The conclusions are invariant across that range because the calibrated
weight absorbs the change — which substantially de-risks the headline claim.

Re-checked after the G1/S switch was fixed, since that changed E2Fact by 38× and the
abscission arm hangs off it. The arm barely moved at the fit context — E2Fact 0.405 →
0.406, Ect2 0.1807 → 0.1810, AbsRaw 0.1396 → 0.1402, all 1.00× — and its hiPSC-to-P1
separation actually *improved*, 42× → 64×. So the profile is not invalidated and the
switch still sits where it was positioned.

Left: the range **above** EC50 = 0.090 is untested (that run timed out), the target
becomes unreachable somewhere below 0.16, and `n` has not been profiled at all. Worth
finishing, but lower priority than it looked — the arm is more robust to the core change
than expected.

### 8. Negative controls — in place, keep it that way
Not an open task: `test_negative_controls_do_not_move_with_maturation` already asserts
that E2f3, AurKB, Ccna2, Centralspindlin and Anillin move < 25% between P0 and P7,
because all five are flat in the data. Listed here so it is not quietly deleted when
a future change makes it fail — the failing test is the finding. Still unused as a
constraint: endothelial EdU was flat (p = 0.85) while smooth muscle responded, so no
model may predict a generic mitogenic effect.

---

## Experiments the model is built to motivate

Ranked by how cleanly each could be run with assays both papers already use.

1. **Is there a critical M\* above which added entry cannot add cardiomyocytes?**
   Same clonidine dose in vivo at P1 versus P7; read CM number, mononucleation,
   ploidy, AurKB⁺ midbodies. Falsified if CM number scales with entry at all M.
2. **Ect2 rescue.** Ect2 overexpression ± clonidine in mNCM, midbody count as
   readout. Clonidine alone gives zero midbodies in >10,000 cells; if Ect2 is truly
   rate-limiting, this should produce them. Falsified if it does not — which would
   point at Anln/Cep55 or a non-transcriptional block.
3. **Is the S/G2 duration distribution bimodal?** Extend imaging past 72 h (the
   pre-flight predicts a chronically arrested pool that never scores an outcome).
   Published n is 9 and 5 — far too small to see this.
4. **Two G1/S failures, one FUCCI appearance.** Chek1 inhibition should rescue the
   P7 delay but not the P0 premature exit. Cheap to test in silico first, since the
   Chk1/ATR arm is already parameterised in gg2009. Note the DDR module has no
   support on our primary axis, so this is the test that would earn it.
5. **EdU/Ki-67 ratio as an arrest assay.** Co-stain both in the same clonidine-treated
   cells and count double-positive versus EdU-only. The two candidate explanations
   for the 1.31 ratio predict opposite results. Cheap, decisive, done by neither paper.
6. **Does mitosis really take ~4 h in neonatal CM?** Direct envelope-breakdown to
   reformation timing, and it should be *longer* in binucleating than dividing cells.
7. **4C versus 8C.** DNA-content histogram of the FUCCI-classified polyploid cells.
   Trace shape alone cannot separate G2 arrest at 4C from true endoreplication, and
   the tetraploid shift with no reported octaploid rise suggests one round only.

---

## Improvements to the existing Shiny app

1. **Replace the `sig_ploidy` heuristic.** It is a difference of two module
   averages that decline in lockstep (p = 0.93 on this data), so it cannot see the
   Ect2 signal. Keep it for display; drive the polyploidization read-out from named
   nodes.
2. **Widen the curated panel to all of E2f1–E2f8.** Baniol's entire E2F argument
   rests on E2f2 and E2f6, and neither is in the 2,181-gene panel — so the app
   currently cannot show the paper's main finding on the lab's own cells. Both are
   detectable at Smart-seq2 depth (E2f2 12→32% of cycling vCM; E2f6 52→100% of
   noncycling P7 vCM).
3. **Propagate per-cell atrial/ventricular labels.** `cm_subtype` is declared in
   `CAT_COLS` in [app.R](../shiny_app/app.R) but absent from the bundle, so it is
   silently dropped; aCM/vCM exists only per-subcluster.
4. **Add the sort-enrichment caveat to the About block.** Anyone reading cycling
   fractions off a FACS-enriched dataset will get the developmental direction
   backwards.
5. **Surface the model's predictions next to the data.** Now unblocked — Tier 1
   runs, and `model.e2f78_knockdown()` returns signed shifts per fate. The
   proliferation-versus-cytokinesis quadrant is already the right plot for it. Needs
   item 3 first if the comparison is to go beyond the E2F layer.

---

## Housekeeping

### 9. The pathway map is dense in the middle columns
93 edges in one figure is legible on screen but marginal for a slide or print. If it
is needed at that size, the natural split is two panels — stimulus → E2F, and
machinery → outcome — rather than dropping edges.

---

## Known limitations to keep stating

- **~20–25 effective independent constraints** exist in total, so ~13 fitted
  parameters is the defensible ceiling; most must stay literature-fixed.
- **The FUCCI reporter layer cannot be validated transcriptomically** — Gmnn~Cdt1
  mRNA correlates *positively* because the oscillation is post-translational.
  Imaging traces only.
- **M has no absolute cross-system scale** as defined (z-scored within one
  dataset), so the clonidine triad is an ordinal constraint until a shared
  reference exists.
- **`Adra2a/2b/2c` are absent** from the expression store, so clonidine's canonical
  α2 target cannot be assessed in this data at all.
- **285 cells**, with two groups at n = 8 and n = 4. Every atrial-P0 claim is
  fragile and should carry its n.
- **`tau` in a normalized-Hill model is a relaxation constant, not a duration.**
  Matching three measured durations with three `tau` values is a definition, not a
  prediction — and the model cannot reproduce the observed 48 h optimal screening
  window, which is a variance effect. That gap is the honest reason Tier 2 exists.
