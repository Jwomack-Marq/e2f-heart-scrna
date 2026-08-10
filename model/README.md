# `model/` — cardiomyocyte cell-cycle fate model

An ODE/logic model of the four-way cardiomyocyte cell-cycle outcome — quiescence,
productive division, binucleation, nuclear polyploidization — calibrated against
two FUCCI papers, plus the re-analysis of their data that constrains it.

Spatial components are deliberately out of scope.

| Document | What it is |
| --- | --- |
| [RESULTS.md](RESULTS.md) | What we did, what came out, with figures |
| [TODO.md](TODO.md) | What is left, and what would improve it |

**Papers.** Baniol et al. 2021, *Exp Cell Res* 408:112880 (FUCCI mouse, in vivo +
285-cell scRNA-seq, ENA PRJEB47622) · Murganti et al. 2022, *Front Cardiovasc Med*
9:840147 (TNNT2-FUCCI human iPSC-CM, live-imaging fates + clonidine in three
systems).

## Layout

```
cmcycle/
  logic.py                normalized-Hill logic-ODE engine + adaptive integrator
  spec.py                 spec loader, structural lint, gate calibration
  model.py                the calibrated model and the experiments run on it
  preflight.py            closed-form consistency checks on the calibration set
  baniol.py               re-analysis of the P0/P7 scRNA-seq
  svg.py                  dependency-free SVG chart primitives
  figures.py              regenerates figures/
  data/
    cmfate_species.csv    55 nodes, each with gene anchor and justification
    cmfate_reactions.csv  78 reactions, each with its own evidence column
    cmfate_model.toml     contexts, calibration targets, fate/gate declarations
    cmcycle_targets.csv   43 published measurements, each with its citation
tools/
  validate_palette.py     Python twin of the dataviz palette validator
  extract_ko_bundle.R     one-shot Docker extraction of the KO app_data.rds
  build_report_page.py    builds the shareable HTML report from figures/
data/ko_export/           small derived tables from the KO bundle
tests/                    51 tests
figures/                  10 generated SVGs + results.json + model_results.json
```

## The model in one paragraph

`cmfate` takes eight inputs — a maturation coordinate, an in-vitro flag, mechanical
load, adrenergic tone, two growth-factor arms, oxidative stress and a clonidine dose
— and returns the four cardiomyocyte cell-cycle fates. Those four sum to exactly 1
with no normalization step, via a complementary product over three gate nodes
(S-phase entry, mitotic entry, abscission), which also makes them **closed-form
identifiable** from Murganti's measured fate fractions: three independent
one-dimensional solves rather than a joint fit. Calibrated at the hiPSC-CM context
alone, it then gets **3 of 3** held-out clonidine outcomes — division, then
polyploidization, then binucleation — by changing only maturation and the in-vitro
flag. It over-predicts the *magnitude* of the entry response at high maturation,
which is [TODO](TODO.md) item 1.

## Running it

**No third-party dependencies.** Everything is standard library, so it runs in a
bare interpreter — which is the point of the pre-flight: an inconsistent dataset
combination should be catchable before any solver exists.

```bash
python3 -m cmcycle.preflight     # six checks; two flag real tensions
python3 -m cmcycle.baniol        # the re-analysis
python3 -m cmcycle.figures       # regenerate all 10 figures + both results.json
python3 -c "from cmcycle import spec, model
net = spec.load_calibrated(verbose=True)
for r in model.clonidine_triad(net):
    print(r['context'], r['predicted'], 'OK' if r['hit'] else 'MISS')"
```

The calibration is three bisections in pure Python (~70 s) and is cached beside the
spec, keyed on a hash of the three spec files — edit any of them and the cache is
ignored automatically, because a stale fit silently applied to a changed network is
the failure worth designing out.

Tests need `pytest` only:

```bash
uv venv --python 3.12 && uv pip install pytest
.venv/bin/python -m pytest       # 51 passed
```

### The expression data is not in this repo

`baniol.py` reads a sibling ~40 GB store. Point `CARDIAC_RNASEQ_ROOT` at the
`cardiac-rnaseq-explorer` checkout containing `data/Baniol2021_FUCCI/` and
`expr/Baniol2021_FUCCI/`; it defaults to the lab-server path. Without it the
analysis tests **skip** rather than fail, and the pre-flight still runs.

The store is a gene-major CSR dump with no header — raw little-endian `int32`
`indptr`/`indices` and `float32` `data`. `Store._check()` asserts five invariants on
construction, because a truncated or byte-swapped copy would otherwise read as
plausible.

## Two things worth knowing before using this data

**Cycling fractions here are a sorting artefact.** P0 is essentially unenriched but
P7 is 4.5–5.2× enriched, so the raw cycling fraction *rises* from P0 to P7 in the
data and *falls* in reality. Only within-group expression contrasts are valid.

**Some genes of interest helped define the labels.** The Regev/Tirosh phase lists
contain `Ect2` and `E2f8`. `baniol.circularity_audit()` splits any panel into
circular and clean; the lists are vendored in `baniol.py` so a copy that drifts
from the upstream pipeline is itself a finding. On the KO data this turned out to be
a false alarm — E2f8 carries under 0.6% of the S-set mean and moves the KO−WT cycling
gap by ≤0.5 pp — but the audit is what established that.

## Figures

Generated as inline-able SVG with no plotting library, each carrying a scoped
`<style>` so it follows the viewer's light/dark theme. Colours come from the
dataviz reference palette; categorical use is capped at three slots because only
those clear the all-pairs colour-vision and normal-vision floors in both modes.
`tools/validate_palette.py` computes that rather than asserting it:

```bash
python3 tools/validate_palette.py "#2a78d6,#eb6834,#1baf7a" --mode light --pairs all
python3 tools/validate_palette.py "#3987e5,#d95926,#199e70" --mode dark  --pairs all
```

Every figure is driven by `baniol.run()` or `preflight.run_all()`, so a figure can
never disagree with the number it claims to show.

## Relationship to the rest of the project

The Shiny app in [../shiny_app/](../shiny_app/) is descriptive; this directory is
the mechanistic half. They share no code today. The model's node vocabulary is
chosen so the app's existing proliferation-versus-cytokinesis quadrant is the
natural place to overlay its predictions.

The Tier-1 logic engine lives here rather than in the sibling `cardiac-models`
project, and was written fresh: that implementation has two bugs this model would
have tripped over — `beta` inverts sign whenever `EC50**n >= 0.5`, which
high-threshold gates hit routinely, and `perturbation_matrix` takes one node at a
time so it cannot express the E2f7/E2f8 double knockout. Both are fixed in
`logic.py`, and both have tests.

The Tier-2 mechanistic ODE is specified but not built — it extends Gérard & Goldbeter
2009 rather than starting over. See [TODO.md](TODO.md) item 4.
