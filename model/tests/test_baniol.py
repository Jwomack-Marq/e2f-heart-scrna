"""Regression tests on the Baniol re-analysis.

These pin the numbers quoted in ``model/RESULTS.md``, so the write-up cannot
silently drift from the data. They need the sibling expression store; when that
is absent (CI, a fresh clone) they skip rather than fail.
"""
from __future__ import annotations

import math
import os

import pytest

from cmcycle import baniol

pytestmark = pytest.mark.skipif(
    not os.path.isdir(os.path.join(baniol._root(), "expr", baniol.DATASET)),
    reason=f"expression store not found under {baniol._root()!r}; "
           f"set CARDIAC_RNASEQ_ROOT to run these",
)


@pytest.fixture(scope="module")
def store():
    return baniol.Store()


@pytest.fixture(scope="module")
def res(store):
    return baniol.run(store)


def test_store_integrity(store):
    """Store()._check() asserts the five invariants; this pins the shape."""
    assert store.info == {"n_genes": 20754, "n_cells": 285, "nnz": 2042448,
                          "source": "ingest.py"}
    assert all(store.has(f"E2f{i}") for i in range(1, 9))


def test_group_sizes(res):
    assert res["group_sizes"] == {
        "P0-aCM-cycling": 8, "P0-aCM-noncycling": 4,
        "P0-vCM-cycling": 26, "P0-vCM-noncycling": 67,
        "P7-aCM-cycling": 18, "P7-aCM-noncycling": 64,
        "P7-vCM-cycling": 63, "P7-vCM-noncycling": 35,
    }


def test_ect2_is_the_specific_outlier(res):
    """0.55x, clears Bonferroni, and sits in the genome-wide lower tail."""
    e = res["per_gene"]["Ect2"]
    assert e["ratio"] == pytest.approx(0.547, abs=0.01)
    assert e["p"] < res["bonferroni_alpha"]
    assert e["pct_genomewide"] < 5.0
    assert e["ratio"] == min(v["ratio"] for v in res["per_gene"].values())


def test_rhoa_is_significant_but_not_specific(res):
    """The distinction the drift baseline exists to make.

    RhoA clears Bonferroni on a t test yet sits mid-distribution genome-wide, so
    its low p-value comes from low variance rather than a large effect. Reporting
    it as a specific finding would be wrong.
    """
    rho = res["per_gene"]["Rhoa"]
    assert rho["p"] < res["bonferroni_alpha"]
    assert rho["pct_genomewide"] > 20.0


def test_no_other_panel_gene_beats_the_drift(res):
    """Ect2 alone is below the 5th percentile of background drift."""
    outliers = [g for g, v in res["per_gene"].items() if v["pct_genomewide"] < 5.0]
    assert outliers == ["Ect2"]


def test_global_drift_is_downward_and_modest(res):
    """P7 cells have slightly lower complexity, so almost everything drifts down.
    Without this baseline a 0.9x change reads as meaningful."""
    d = res["drift"]
    assert 0.85 < d["median"] < 0.92
    assert d["q05"] < d["median"] < d["q95"]
    assert d["n_genes"] > 300


def test_module_difference_index_is_flat(res):
    """The heuristic's blind spot, stated as a test."""
    m = res["modules"]
    assert m["cytokinesis"]["ratio"] < 1.0 and m["mitotic"]["ratio"] < 1.0
    assert abs(m["difference_index"]["delta"]) < 0.05
    assert m["difference_index"]["p"] > 0.5


def test_maturation_orders_the_four_well_powered_groups(res):
    """P0-vCM < P7-aCM < P7-vCM, which reproduces a Baniol claim M was not built
    to test. The two P0-atrial groups (n=8, n=4) are excluded as underpowered."""
    g = res["maturation"]["by_group"]
    assert g["P0-vCM-cycling"]["mean"] < g["P7-aCM-cycling"]["mean"]
    assert g["P7-aCM-cycling"]["mean"] < g["P7-vCM-cycling"]["mean"]
    assert g["P0-vCM-noncycling"]["mean"] < g["P7-vCM-noncycling"]["mean"]


def test_maturation_is_graded_not_a_relabelled_stage(res):
    """Between-stage separation exceeds within-stage spread, but the within-stage
    SD is large enough that single cells can still be ordered."""
    w = res["maturation"]["within_stage_vCM"]
    sep = w["P7"]["mean"] - w["P0"]["mean"]
    assert sep > 1.0
    assert all(0.3 < w[s]["sd"] < 0.8 for s in ("P0", "P7"))
    assert sep > max(w[s]["sd"] for s in ("P0", "P7"))


def test_m_predicts_the_two_load_bearing_couplings(res):
    c = res["maturation"]["correlations"]
    assert c["Ect2"]["r"] < -0.45 and abs(c["Ect2"]["t"]) > 4
    assert c["E2f6"]["r"] > 0.30 and c["E2f6"]["t"] > 3
    assert c["Ect2"]["n"] == 89


def test_e2f_three_roles(res):
    e = res["e2f"]
    # cycling-restricted: E2f1/7/8 near-zero in noncycling
    for g in ("E2f1", "E2f7", "E2f8"):
        assert e[g]["P7-vCM-noncycling"]["mean"] < 0.06
        assert e[g]["P7-vCM-cycling"]["mean"] > 0.15
    # E2f6 is the only member present in noncycling cells, and it rises
    assert e["E2f6"]["P7-vCM-noncycling"]["mean"] > 0.30
    assert e["E2f6"]["P7-vCM-noncycling"]["detected_pct"] > 90
    assert e["E2f6"]["maturation_t"] > 3
    # E2f7 flat while E2f8 rises -- they must be separate model nodes
    assert abs(e["E2f7"]["maturation_t"]) < 1.0
    assert e["E2f8"]["maturation_t"] > 1.5
    # E2f2 is the endoreplication candidate: absent unless cycling, and rising
    assert e["E2f2"]["P7-vCM-noncycling"]["mean"] < 0.01
    assert e["E2f2"]["maturation_t"] > 1.5


def test_e2f1_induces_the_atypical_repressors(res):
    """Supports the E2F1 -> E2F7/8 delayed-feedback arm as a real edge."""
    for k in ("E2f1~E2f7", "E2f1~E2f8"):
        assert res["edges"][k]["rho"] > 0.4
    # and E2f6 is NOT part of the activator programme
    assert abs(res["edges"]["E2f1~E2f6"]["rho"]) < 0.25
    # E2f8 is G1/S-restricted: strongly anti-correlated with the S/G2 marker
    assert res["edges"]["E2f8~Ccna2"]["rho"] < -0.4


def test_fucci_reporter_module_cannot_be_validated_transcriptomically(res):
    """Gmnn and Cdt1 mRNA correlate POSITIVELY even though the FUCCI reporters
    they inspire are anti-correlated by phase -- the oscillation is
    post-translational. So the reporter layer must be checked against imaging."""
    assert res["edges"]["Gmnn~Cdt1"]["rho"] > 0


def test_circularity_audit_flags_ect2_and_e2f8(res):
    """Both are inside the phase-calling gene lists. Ect2's presence makes the
    maturation result conservative; E2f8's makes the KO's phase labels suspect."""
    c = res["circularity"]
    assert "Ect2" in c["in_G2M"]
    assert "E2f8" in c["in_S"]
    for g in ("Rhoa", "Racgap1", "Cep55", "Ccnb1", "Ccna2"):
        assert g in c["clean"], f"{g} should be a clean comparator"


def test_sort_enrichment_is_asymmetric(res):
    """P0 unenriched, P7 heavily enriched -- so the cycling fraction rises in the
    data and falls in reality."""
    se = res["sort_enrichment"]
    assert se["P0"]["enrichment"] == pytest.approx(1.0, abs=0.1)
    assert min(se["P7"]["enrichment_range"]) > 4.0
    assert se["P7"]["observed"] > se["P0"]["observed"]      # the trap, made explicit


def test_anchor_corrections(res):
    """Five node assignments the data overturns."""
    by = {a["node"]: a for a in res["anchors"]}
    # Myt1 rises with maturation; Wee1 falls -- the brake must be Pkmyt1
    assert by["mitotic brake"]["obvious"]["t"] < 0 < by["mitotic brake"]["use"]["t"]
    # p38gamma rises; p38alpha does not
    assert by["p38"]["obvious"]["t"] < 0 < by["p38"]["use"]["t"]
    # Cdh1 (E-cadherin) is near-absent; Fzr1 is the APC/C gene and is ubiquitous
    assert by["APC/C co-activator"]["obvious"]["det_P7"] < 10
    assert by["APC/C co-activator"]["use"]["det_P7"] > 90
    # Ccnd1 undetectable, Ccnd2 ubiquitous
    assert by["Cyclin D"]["obvious"]["det_P7"] < 15
    assert by["Cyclin D"]["use"]["det_P7"] > 90
    # alpha1A absent, alpha1B usable -- Murganti cited alpha1B and were right
    assert by["alpha1-adrenergic"]["obvious"]["det_P7"] < 10
    assert by["alpha1-adrenergic"]["use"]["det_P7"] > 40


def test_alpha2_receptors_are_unmeasurable_here(res):
    """Clonidine's canonical target is absent from the store, which bounds what
    the receptor layer can claim."""
    for g in ("Adra2a", "Adra2b", "Adra2c"):
        assert res["adrenergic"][g].get("absent") is True


def test_adrenergic_balance_tips_toward_the_brake(res):
    """Nischarin falls while alpha1B and beta1 rise -- a mechanism the papers did
    not state for why clonidine's effect should become less productive."""
    a = res["adrenergic"]
    assert a["Nisch"]["P7"] < a["Nisch"]["P0"]
    assert a["Adra1b"]["P7"] > a["Adra1b"]["P0"]
    assert a["Adrb1"]["P7"] > a["Adrb1"]["P0"]
    assert a["Nisch"]["det_P0"] > 95 and a["Nisch"]["det_P7"] > 95


def test_ddr_module_has_no_support_on_the_maturation_axis(res):
    """Baniol attributed the P7 G1/S delay to a DNA-damage response. On this axis
    all four genes are flat, so the DDR module is hypothesis-only."""
    assert all(abs(v["t"]) < 1.5 for v in res["ddr"].values())


def test_welch_matches_a_known_case():
    """Guard the hand-rolled t test against a worked example."""
    a = [1.0, 2.0, 3.0, 4.0, 5.0]
    b = [2.0, 3.0, 4.0, 5.0, 6.0]
    m1, m2, t, p = baniol.welch(a, b)
    assert m1 == 3.0 and m2 == 4.0
    assert t == pytest.approx(1.0, abs=1e-9)
    assert p == pytest.approx(0.3465935, abs=1e-4)


def test_pearson_and_spearman_edges():
    xs = [1.0, 2.0, 3.0, 4.0]
    assert baniol.pearson(xs, xs)[0] == pytest.approx(1.0)
    assert baniol.pearson(xs, [4.0, 3.0, 2.0, 1.0])[0] == pytest.approx(-1.0)
    # Spearman is rank-based, so a monotone nonlinear map is still +1
    assert baniol.spearman(xs, [1.0, 8.0, 27.0, 64.0])[0] == pytest.approx(1.0)
