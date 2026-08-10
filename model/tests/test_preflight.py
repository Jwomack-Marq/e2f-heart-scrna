"""Pre-flight consistency checks over the cmcycle calibration set.

These are stdlib-only and must stay that way -- the whole point of the pre-flight
is that it runs before any solver, in any interpreter, so that an inconsistent
dataset combination is caught before 30-odd differential equations are written.
"""
from __future__ import annotations

import math

import pytest

from cmcycle import preflight as pf


def test_targets_load_and_are_cited():
    targets = pf.load_targets()
    assert len(targets) > 30
    assert all(t.cite for t in targets), "every target needs a source"
    assert all(t.context for t in targets)


def test_error_kind_is_explicit_and_round_trips():
    """SD/SEM must never be inferred. Conflating them rescales by sqrt(n)."""
    for t in pf.load_targets():
        if t.err is None:
            continue
        assert t.err_kind in ("sd", "sem")
        if t.n:
            assert t.sd is not None and t.sem is not None
            assert t.sd == pytest.approx(t.sem * math.sqrt(t.n))
        if t.err_kind == "sd":
            assert t.sd == t.err
        else:
            assert t.sem == t.err


def test_percentage_sd_cannot_exceed_100():
    """The n-vs-n_obs guard. Murganti's fate fractions are 570 cells from 20
    movies with SEM across movies; using 570 for the conversion yields 253%."""
    for t in pf.load_targets():
        if t.unit == "pct" and t.sd is not None:
            assert t.sd <= 100, f"{t.name}/{t.context} derived SD {t.sd:.0f}%"


def test_fate_fractions_sum_to_100():
    t = pf.load_targets()
    total = sum(pf.get(t, f"fate_{k}", "hiPSC_CM").value for k in
                ("noncycling", "division", "binucleation", "polyploidization"))
    assert total == pytest.approx(100.0, abs=0.05)


def test_hipsc_reservoir_is_an_excess_not_a_shortfall():
    """The imaging-derived event rate under-predicts the flow mAG+ fraction.

    Direction matters: an excess implies cells that stay geminin-positive without
    resolving into a scored outcome, i.e. sustained G2 arrest. A shortfall would
    instead mean the flow gate was too strict.
    """
    c = pf.check_hipsc_reservoir(pf.load_targets())
    assert c.observed > c.predicted
    assert c.ratio == pytest.approx(2.32, abs=0.05)


def test_mouse_edu_agrees_within_a_third():
    """Baniol's FUCCI snapshots predict Murganti's EdU index, nothing fitted."""
    c = pf.check_mouse_edu(pf.load_targets())
    assert 0.7 < c.ratio < 1.4, f"got {c.ratio:.2f}"


def test_edu_fold_exceeds_ki67_fold_which_saturation_forbids():
    """A pure entry-rate rise always gives fold(EdU) < fold(Ki-67).

    EdU is 1-exp(-rT) (sublinear in r); Ki-67 is ~r*D_K (linear). The mNCM data
    has the opposite ordering, which is positive evidence that Ki-67 is lost from
    arrested cells while the EdU label persists.
    """
    c = pf.check_ki67_edu_sign(pf.load_targets())
    assert c.observed > c.predicted

    # the forbidden direction, demonstrated
    for fold_r in (1.5, 2.0, 3.0):
        base = 0.10
        edu_fold = (1 - math.exp(-base * fold_r)) / (1 - math.exp(-base))
        assert edu_fold < fold_r


def test_saturation_correction_increases_the_inferred_rate_fold():
    naive = 36.15 / 25.54
    corrected = pf.entry_rate_fold(25.54, 36.15)
    assert corrected > naive
    assert corrected == pytest.approx(1.52, abs=0.01)


def test_fate_budget_identity_holds_on_the_mouse_arm():
    """dEdU ~ dBinuc + dPoly + 2*dDiv, converting per-cell to per-nucleus."""
    c = pf.check_fate_budget(pf.load_targets())
    assert c.verdict == "CONSISTENT"
    assert 0.6 < c.predicted / c.observed < 1.0


def test_mitotic_dwell_time_is_several_hours():
    """~4 h, not the ~1 h of cycling somatic cells."""
    c = pf.check_mitotic_duration(pf.load_targets())
    assert 3.0 < c.predicted < 5.5


def test_run_all_flags_exactly_the_two_known_tensions():
    checks = pf.run_all()
    flagged = {c.name for c in checks if c.verdict.startswith("INCONSISTENT")}
    assert flagged == {"hiPSC arrested reservoir", "Ki-67 vs EdU fold sign"}
