"""Tests for the cmfate logic network: engine, spec invariants, and behaviour.

The behavioural tests are deliberately written as the claims they check, so a
regression reads as "the model stopped predicting X" rather than "a number moved".
"""
from __future__ import annotations

import math

import pytest

from cmcycle import logic, model, spec


# --------------------------------------------------------------------------- #
# engine
# --------------------------------------------------------------------------- #
def test_ec50_ceiling_is_enforced():
    """EC50**n >= 0.5 makes beta negative and inverts act() for every input.

    The reference implementation clamps only the upper bound, so the sign error
    propagates silently through `w - act(x) > w` and pushes activities above 1.
    """
    logic.Reaction(target=0, reactants=[(1, False)], w=1.0, n=1.4, ec50=0.55, rid="ok")
    with pytest.raises(ValueError, match="beta negative"):
        logic.Reaction(target=0, reactants=[(1, False)], w=1.0, n=1.4, ec50=0.70, rid="bad")
    assert logic.ec50_ceiling(1.4) == pytest.approx(0.6095, abs=1e-4)
    assert logic.ec50_ceiling(3.0) == pytest.approx(0.7937, abs=1e-4)


def test_reaction_rejects_degenerate_parameters():
    for kw in ({"n": 0.0, "ec50": 0.4}, {"n": 1.4, "ec50": 0.0}, {"n": 1.4, "ec50": 1.0}):
        with pytest.raises(ValueError):
            logic.Reaction(target=0, reactants=[], w=1.0, rid="x", **kw)


def test_activation_stays_in_bounds():
    r = logic.Reaction(target=0, reactants=[(1, False)], w=0.8, n=1.4, ec50=0.5)
    net = logic.Network(species=["a", "b"], yinit=[0, 0], ymax=[1, 1], tau=[1, 1],
                        reactions=[r], idx={"a": 0, "b": 1})
    for x in (-1.0, 0.0, 0.3, 1.0, 5.0):
        assert 0.0 <= net._act(x, r) <= r.w


def test_perturbation_matrix_accepts_a_tuple():
    """A double knockdown must be one column: the lab's experiment is a double KO,
    and the reference implementation could only do one node at a time."""
    net = spec.load_calibrated()
    mat, outs, cols = net.perturbation_matrix(
        perturb_nodes=["E2F7", ("E2F7", "E2F8")],
        output_nodes=["Division"], inputs=model.contexts(net)["mouse_p0_invivo"])
    assert cols == ["E2F7", "E2F7+E2F8"]
    assert len(mat) == 1 and len(mat[0]) == 2


def test_steady_state_raises_rather_than_returning_half_settled():
    """Node 'b' must be driven by 'a', or 'a' is a sink and converges trivially --
    sinks are excluded from the convergence test because their value is algebraic."""
    rin = logic.Reaction(target=0, reactants=[], w=1.0, n=1.4, ec50=0.5, is_input=True)
    rab = logic.Reaction(target=1, reactants=[(0, False)], w=1.0, n=1.4, ec50=0.5)
    net = logic.Network(species=["a", "b"], yinit=[0, 0], ymax=[1, 1], tau=[1, 1],
                        reactions=[rin, rab], idx={"a": 0, "b": 1})
    assert net.sink_indices() == {1}                 # only b is a pure output
    with pytest.raises(RuntimeError, match="steady state not reached"):
        net.steady_state(t_max=1e-3, tol=1e-30)


def test_sink_nodes_are_not_integrated():
    """The fate nodes carry tau up to 35 h; integrating them to convergence
    dominated the run time while changing no steady-state value."""
    net = spec.load_calibrated()
    sinks = {net.species[i] for i in net.sink_indices()}
    assert set(model.FATES) <= sinks
    assert "SPhase" not in sinks


# --------------------------------------------------------------------------- #
# spec invariants
# --------------------------------------------------------------------------- #
def test_spec_lints_clean():
    assert spec.lint(spec.load()) == []


def test_loader_is_strict_about_typos():
    with pytest.raises(ValueError, match="no '=>'"):
        spec.parse_rule("A & B")
    with pytest.raises(ValueError, match="empty reactant"):
        spec.parse_rule("A & => B")
    assert spec.parse_rule("=> A") == ("A", [], True)
    assert spec.parse_rule("C & !D => E") == ("E", [("C", False), ("D", True)], False)


def test_fate_partition_sums_to_one_everywhere():
    """The complementary product makes Q+D+B+P = 1 algebraically, for any input."""
    net = spec.load_calibrated()
    for key in net.meta["contexts"]:
        for clon in (0.0, 0.5, 1.0):
            f = model.state(net, key, clonidine=clon)
            total = sum(f[k] for k in model.FATES)
            assert total == pytest.approx(1.0, abs=1e-6), f"{key} clon={clon}"


def test_calibration_hits_murganti_exactly():
    """Three measured fractions, inverted, give three independent 1-D solves."""
    net = spec.load_calibrated()
    f = model.state(net, net.meta["calibration"]["context"])
    for k, v in (("Quiescent", 0.9035), ("Division", 0.0509),
                 ("Binucleation", 0.0316), ("Polyploidization", 0.0140)):
        assert f[k] == pytest.approx(v, abs=5e-4)


def test_calibration_cache_is_keyed_on_the_spec():
    """A stale fit silently applied to a changed network is the failure to design
    out, so the cache carries a hash of the spec files."""
    import json
    from importlib.resources import files
    blob = json.loads(files("cmcycle").joinpath("data", spec.CACHE).read_text())
    assert blob["fingerprint"] == spec._spec_fingerprint()


# --------------------------------------------------------------------------- #
# behaviour
# --------------------------------------------------------------------------- #
def test_the_clonidine_triad_is_reproduced():
    """One drug, three contexts, one input change. This is the held-out test."""
    rows = model.clonidine_triad(spec.load_calibrated())
    misses = [r["context"] for r in rows if not r["hit"]]
    assert not misses, f"wrong dominant fate in {misses}"


def test_entry_response_is_over_predicted_at_high_maturation():
    """A known, reported limitation rather than a passing claim.

    Clonidine's leverage on entry scales with baseline PKA, which the model ties to
    beta-adrenergic tone; that tone rises steeply with maturation, so the entry
    fold is far too large in the mature contexts. Pinned here so a fix is visible
    as a test change rather than a silent improvement.
    """
    rows = {r["context"]: r for r in model.clonidine_triad(spec.load_calibrated())}
    assert 1.4 < rows["hipsc_cm"]["entry_fold"] < 3.0          # observed ~2.4x
    assert rows["mncm_invitro"]["entry_fold"] > 4.0            # observed ~2.1x -- too high
    assert rows["mouse_p1_invivo"]["entry_fold"] > 4.0         # observed ~1.5x -- too high


def test_productive_share_falls_monotonically_with_maturation():
    """The headline structure: entry is necessary but not sufficient."""
    net = spec.load_calibrated()
    rows = model.maturation_series(net, values=[0.1, 0.3, 0.5, 0.7, 0.9])
    shares = [r["productive_share"] for r in rows]
    assert all(a >= b - 1e-9 for a, b in zip(shares, shares[1:])), shares
    # This sweep holds load and beta-AR tone at their in-vivo values, so the
    # absolute share is low throughout; the monotonic fall is the claim.
    assert shares[0] > 10 * shares[-1]
    ect2 = [r["Ect2"] for r in rows]
    assert all(a >= b - 1e-9 for a, b in zip(ect2, ect2[1:])), ect2
    # The hiPSC CONTEXT differs in more than M and does reach a high share -- which
    # is the two-factor point: maturation alone does not set the outcome.
    h = model.state(net, "hipsc_cm")
    cyc = sum(h[f] for f in model.CYCLING_FATES)
    assert h["Division"] / cyc > 0.4


def test_ect2_is_rate_limiting_for_division():
    """Knocking Ect2 down must cost division far more than it costs entry -- the
    model's central claim, and the reason RhoA is obligatory in the furrow rule."""
    net = spec.load_calibrated()
    inp = model.contexts(net)["mouse_p0_invivo"]
    base = net.ss(inputs=inp)
    ko = net.ss(inputs=inp, knockdowns=["Ect2"])
    assert ko["Division"] < 0.5 * base["Division"]
    assert ko["SPhase"] > 0.9 * base["SPhase"]


def test_negative_controls_do_not_move_with_maturation():
    """Free controls from the data: E2f3, AurKB, Ccna2 and the comparator arm are
    flat on the maturation axis, so the model must not swing them."""
    net = spec.load_calibrated()
    lo = model.state(net, "mouse_p0_invivo")
    hi = model.state(net, "mouse_p7_invivo")
    for n in ("E2F3", "AurKB", "CycA", "Centralspindlin", "Anillin"):
        rel = abs(hi[n] - lo[n]) / max(lo[n], 1e-9)
        assert rel < 0.25, f"{n} moved {rel:.0%} between P0 and P7"


def test_e2f78_knockdown_raises_entry_and_ect2():
    """Removing the atypical repressors de-represses E2F activity, and E2F7/8 also
    repress Ect2 (Pandit 2012) -- so both entry and cytokinesis competence rise."""
    for row in model.e2f78_knockdown(spec.load_calibrated()):
        d = row["double"]
        assert d["E2Fact"] > 0
        assert d["SPhase"] > 0
        assert d["Ect2"] > 0


def test_the_double_knockdown_is_not_the_sum_of_the_singles():
    """Epistasis is a prediction, and the lab's data is a double knockout."""
    rows = model.e2f78_knockdown(spec.load_calibrated())
    assert any(abs(r["epistasis"]["E2Fact"]) > 1e-4 for r in rows)


def test_a_hit_needs_both_a_fold_and_an_absolute_effect():
    """The relative test alone let a 0.4-percentage-point change at P7 score as a
    '1.5x hit'. Nine of eleven hits were that, which is what made the screen look
    twice as permissive as the bench."""
    net = spec.load_calibrated()
    rows, _, base = model._sweep(net, "mouse_p7_invivo")
    assert base["entry"] < 0.02, "P7 baseline entry should be tiny -- that is the trap"
    for r in rows:
        if r["entry_hit"]:
            assert r["d_entry"] >= model.HIT_FLOOR
            assert r["d_entry"] + base["entry"] >= model.HIT_FOLD * base["entry"]


def test_specificity_is_compared_like_for_like():
    """Murganti's screen was in hiPSC-CM, so ours has to be too. Comparing a P7
    screen against their number was never valid -- baseline entry differs 36-fold."""
    net = spec.load_calibrated()
    sp = model.screen_specificity(net)
    assert sp["context"] == "hipsc_cm"
    assert sp["baseline_entry"] > 0.05
    p7 = model.screen_specificity(net, context="mouse_p7_invivo")
    assert sp["baseline_entry"] / p7["baseline_entry"] > 10


def test_an_entry_scored_screen_mostly_finds_things_that_do_not_help():
    """The claim the screen exists to make. Murganti scored percent-mVenus-positive,
    i.e. S-phase entry, and in this partition every fate scales with entry -- so most
    entry hits do not raise the division share. Their own 2 of 6 surviving validation
    is the wet version of this number."""
    sp = model.screen_specificity(spec.load_calibrated())
    assert sp["entry_hits"] > 0
    assert sp["entry_hits_that_help"] < 0.5
    assert abs(sp["entry_hit_rate"] - sp["wet_hit_rate"]) < 0.06
    assert abs(sp["entry_hits_that_help"] - sp["wet_survival"]) < 0.20
