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


def test_negative_controls_move_only_as_much_as_e2f_activity_does():
    """The five comparators are flat in the data (|t| <= 1.73 for E2f3, AurKB, Ccna2,
    Racgap1/Kif23, Anln), so the model must not swing them *selectively*.

    Fixing the G1/S switch cost 28% here, and the bound was raised from 25% to
    accommodate it -- which needs justifying rather than waving through:

    * the drift is **uniform**: all five move 27-28%, which is the signature of one
      shared upstream cause rather than five spurious edges. It tracks E2Fact falling
      0.358 -> 0.255 (29%) between P0 and P7 -- i.e. it is the switch now working;
    * the measurements are **conditional** -- P0-cycling versus P7-cycling cells --
      while these model nodes are unconditional per-cell activities, so the comparison
      was never exact. This is the same conditional-versus-prevalence distinction that
      required CycA to hang off E2Fact rather than SPhase;
    * the test still catches what it is for: any comparator moving *differently* from
      the others, which is what a wrongly-wired edge would look like.

    The residual tension is real and recorded in TODO: truly flat targets under a
    travelling E2F activity would need their reactions to saturate, so that a 29% fall
    in input gives a small fall in output. That is a concrete refinement, not a fudge.
    """
    net = spec.load_calibrated()
    lo = model.state(net, "mouse_p0_invivo")
    hi = model.state(net, "mouse_p7_invivo")
    moves = {n: abs(hi[n] - lo[n]) / max(lo[n], 1e-9)
             for n in ("E2F3", "AurKB", "CycA", "Centralspindlin", "Anillin")}
    for n, rel in moves.items():
        assert rel < 0.35, f"{n} moved {rel:.0%} between P0 and P7"
    # uniformity is the actual assertion: a selective swing means a mis-wired edge
    assert max(moves.values()) - min(moves.values()) < 0.05, moves
    # and it must be no larger than the E2F activity change driving it
    e2f_drop = abs(hi["E2Fact"] - lo["E2Fact"]) / lo["E2Fact"]
    assert max(moves.values()) <= e2f_drop + 0.02, (moves, e2f_drop)


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


def test_compiled_cache_notices_a_structural_edit():
    """Keying the compiled-reaction cache on weights alone was a correctness bug:
    replacing a reaction to change its reactants left the key unchanged, so the edit
    was silently ignored. It produced a diagnosis that had to be thrown away."""
    net = spec.load_calibrated()
    inp = model.contexts(net)["hipsc_cm"]
    before = net.ss(inputs=inp)["Rb"]
    j = net.meta["rid"]["r044"]
    r = net.reactions[j]
    # drop the CycE leg of `!CycD & !CycE => Rb`, changing nothing else
    net.reactions[j] = logic.Reaction(
        target=r.target, reactants=[(net.idx["CycD"], True)],
        w=r.w, n=r.n, ec50=r.ec50, is_input=r.is_input, rid="r044")
    after = net.ss(inputs=inp)["Rb"]
    assert after != pytest.approx(before, abs=1e-6), (
        "structural edit was ignored -- the compiled cache did not invalidate")


def test_the_g1s_switch_actually_travels():
    """Was a pinned DEFECT; now the fix. The restriction point used to contribute
    almost nothing to entry -- E2F1 spanned 1.4x and CycE 1.2x across every context
    while entry spanned 142x, so all of the variation came from the CKI brakes inside
    one reaction, downstream of the switch.

    Three changes fixed it: the three restriction-point reactions were gate-shaped
    (they were the only graded ones in a model whose two downstream gates are steep),
    an OR'd CycE leg that cancelled the switch's travel was removed, and CycD's drive
    was scaled to put the fold inside the context range.
    """
    net = spec.load_calibrated()
    span = lambda n: (lambda v: max(v) / max(min(v), 1e-12))(
        [model.state(net, c)[n] for c in net.meta["contexts"]])
    assert span("E2F1") > 10.0, "the switch has lost its travel again"
    assert span("CycE") > 100.0
    # and Rb now runs from open at hiPSC to closed in the adult
    assert model.state(net, "hipsc_cm")["Rb"] < 0.25
    assert model.state(net, "adult")["Rb"] > 0.8


def test_entry_is_now_driven_through_the_switch():
    """The consequence of the fix, and the reason it mattered. Before, no upstream
    route reached entry -- ERK and Autophagy 1.00x, CycD 1.10x -- which is why
    clonidine's effect had to be wired directly onto the S-phase reaction as `!PKA`.
    That shortcut is gone, so upstream drive must now get through."""
    net = spec.load_calibrated()
    rids = {net.reactions[net.meta["rid"]["r063"]].rid}
    assert rids == {"r063"}
    sphase = net.reactions[net.meta["rid"]["r063"]]
    names = {net.species[i] for i, _ in sphase.reactants}
    assert "PKA" not in names, "the !PKA shortcut is back on the S-phase reaction"
    inp = model.contexts(net)["mouse_p7_invivo"]
    base = net.ss(inputs=inp)["SPhase"]
    got = net.ss(inputs=inp, overexpress=["CycD"])["SPhase"]
    assert got / base > 1.5, "upstream drive still cannot reach entry"


def test_the_clonidine_entry_folds_are_close_at_the_mature_contexts():
    """What the fix bought quantitatively: mean fold error across the three contexts
    fell from 236% to 26%, with the two mature contexts nearly exact. The residual is
    hiPSC-CM, where the switch is already open so relieving a brake cannot open it
    further -- clonidine's real 2.44x there must come through something else."""
    rows = {r["context"]: r for r in model.clonidine_triad(spec.load_calibrated())}
    assert rows["mncm_invitro"]["entry_fold"] == pytest.approx(2.12, rel=0.15)
    assert rows["mouse_p1_invivo"]["entry_fold"] == pytest.approx(1.52, rel=0.30)
    assert rows["hipsc_cm"]["entry_fold"] < 1.3          # the remaining miss


def test_the_mouse_in_vivo_cycling_level_is_reproduced_unfitted():
    """A genuine held-out validation, and a lesson about comparing observables.

    Baniol's Fig 1D gives FUCCI *states*, and the naive cycling fraction (1 - mKO2+)
    is 32.5% at P0. The model cannot reach that at any weight, which looked like a
    2x miss. But their own Suppl 1G shows only 22.7% of the G1/S double-positives are
    Ki67+ at P0 -- the rest have prematurely exited, which is a distinction they make
    explicitly. Correcting with their own co-staining:

        mAG+ 12.36% (all Ki67+) + 22.7% of 19.1% + mitotic 1%  =  17.7%

    against a model ceiling of 17.2%. Nothing was fitted to this.
    """
    net = spec.load_calibrated()
    corrected_p0 = 0.1236 + 0.191 * 0.227 + 0.010
    r = net.reactions[net.meta["rid"]["r063"]]
    saved, r.w = r.w, 1.0
    try:
        ceiling = 1 - net.ss(inputs=model.contexts(net)["mouse_p0_invivo"])["Quiescent"]
    finally:
        r.w = saved
    assert ceiling == pytest.approx(corrected_p0, rel=0.15)
    # and the naive figure is out of reach, which is the point of the correction
    assert ceiling < 0.325


def test_the_maturation_slope_of_entry_is_too_steep():
    """A pinned residual. Corrected for Ki-67 as above, Baniol's P0-to-P7 fall in
    cycling cardiomyocytes is 17.7% -> 5.5%, a factor of 3.2. The model falls ~13x.

    Same root cause as the CycE over-swing: SPhase multiplies four maturation-dependent
    factors (CycE and three brakes), and after the G1/S fix CycE travels ~36,000x, which
    compounds into far too steep a slope. Flattening CycE would address both.
    """
    net = spec.load_calibrated()
    r = net.reactions[net.meta["rid"]["r063"]]
    saved = r.w
    tops = {}
    try:
        r.w = 1.0
        for c in ("mouse_p0_invivo", "mouse_p7_invivo"):
            tops[c] = 1 - net.ss(inputs=model.contexts(net)[c])["Quiescent"]
    finally:
        r.w = saved
    slope = tops["mouse_p0_invivo"] / tops["mouse_p7_invivo"]
    observed = (0.1236 + 0.191 * 0.227 + 0.010) / (0.0083 + 0.079 * 0.588)
    assert observed == pytest.approx(3.2, abs=0.3)
    assert slope > 2 * observed, "slope has been flattened -- update this test"
