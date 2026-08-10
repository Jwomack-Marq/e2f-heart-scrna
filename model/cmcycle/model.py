"""The calibrated ``cmfate`` model, and the experiments run on it.

Everything here takes the network from :func:`cmcycle.spec.load_calibrated`, which
is fitted at one context (hiPSC-CM) and then asked about others. Each function
returns plain dicts so the figures and the report cannot disagree with the model.

What is fitted, and what is not
-------------------------------
FITTED (4 quantities, all at the hiPSC context):
  * three gate weights, solved by bisection against Murganti's fate fractions
    inverted through the fate partition;
  * the abscission switch position (r080's n and EC50), chosen so the switch sits
    inside the Ect2/RhoA arm's operating range.

ADJUSTED during development, so not independent evidence:
  * the E2F repressor-pool weights (an OR pool of three saturates);
  * the abscission arm's parameter class (attenuating -> amplifying);
  * removal of two structures that defeated the model's own claim -- an OR'd floor
    on Ect2, and an OR'd Midbody route that bypassed RhoA entirely.

PREDICTED (not tuned against):
  * which fate dominates in each of the three clonidine contexts;
  * the direction of the E2f7/E2f8 double-knockdown effect and its
    maturation-dependence;
  * the flatness of the declared negative controls.
"""
from __future__ import annotations

from . import spec

FATES = ("Quiescent", "Division", "Binucleation", "Polyploidization")
CYCLING_FATES = ("Division", "Binucleation", "Polyploidization")

#: What the papers observed, per context. Held out of the calibration.
OBSERVED_TRIAD = {
    "hipsc_cm": {"fate": "Division", "ki67_fold": 22.0 / 9.0,
                 "note": "AurKB+ midbodies doubled; ploidy and binucleation flat"},
    "mncm_invitro": {"fate": "Polyploidization", "ki67_fold": 23.23 / 10.97,
                     "note": "tetraploid 26.5->38.9%; ZERO midbodies in >10,000 cells"},
    "mouse_p1_invivo": {"fate": "Binucleation", "ki67_fold": 1.52,
                        "note": "mononucleated 22.3->14.8%; ploidy flat; EdU fold is "
                                "saturation-corrected"},
}


def load(**kw):
    """The calibrated network."""
    return spec.load_calibrated(**kw)


def contexts(net) -> dict:
    return {k: v["inputs"] for k, v in net.meta["contexts"].items()}


def state(net, context, clonidine=0.0, **kw) -> dict:
    """Steady state in a named context, optionally dosed."""
    inp = dict(net.meta["contexts"][context]["inputs"])
    if clonidine:
        inp["Clonidine"] = float(clonidine)
    return net.ss(inputs=inp, **kw)


def maturation_series(net, values=None, clonidine=0.0) -> list:
    """Sweep the maturation coordinate at fixed in-vivo context.

    This is the model's headline structure: the fate the cycling flux lands in
    changes with maturation while entry falls monotonically.
    """
    base = dict(net.meta["contexts"]["mouse_p1_invivo"]["inputs"])
    values = values or [round(0.05 * i, 3) for i in range(1, 20)]
    out = []
    for m in values:
        inp = {**base, "Maturation": m}
        if clonidine:
            inp["Clonidine"] = float(clonidine)
        s = net.ss(inputs=inp)
        row = {"Maturation": m, **{f: s[f] for f in FATES},
               **{g: s[g] for g in net.meta["gate_nodes"]},
               "Ect2": s["Ect2"], "Ki67": s["Ki67"]}
        cyc = sum(row[f] for f in CYCLING_FATES)
        row["productive_share"] = row["Division"] / cyc if cyc else 0.0
        out.append(row)
    return out


def clonidine_triad(net) -> list:
    """One drug, three contexts, one input change. The held-out test.

    Returns a row per context with the predicted dominant fate, the observed one,
    and the entry-response fold beside the measured value.
    """
    rows = []
    for key, obs in OBSERVED_TRIAD.items():
        inp = dict(net.meta["contexts"][key]["inputs"])
        s0 = net.ss(inputs=inp)
        s1 = net.ss(inputs={**inp, "Clonidine": 1.0})
        delta = {f: s1[f] - s0[f] for f in CYCLING_FATES}
        pred = max(delta, key=delta.get)
        rows.append({
            "context": key,
            "label": net.meta["contexts"][key]["label"],
            "Maturation": inp["Maturation"], "InVitro": inp.get("InVitro", 0.0),
            "gates": {g: s0[g] for g in net.meta["gate_nodes"]},
            "base": {f: s0[f] for f in FATES},
            "dosed": {f: s1[f] for f in FATES},
            "delta": delta,
            "predicted": pred, "observed": obs["fate"], "hit": pred == obs["fate"],
            "entry_fold": (s1["SPhase"] / s0["SPhase"]) if s0["SPhase"] else float("nan"),
            "ki67_fold": (s1["Ki67"] / s0["Ki67"]) if s0["Ki67"] else float("nan"),
            "ki67_fold_observed": obs["ki67_fold"],
            "note": obs["note"],
        })
    return rows


def e2f78_knockdown(net, contexts_=("mouse_p0_invivo", "mouse_p7_invivo")) -> list:
    """In-silico E2f7/E2f8 double knockdown -- the lab's own experiment.

    Reports the singles and the double so the epistasis is visible: the difference
    between the double and the sum of the singles is a model prediction, and the
    lab's data is a double knockout.
    """
    rows = []
    for key in contexts_:
        inp = dict(net.meta["contexts"][key]["inputs"])
        mat, _, cols = net.perturbation_matrix(
            perturb_nodes=["E2F7", "E2F8", ("E2F7", "E2F8")],
            output_nodes=list(FATES) + ["SPhase", "Ect2", "E2Fact", "Ki67"],
            inputs=inp, mode="knockdown")
        keys = list(FATES) + ["SPhase", "Ect2", "E2Fact", "Ki67"]
        by = {c: {k: mat[i][j] for i, k in enumerate(keys)} for j, c in enumerate(cols)}
        dbl, s7, s8 = by["E2F7+E2F8"], by["E2F7"], by["E2F8"]
        rows.append({
            "context": key, "Maturation": inp["Maturation"],
            "single_E2F7": s7, "single_E2F8": s8, "double": dbl,
            "epistasis": {k: dbl[k] - (s7[k] + s8[k]) for k in keys},
        })
    return rows


#: A perturbation counts as a hit only if it clears BOTH a relative fold and an
#: absolute floor. The relative test alone is what inflated the first version of
#: these numbers: at P7 the baseline entry is ~1%, so a 0.4-percentage-point change
#: scored as a "1.5x hit". Nine of eleven hits were that.
#:
#: The absolute floor is a stated convention, not a measurement -- Murganti do not
#: publish a numeric hit criterion. 1 percentage point is chosen to sit comfortably
#: above their screen's own noise: CV 12.1% at its 48 h optimum, on a baseline
#: mVenus+ fraction near 10%, is roughly 1.2 pp of scatter.
HIT_FOLD = 1.5
HIT_FLOOR = 0.01


def _sweep(net, context, fold=HIT_FOLD, floor=HIT_FLOOR):
    """One pass over every single-node knockdown and overexpression.

    Both the ranking and the specificity check need the same sweep, and it is the
    slowest thing here, so it runs once.
    """
    inp = dict(net.meta["contexts"][context]["inputs"])
    # The gate nodes and their raw pre-gate inputs are excluded: "overexpress
    # Abscission" trivially maximises the division share because Abscission IS the
    # third factor of that share. Ranking them would be circular, not a finding.
    circular = set(net.meta["gate_nodes"]) | {"MitCompRaw", "AbsRaw"}
    skip = set(net.input_nodes()) | set(FATES) | circular
    nodes = [s for s in net.species if s not in skip]
    base = net.ss(inputs=inp)

    def share(st):
        cyc = sum(st[f] for f in CYCLING_FATES)
        return (st["Division"] / cyc) if cyc else 0.0

    b_share, b_entry = share(base), 1 - base["Quiescent"]
    rows = []
    for mode, kwname in (("knockdown", "knockdowns"), ("overexpress", "overexpress")):
        for n in nodes:
            try:
                st = net.ss(inputs=inp, **{kwname: [n]})
            except RuntimeError:
                continue
            entry, sh = 1 - st["Quiescent"], share(st)
            rows.append({
                "node": n, "mode": mode,
                "d_share": sh - b_share,
                "d_entry": entry - b_entry,
                "d_division": st["Division"] - base["Division"],
                # both criteria, so a tiny absolute change against a tiny baseline
                # cannot masquerade as a hit
                "entry_hit": entry >= fold * b_entry and (entry - b_entry) >= floor,
                "share_hit": (sh - b_share) >= floor,
            })
    return rows, len(nodes) * 2, {"entry": b_entry, "share": b_share}


def regeneration_screen(net, context="mouse_p7_invivo", top=12) -> list:
    """Rank every node by its ability to convert non-productive cycling into division.

    Scored on the CONDITIONAL division share, D/(D+B+P), not on Delta-Division and
    not on entry. That distinction is the point: Murganti's wet screen scored
    percent-mVenus-positive, i.e. S-phase entry, and in this partition every fate
    scales with entry -- so a screen scored that way is confounded with respect to
    productive division.

    Ranked at P7 by default, because "what would help a mature cardiomyocyte" is the
    useful question. Note that :func:`screen_specificity` deliberately runs somewhere
    else -- see its docstring.
    """
    rows, _, _ = _sweep(net, context)
    rows.sort(key=lambda r: -r["d_share"])
    return rows[:top]


def screen_specificity(net, context="hipsc_cm", fold=HIT_FOLD,
                       floor=HIT_FLOOR) -> dict:
    """Does an entry-scored screen find the same things as a division-scored one?

    Run in **hiPSC-CM by default, because that is where Murganti ran theirs**. An
    earlier version compared a P7 screen against their hiPSC-CM numbers, which was
    never a like-for-like comparison: baseline entry differs by an order of
    magnitude between the two, and the hit criterion is partly relative.

    Reports the two hit sets and their disagreement rather than a single
    "conversion" rate. The previous conversion figure counted ``d_share > 0``, and
    almost every value was epsilon-positive, so it measured nothing.

    Wet comparator, like-for-like: Murganti found 6/94 = 6.4% hits in hiPSC-CM, of
    which 2/6 = 33% still worked in the more mature mNCM.
    """
    from .baniol import spearman

    rows, n_tested, base = _sweep(net, context, fold=fold, floor=floor)
    entry_hits = [r for r in rows if r["entry_hit"]]
    share_hits = [r for r in rows if r["share_hit"]]
    both = [r for r in entry_hits if r["share_hit"]]

    # rank discordance: if an entry-scored screen were a good proxy for a
    # division-scored one, these two orderings would agree
    rho, _ = spearman([r["d_entry"] for r in rows], [r["d_share"] for r in rows])

    return {
        "context": context,
        "n_tested": n_tested,
        "baseline_entry": base["entry"],
        "baseline_share": base["share"],
        "criterion": {"fold": fold, "floor_pp": floor * 100},
        "entry_hits": len(entry_hits),
        "entry_hit_rate": len(entry_hits) / n_tested if n_tested else 0.0,
        "share_hits": len(share_hits),
        "both": len(both),
        "entry_hits_that_help": (len(both) / len(entry_hits)) if entry_hits else 0.0,
        "rank_agreement_rho": rho,
        "wet_hit_rate": 6 / 94,
        "wet_survival": 2 / 6,
    }


def summary(net, screen_context="mouse_p7_invivo") -> dict:
    """Everything the report quotes, in one call, with the node sweep run once."""
    rows, n_tested, _ = _sweep(net, screen_context)
    ranked = sorted(rows, key=lambda r: -r["d_share"])
    return {
        "n_species": len(net.species),
        "n_reactions": len(net.reactions),
        "inputs": net.input_nodes(),
        "calibrated": net.meta["calibrated"],
        "fit_context": net.meta["calibration"]["context"],
        "gate_targets": net.meta["calibration"]["targets"],
        "triad": clonidine_triad(net),
        "maturation": maturation_series(net),
        "ko": e2f78_knockdown(net),
        "screen": ranked[:12],
        "screen_worst": ranked[-6:],
        "screen_context": screen_context,
        # deliberately a different context from the ranking -- see the docstring
        "specificity": screen_specificity(net),
    }
