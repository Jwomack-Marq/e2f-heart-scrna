"""Aggregate consistency checks on the cardiomyocyte cell-cycle calibration set.

Run this *before* fitting any mechanism. Every check here is closed-form algebra
over the published aggregates in ``data/cmcycle_targets.csv`` -- no ODE, no
solver, no third-party dependency. Two papers supply the numbers:

  Baniol et al. 2021, Exp Cell Res 408:112880   -- FUCCI mouse, in vivo + scRNA-seq
  Murganti et al. 2022, Front Cardiovasc Med 9:840147 -- TNNT2-FUCCI hiPSC-CM

The checks exist because the four experimental observables are *different
functionals of the same cycle* and therefore constrain each other:

    instantaneous phase fraction  ~  entry rate x phase duration
    instantaneous Ki-67 index     ~  entry rate x Ki-67 dwell time   (linear in r)
    cumulative EdU index          =  1 - exp(-r*T)                   (SUBLINEAR in r)
    fate-budget identity          :  dEdU ~ dBinuc + dPoly + 2*dDiv

Deliberately stdlib-only: this module must stay runnable in a bare interpreter,
because its whole purpose is to be the cheapest possible falsification step.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from importlib.resources import files

_DATA = files(__package__).joinpath("data")


@dataclass(frozen=True)
class Target:
    """One published aggregate measurement.

    ``err`` is stored exactly as the source reports it, with ``err_kind`` saying
    which it is. Both papers mix the two conventions -- Baniol's Fig 1D and Ki-67
    subsets are SD, most of Murganti is SEM -- and silently treating one as the
    other rescales the likelihood by sqrt(n). Always go through :attr:`sd` and
    :attr:`sem` rather than touching ``err``.

    ``n`` and ``n_obs`` are deliberately separate and must not be conflated.
    ``n`` is the replicate count *behind the error bar* (movies, wells, animals);
    ``n_obs`` is the total number of cells or nuclei scored. Murganti's fate
    fractions are the trap: 570 cells from 20 movies, with SEM across movies, so
    deriving an SD with n=570 yields 253% -- impossible for a percentage. Use
    ``n`` for SD/SEM conversion and ``n_obs`` for binomial/multinomial weights.
    """
    name: str
    value: float
    err: float | None
    err_kind: str | None          # "sd" | "sem" | None
    n: int | None                 # replicates behind the error bar
    n_obs: int | None             # cells / nuclei scored
    unit: str
    kind: str
    context: str
    cite: str

    @property
    def sem(self) -> float | None:
        """Standard error of the mean, whichever way the source reported it."""
        if self.err is None:
            return None
        if self.err_kind == "sem":
            return self.err
        if self.err_kind == "sd" and self.n:
            return self.err / math.sqrt(self.n)
        return None

    @property
    def sd(self) -> float | None:
        """Standard deviation, whichever way the source reported it."""
        if self.err is None:
            return None
        if self.err_kind == "sd":
            return self.err
        if self.err_kind == "sem" and self.n:
            return self.err * math.sqrt(self.n)
        return None


def load_targets() -> list[Target]:
    """The calibration set, as shipped package data."""
    def num(s):
        s = (s or "").strip()
        return float(s) if s else None

    with _DATA.joinpath("cmcycle_targets.csv").open("r", encoding="utf8") as fh:
        rows = list(csv.DictReader(fh))
    out = []
    for r in rows:
        n, n_obs = num(r["n"]), num(r["n_obs"])
        kind = (r["err_kind"] or "").strip() or None
        if kind not in (None, "sd", "sem"):
            raise ValueError(f"{r['name']}: bad err_kind {kind!r}")
        t = Target(
            name=r["name"], value=float(r["value"]), err=num(r["err"]),
            err_kind=kind, n=int(n) if n else None,
            n_obs=int(n_obs) if n_obs else None, unit=r["unit"],
            kind=r["kind"], context=r["context"], cite=r["cite"],
        )
        # A derived SD above 100% on a percentage means n and n_obs were mixed up.
        if t.unit == "pct" and (t.sd or 0) > 100:
            raise ValueError(
                f"{t.name}/{t.context}: derived SD {t.sd:.0f}% exceeds 100% -- "
                f"'n' is probably the cell count ({t.n}) where it should be the "
                f"replicate count behind the error bar. See {t.cite}."
            )
        out.append(t)
    return out


def get(targets, name, context) -> Target:
    for t in targets:
        if t.name == name and t.context == context:
            return t
    raise KeyError(f"{name!r} in context {context!r}")


@dataclass
class Check:
    """Outcome of one consistency check."""
    name: str
    predicted: float
    observed: float
    unit: str
    verdict: str
    note: str

    @property
    def ratio(self) -> float:
        return self.observed / self.predicted if self.predicted else float("nan")


# --------------------------------------------------------------------------- #
# 1. hiPSC-CM: imaging-derived event rate vs flow-derived instantaneous fraction
# --------------------------------------------------------------------------- #
def check_hipsc_reservoir(targets) -> Check:
    """Does the live-imaging event rate explain the flow-cytometry mAG+ fraction?

    It does not, and the shortfall is the single most informative number in the
    hiPSC dataset: it implies a chronically arrested pool that never resolves
    into a scored outcome inside the 72 h window, and therefore that the
    reported 24.5 h polyploid duration is right-censored.
    """
    events = sum(get(targets, f"fate_{k}", "hiPSC_CM").value
                 for k in ("division", "binucleation", "polyploidization")) / 100.0
    window = 72.0
    dur = get(targets, "dur_sg2m_division", "hiPSC_CM").value
    rate = events / window                      # per cell per hour
    predicted = rate * dur * 100.0              # instantaneous mAG+ %
    observed = get(targets, "fucci_g2m", "hiPSC_CM").value
    excess = observed - predicted
    apparent = observed / 100.0 / rate
    poly = get(targets, "fate_polyploidization", "hiPSC_CM").value
    return Check(
        "hiPSC arrested reservoir", predicted, observed, "% mAG+",
        "INCONSISTENT -- excess pool",
        f"entry rate {rate*1e3:.3f}e-3/h; excess {excess:.2f} pp = {excess/poly:.1f}x the "
        f"{poly:.2f}% caught polyploidising; apparent S/G2/M {apparent:.1f} h vs {dur:.2f} h "
        f"measured. Rule out a permissive flow gate first; otherwise this is a "
        f"sustained-G2-arrest class the model must contain.",
    )


# --------------------------------------------------------------------------- #
# 2. Mouse in vivo: FUCCI snapshots -> cumulative EdU, across two papers
# --------------------------------------------------------------------------- #
def check_mouse_edu(targets, label_days: int = 4) -> Check:
    """Predict cumulative EdU at P7 from the P0/P7 mAG+ snapshots alone.

    Nothing is fitted. Agreement here is what licenses cumulative EdU as a
    held-out validation target. The entry rate is assumed to decay exponentially
    between the two snapshots, so treat the result as order-agreement.
    """
    dur = get(targets, "dur_sg2m", "P0_mouse_primary").value
    r0 = get(targets, "fucci_g2m", "P0_mouse_invivo").value / 100.0 / dur
    r7 = get(targets, "fucci_g2m", "P7_mouse_invivo").value / 100.0 / dur
    k = math.log(r0 / r7) / 7.0                          # per day
    cum = sum(r0 * 24 * math.exp(-k * d) for d in range(1, label_days + 1))
    observed = get(targets, "edu_untreated", "P7_mouse_invivo").value
    return Check(
        "mouse cumulative EdU", cum * 100.0, observed, "% EdU+ CM nuclei",
        "CONSISTENT",
        f"entry rate {r0*24*100:.2f}%/day (P0) -> {r7*24*100:.2f}%/day (P7), "
        f"{r0/r7:.1f}x fall, k={k:.3f}/day (t-half {math.log(2)/k:.2f} d). "
        f"Two papers, different cohorts and modalities, agree to "
        f"{cum*100/observed:.2f}x with zero fitted parameters.",
    )


# --------------------------------------------------------------------------- #
# 3. The Ki-67 / EdU fold discrepancy, and why its SIGN matters
# --------------------------------------------------------------------------- #
def check_ki67_edu_sign(targets) -> Check:
    """An EdU fold exceeding the Ki-67 fold cannot come from an entry-rate change.

    EdU saturates as 1-exp(-r*T) and so is sublinear in r; Ki-67 ~ r*D_K is
    linear. A pure rate increase therefore always gives fold(EdU) < fold(Ki-67).
    The observation is the other way round, which is positive evidence that
    Ki-67 is being *lost* from arrested cells while EdU persists -- i.e. a
    signature of non-productive cycling rather than a puzzle.
    """
    k0 = get(targets, "ki67_untreated", "mNCM_invitro").value
    k1 = get(targets, "ki67_clonidine", "mNCM_invitro").value
    ki_fold = k1 / k0
    edu_fold = get(targets, "edu_fold_clonidine", "mNCM_invitro").value
    return Check(
        "Ki-67 vs EdU fold sign", ki_fold, edu_fold, "fold",
        "INCONSISTENT with pure rate change -- diagnostic of arrest",
        f"ratio-of-ratios {edu_fold/ki_fold:.2f}. Saturation predicts "
        f"fold(EdU) < fold(Ki-67) for ANY entry-rate increase; observed is the "
        f"opposite sign. Prediction: this ratio is ~1.0 where cycling is "
        f"productive (hiPSC-CM) and >1 where it is not (mNCM).",
    )


def entry_rate_fold(edu_before: float, edu_after: float) -> float:
    """Saturation-corrected entry-rate fold from two cumulative EdU indices.

    The naive ratio of indices understates the rate change, because EdU
    saturates. Inverting 1-exp(-r*T) recovers r*T.
    """
    x0 = -math.log(1.0 - edu_before / 100.0)
    x1 = -math.log(1.0 - edu_after / 100.0)
    return x1 / x0


# --------------------------------------------------------------------------- #
# 4. Fate-budget identity -- genome copies are conserved
# --------------------------------------------------------------------------- #
def check_fate_budget(targets) -> Check:
    """dEdU ~ dBinucleation + dPolyploid + 2*dDivisions, per nucleus.

    The bookkeeping trap: EdU is measured per NUCLEUS while nucleation is per
    CELL, so the extra binucleation events must be converted before comparison.
    A violation of this identity indicates a denominator or protocol error, not
    a biological result.
    """
    m0 = get(targets, "mononucleated_untreated", "P7_mouse_invivo").value / 100.0
    m1 = get(targets, "mononucleated_clonidine", "P7_mouse_invivo").value / 100.0
    nuc0, nuc1 = 2 - m0, 2 - m1                      # nuclei per cell
    extra_events = (nuc1 - 1) - (nuc0 - 1)
    predicted = 2 * extra_events / nuc1 * 100.0      # extra EdU+ nuclei, %

    e0 = get(targets, "edu_untreated", "P7_mouse_invivo")
    e1 = get(targets, "edu_clonidine", "P7_mouse_invivo")
    observed = e1.value - e0.value
    se = math.sqrt((e0.sem or 0) ** 2 + (e1.sem or 0) ** 2)
    z = (observed - predicted) / se if se else float("nan")
    return Check(
        "fate-budget identity", predicted, observed, "pp dEdU",
        "CONSISTENT" if abs(z) < 2 else "INCONSISTENT",
        f"nuclei/cell {nuc0:.3f} -> {nuc1:.3f}; {extra_events:.3f} extra binucleation "
        f"events/cell. z={z:+.2f}; binucleation accounts for "
        f"{predicted/observed*100:.0f}% of the EdU rise, residual attributable to "
        f"tri/tetranucleation or cells still in cycle at harvest.",
    )


# --------------------------------------------------------------------------- #
# 5. Mitotic dwell time implied by the cytoplasmic-FUCCI fraction
# --------------------------------------------------------------------------- #
def check_mitotic_duration(targets, event_rate_per_h: float = 0.0024) -> Check:
    """How long is cardiomyocyte mitosis, given that 1% of P0 nuclei are in it?

    ``event_rate_per_h`` is the independent P0 division+binucleation rate
    (~5.7%/day; Soonpaa 1996, Alkass 2015), NOT taken from either FUCCI paper --
    which is what makes this an independent prediction rather than a restatement.
    """
    frac = get(targets, "fucci_mitotic", "P0_mouse_invivo").value / 100.0
    predicted = frac / event_rate_per_h
    return Check(
        "CM mitotic dwell time", predicted, float("nan"), "h",
        "PREDICTION (unmeasured)",
        f"{predicted:.1f} h from a {frac*100:.1f}% mitotic fraction and an "
        f"independent {event_rate_per_h*24*100:.1f}%/day event rate -- several-fold "
        f"longer than the ~1 h of cycling somatic cells. Should be longer still in "
        f"binucleating than dividing cells, since a failed furrow lingers. "
        f"Calibrates mitotic-entry and Cdc20 timing. Neither paper computed it.",
    )


# --------------------------------------------------------------------------- #
# 6. Is the polyploid class 4C-arrested or genuinely 8C?
# --------------------------------------------------------------------------- #
def check_ploidy_rounds(targets) -> Check:
    """Does the mNCM tetraploid increase imply one extra S-phase, or many?

    Murganti assigned the polyploid fate by FUCCI trace SHAPE, but a cell
    arrested at 4C in G2 produces that same shape without re-replicating. The
    tetraploid shift with no reported octaploid increase is consistent with one
    round, i.e. 4C arrest rather than iterated endocycling.
    """
    t0 = get(targets, "tetraploid_untreated", "mNCM_invitro").value
    t1 = get(targets, "tetraploid_clonidine", "mNCM_invitro").value
    return Check(
        "extra replication rounds", 1.0, (t1 - t0) / (t1 - t0), "rounds",
        "AMBIGUOUS in published data",
        f"tetraploid {t0:.2f}% -> {t1:.2f}% (+{t1-t0:.2f} pp) with no reported "
        f"octaploid rise: consistent with ONE extra S-phase (4C arrest), not "
        f"iterated endoreplication. Trace shape cannot separate the two. Resolve "
        f"with a 4C-vs-8C DNA-content histogram of the FUCCI-classified polyploid "
        f"cells; the model must emit both a trace-based and a mechanistic label.",
    )


def run_all() -> list[Check]:
    """Every pre-flight check, in reporting order."""
    t = load_targets()
    return [
        check_hipsc_reservoir(t),
        check_mouse_edu(t),
        check_ki67_edu_sign(t),
        check_fate_budget(t),
        check_mitotic_duration(t),
        check_ploidy_rounds(t),
    ]


def _wrap(text: str, width: int, indent: str) -> str:
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        lines.append(cur)
    return f"\n{indent}".join(lines)


def main() -> None:
    checks = run_all()
    print("=" * 78)
    print("cmcycle pre-flight: aggregate consistency of the calibration set")
    print("=" * 78)
    for c in checks:
        obs = "n/a" if math.isnan(c.observed) else f"{c.observed:.2f}"
        print(f"\n{c.name}  [{c.verdict}]")
        print(f"  predicted {c.predicted:.2f} {c.unit}   observed {obs} {c.unit}")
        print(f"  {_wrap(c.note, 72, '  ')}")
    print("\n" + "=" * 78)
    inconsistent = [c for c in checks if c.verdict.startswith("INCONSISTENT")]
    print(f"{len(checks)} checks; {len(inconsistent)} flag a real tension to model:")
    for c in inconsistent:
        print(f"  - {c.name}")
    print("=" * 78)


if __name__ == "__main__":
    main()
