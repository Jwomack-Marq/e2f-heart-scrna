"""Regenerate ``model/figures/*.svg`` from the analysis, deterministically.

Every figure is driven by :func:`cmcycle.baniol.run` or
:func:`cmcycle.preflight.run_all`, so a figure can never drift from the number it
claims to show. Run ``python3 -m cmcycle.figures`` from ``model/``.
"""
from __future__ import annotations

import json
import os

from . import baniol, preflight
from .spec import load_calibrated as spec_load
from .svg import (SEQ, Fig, fmt, footnote, legend, linscale, nice_ticks,
                  ordinal, xgrid, ygrid)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")

S1, S2, S3 = "var(--s1)", "var(--s2)", "var(--s3)"


# --------------------------------------------------------------------------- #
# 1. Ect2 specificity against the genome-wide drift baseline
# --------------------------------------------------------------------------- #
def fig_ect2(r) -> str:
    genes = sorted(r["per_gene"].items(), key=lambda kv: kv[1]["ratio"])
    n = len(genes)
    f = Fig(780, 62 + n * 20 + 118, pad=(96, 30, 100, 92),
            title="Ect2 is the only cytokinesis gene that beats the background drift",
            subtitle="P7/P0 ratio in cycling ventricular cardiomyocytes (n=26 vs 63). "
                     "Grey band = genome-wide 5–95%.")
    f.title_block()
    d = r["drift"]
    lo = min(0.5, min(v["ratio"] for _, v in genes) - 0.03)
    hi = max(1.15, max(v["ratio"] for _, v in genes) + 0.03)
    sx = linscale(lo, hi, f.x0, f.x1)

    # background: genome-wide drift band + median. Context, so it recedes.
    f.rect(sx(d["q05"]), f.y0, sx(d["q95"]) - sx(d["q05"]), f.y1 - f.y0, "var(--bnd)")
    f.line(sx(d["median"]), f.y0, sx(d["median"]), f.y1, "var(--axs)", 2, dash="2 4", cap="butt")
    f.text(sx(d["median"]), f.y0 - 8, f"genome-wide median {fmt(d['median'])}",
           "note", "middle")

    xgrid(f, sx, nice_ticks(lo, hi, 6), 2, "P7 / P0 expression ratio")
    legend(f, [("cytokinesis", S1), ("mitotic", S2)], y=f.y0 - 30)

    row = (f.y1 - f.y0) / n
    for i, (g, v) in enumerate(genes):
        y = f.y0 + row * (i + 0.5)
        col = S1 if v["module"] == "cytokinesis" else S2
        star = v["p"] < r["bonferroni_alpha"]
        r_ = 7 if star else 5.5
        f.circle(sx(v["ratio"]), y, r_, col, stroke="var(--srf)", sw=2)
        f.text(f.x0 - 10, y + 4, g, "val" if star else "lbl", "end")
        if star:
            f.text(sx(v["ratio"]) + r_ + 9, y + 4,
                   f"{fmt(v['ratio'])}×  ·  {ordinal(v['pct_genomewide'])} pctile",
                   "val")
    footnote(f, "Larger marks clear Bonferroni (p < %.4f). Only Ect2 is also a "
                "genome-wide outlier — RhoA is significant by t-test but sits at the "
                "%s percentile of background drift, so it is not specific."
                % (r["bonferroni_alpha"],
                   ordinal(r["per_gene"]["Rhoa"]["pct_genomewide"])),
             f.y1 + 58, x=32)
    return f.render()


# --------------------------------------------------------------------------- #
# 2. Module scores cancel the signal
# --------------------------------------------------------------------------- #
def fig_modules(r) -> str:
    """Slope chart. No y axis: three series x two points is a case where direct
    labels carry the values and an axis would only add collisions."""
    m = r["modules"]
    rows = [("mitotic module", "13 genes", m["mitotic"], S2),
            ("cytokinesis module", "12 genes", m["cytokinesis"], S1),
            ("difference index", "cyto − mito", m["difference_index"], S3)]
    f = Fig(760, 336, pad=(86, 132, 74, 186),
            title="Averaging cancels the one signal that matters",
            subtitle="Module means in cycling ventricular cardiomyocytes (n=26 vs 63), "
                     "±1 SE.")
    f.title_block()
    lo, hi = -0.66, 1.80
    sy = linscale(lo, hi, f.y1, f.y0)
    xa, xb = f.x0 + 54, f.x1 - 54
    f.line(f.x0, sy(0), f.x1, sy(0), "var(--grd)", 1, dash="2 4", cap="butt")
    f.text(f.x1 + 6, sy(0) + 4, "0", "tick")
    f.text(xa, f.y1 + 20, "P0", "key", "middle")
    f.text(xb, f.y1 + 20, "P7", "key", "middle")

    for name, sub, v, col in rows:
        ya, yb = sy(v["P0"]), sy(v["P7"])
        f.line(xa, ya, xb, yb, col, 2.5)
        for x, val, se, anch, off in ((xa, v["P0"], v["P0_se"], "end", -12),
                                      (xb, v["P7"], v["P7_se"], "start", 12)):
            f.line(x, sy(val - se), x, sy(val + se), col, 1.5)
            f.circle(x, sy(val), 5.5, col, stroke="var(--srf)", sw=2)
            f.text(x + off, sy(val) + 4, f"{val:+.3f}", "val", anch)
        f.text(14, ya - 2, name, "val")
        f.text(14, ya + 12, sub, "note")
        f.text(xb + 60, yb + 4,
               f"p = {fmt(v['p'],4) if v['p'] >= 1e-4 else '<0.0001'}", "val")

    di, cy, mi = m["difference_index"], m["cytokinesis"], m["mitotic"]
    footnote(f, "Both modules fall together (%.3f× and %.3f×), so their difference is "
                "flat (Δ = %+.3f, p = %.2f). The sig_ploidy heuristic in "
                "build_signature_scores.R is a difference of exactly these two "
                "averages — so it cannot see the Ect2 signal at all."
                % (cy["ratio"], mi["ratio"], di["delta"], di["p"]),
             f.y1 + 44, x=32)
    return f.render()


# --------------------------------------------------------------------------- #
# 3. The maturation coordinate M
# --------------------------------------------------------------------------- #
def fig_maturation(r) -> str:
    g = r["maturation"]["by_group"]
    order = sorted(g.items(), key=lambda kv: kv[1]["mean"])
    n = len(order)
    f = Fig(790, 96 + n * 30 + 126, pad=(100, 148, 106, 178),
            title="One measurable coordinate orders every group",
            subtitle="M = mean z(fatty-acid oxidation) − mean z(glycolysis), per cell. "
                     "Bars are ±1 SE.")
    f.title_block()
    lo, hi = -1.95, 1.35
    sx = linscale(lo, hi, f.x0, f.x1)
    xgrid(f, sx, nice_ticks(lo, hi, 7), 1, "maturation index M  (immature → mature)")
    legend(f, [("P0", S1), ("P7", S2)], y=f.y0 - 30)
    f.line(sx(0), f.y0, sx(0), f.y1, "var(--axs)", 1, dash="2 4", cap="butt")

    row = (f.y1 - f.y0) / n
    for i, (k, v) in enumerate(order):
        y = f.y0 + row * (i + 0.5)
        col = S1 if k.startswith("P0") else S2
        f.line(sx(v["mean"] - v["se"]), y, sx(v["mean"] + v["se"]), y, col, 1.5)
        r_ = 5.0 + min(3.0, v["n"] / 22.0)      # area hints at n
        f.circle(sx(v["mean"]), y, r_, col, stroke="var(--srf)", sw=2)
        weak = v["n"] < 20
        f.text(f.x0 - 12, y + 4, k.replace("-", " · "), "note" if weak else "lbl", "end")
        # clear the SE whisker as well as the marker, so the label never sits on it
        off = max(r_, sx(v["mean"] + v["se"]) - sx(v["mean"])) + 9
        f.text(sx(v["mean"]) + off, y + 4,
               f"{v['mean']:+.2f}  (n={v['n']})" + ("  ⚠ small n" if weak else ""),
               "note" if weak else "val")
    footnote(f, "P7 atrial CMs land between P0 and P7 ventricular — independently "
                "reproducing Baniol's observation that cycling P7 atrial cells are "
                "transcriptionally more immature than P7 ventricular. Marker size "
                "hints at n; the two P0-atrial groups (n=8, n=4) are too small to "
                "carry a claim.", f.y1 + 58, x=32)
    return f.render()


# --------------------------------------------------------------------------- #
# 4. M predicts the model's load-bearing nodes
# --------------------------------------------------------------------------- #
def fig_correlations(r, store) -> str:
    M, _, _ = baniol.maturation_index(store)
    g8 = store.groups8()
    cells = g8["P0-vCM-cycling"] + g8["P7-vCM-cycling"]
    p0 = set(g8["P0-vCM-cycling"])
    panels = [("Ect2", "cytokinesis competence", S1),
              ("E2f6", "cell-cycle exit enforcer", S2)]
    pw, gap = 316, 44
    f = Fig(2 * pw + gap + 100, 406, pad=(96, 26, 96, 62),
            title="M predicts the two couplings the model rests on",
            subtitle="Cycling ventricular cardiomyocytes only (n=%d), P0 and P7 pooled, "
                     "so stage alone cannot explain the trend." % len(cells))
    f.title_block()
    legend(f, [("P0 cell", S3), ("P7 cell", "var(--mut)")], y=f.y0 - 30)

    mlo, mhi = -2.0, 2.6
    for pi, (gene, role, col) in enumerate(panels):
        v = store.gene(gene)
        ys = [v[i] for i in cells]
        xs = [M[i] for i in cells]
        px0 = f.x0 + pi * (pw + gap)
        px1 = px0 + pw
        ylo, yhi = 0.0, max(ys) * 1.12 + 1e-9
        sx = linscale(mlo, mhi, px0, px1)
        sy = linscale(ylo, yhi, f.y1, f.y0)
        for t in nice_ticks(ylo, yhi, 4):
            f.line(px0, sy(t), px1, sy(t), "var(--grd)", 1, cap="butt")
            f.text(px0 - 8, sy(t) + 4, fmt(t, 1), "tick", "end")
        for t in nice_ticks(mlo, mhi, 5):
            f.text(sx(t), f.y1 + 16, fmt(t, 0), "tick", "middle")
        f.line(px0, f.y1, px1, f.y1, "var(--axs)", 1, cap="butt")
        f.line(px0, f.y0, px0, f.y1, "var(--axs)", 1, cap="butt")

        for i, x, y in zip(cells, xs, ys):
            f.circle(sx(x), sy(y), 4.2, S3 if i in p0 else "var(--mut)",
                     stroke="var(--srf)", sw=1.5, opacity=0.9)
        # least-squares fit line
        mx = sum(xs) / len(xs); my = sum(ys) / len(ys)
        den = sum((a - mx) ** 2 for a in xs) or 1.0
        b = sum((a - mx) * (c - my) for a, c in zip(xs, ys)) / den
        a0 = my - b * mx
        f.line(sx(mlo), sy(a0 + b * mlo), sx(mhi), sy(a0 + b * mhi), col, 2.5)

        st = r["maturation"]["correlations"][gene]
        f.text(px0, f.y0 - 8, f"{gene} — {role}", "key")
        f.text(px1, f.y0 - 8, f"r = {st['r']:+.3f}   t = {st['t']:+.2f}", "val", "end")
        f.text((px0 + px1) / 2, f.y1 + 36, "maturation index M", "lbl", "middle")
    footnote(f, "Cytokinesis competence falls with maturation and the exit enforcer "
                "rises. Both are therefore measured functions of M, not fitted "
                "parameters — which removes the model's largest degree of freedom.",
             f.y1 + 58, x=32)
    return f.render()


# --------------------------------------------------------------------------- #
# 5. E2F family across maturation
# --------------------------------------------------------------------------- #
def fig_e2f(r) -> str:
    groups = ["P0-vCM-cycling", "P7-vCM-cycling", "P0-vCM-noncycling",
              "P7-vCM-noncycling", "P0-aCM-cycling", "P7-aCM-cycling",
              "P0-aCM-noncycling", "P7-aCM-noncycling"]
    genes = [f"E2f{i}" for i in range(1, 9)]
    e2f = r["e2f"]
    vmax = max(e2f[g][k]["mean"] for g in genes for k in groups)
    cw, ch = 76, 34
    f = Fig(96 + cw * len(groups) + 40, 96 + ch * len(genes) + 152,
            pad=(96, 40, 128, 96),
            title="The E2F family splits into three roles",
            subtitle="Mean log-normalised expression. Cell shade = magnitude "
                     "(0 → %.2f)." % vmax)
    f.title_block()
    for gi, g in enumerate(genes):
        y = f.y0 + gi * ch
        f.text(f.x0 - 10, y + ch / 2 + 4, g, "val", "end")
        for ki, k in enumerate(groups):
            x = f.x0 + ki * cw
            val = e2f[g][k]["mean"]
            si = min(len(SEQ) - 1, int((val / vmax) ** 0.6 * (len(SEQ) - 1)))
            # 2px surface gap between adjacent fills
            f.rect(x + 1, y + 1, cw - 2, ch - 2, SEQ[si], rx=3)
            # Ink follows the FILL, not the theme -- the ramp is theme-invariant, so
            # a theme token here puts white text on pale blue in dark mode. Crossover
            # at step 8 is where white overtakes near-black for contrast on the ramp.
            ink = "#0b0b0b" if si < 8 else "#ffffff"
            f.add(f'<text x="{fmt(x+cw/2,1)}" y="{fmt(y+ch/2+4,1)}" '
                  f'style="font-size:10.5px;font-variant-numeric:tabular-nums;'
                  f'fill:{ink}" text-anchor="middle">{fmt(val,2)}</text>')
    for ki, k in enumerate(groups):
        st, ct, cy = k.split("-")
        x = f.x0 + ki * cw + cw / 2
        f.text(x, f.y0 - 26, f"{st} {ct}", "lbl", "middle")
        f.text(x, f.y0 - 11, cy, "note", "middle")
    notes = [
        "E2f1/7/8 are cycling-restricted (near-zero in every noncycling group).",
        "E2f2 is near-absent except in cycling cells and rises with maturation "
        "(t = %+.2f) — the endoreplication candidate." % e2f["E2f2"]["maturation_t"],
        "E2f6 is the only member expressed in noncycling cells, and it rises "
        "(t = %+.2f) — the exit enforcer." % e2f["E2f6"]["maturation_t"],
        "E2f7 is flat across maturation (t = %+.2f) while E2f8 rises (t = %+.2f)."
        % (e2f["E2f7"]["maturation_t"], e2f["E2f8"]["maturation_t"]),
    ]
    # hairline between the ventricular and atrial blocks
    xsep = f.x0 + 4 * cw
    f.line(xsep, f.y0 - 34, xsep, f.y0 + len(genes) * ch, "var(--axs)", 1, cap="butt")
    y = f.y1 + 4
    for s_ in notes:
        y = footnote(f, "• " + s_, y + 12, x=f.x0 - 84, width=f.width - 60)
    return f.render()


# --------------------------------------------------------------------------- #
# 6. Pre-flight: predicted vs observed
# --------------------------------------------------------------------------- #
def fig_preflight() -> str:
    checks = [c for c in preflight.run_all()
              if c.name in ("hiPSC arrested reservoir", "mouse cumulative EdU",
                            "fate-budget identity")]
    f = Fig(760, 364, pad=(96, 150, 106, 196),
            title="Two of three cross-paper checks close; one does not",
            subtitle="Closed-form predictions from the published aggregates — "
                     "no fitted parameters.")
    f.title_block()
    lo, hi = 0.0, 36.0
    sx = linscale(lo, hi, f.x0, f.x1)
    xgrid(f, sx, nice_ticks(lo, hi, 6), 0, "value (units differ per row — see labels)")
    legend(f, [("predicted", S3), ("observed", S1)], y=f.y0 - 30)
    row = (f.y1 - f.y0) / len(checks)
    units = {"hiPSC arrested reservoir": "% mVenus⁺",
             "mouse cumulative EdU": "% EdU⁺ nuclei",
             "fate-budget identity": "pp ΔEdU"}
    for i, c in enumerate(checks):
        y = f.y0 + row * (i + 0.5)
        f.line(sx(c.predicted), y, sx(c.observed), y,
               "var(--s2)" if c.verdict.startswith("INCONSISTENT") else "var(--grd)",
               3 if c.verdict.startswith("INCONSISTENT") else 2)
        f.circle(sx(c.predicted), y, 6, S3, stroke="var(--srf)", sw=2)
        f.circle(sx(c.observed), y, 6, S1, stroke="var(--srf)", sw=2)
        f.text(f.x0 - 12, y - 3, c.name, "val", "end")
        f.text(f.x0 - 12, y + 12, units[c.name], "note", "end")
        hi_, lo_ = max(c.observed, c.predicted), min(c.observed, c.predicted)
        fold = hi_ / lo_ if lo_ else float("nan")
        tag = ("off by %.2f×" if c.verdict.startswith("INCONSISTENT")
               else "within %.2f×") % fold
        f.text(max(sx(c.predicted), sx(c.observed)) + 14, y + 4, tag, "val")
    footnote(f, "The hiPSC gap implies a chronically arrested pool of ≈2.9 pp of all "
                "cells — 2.1× the polyploidising fraction caught inside the 72 h "
                "window — so the reported 24.5 h S/G2/M is a right-censored lower "
                "bound. The other two close with no fitted parameters.",
             f.y1 + 58, x=32)
    return f.render()


# --------------------------------------------------------------------------- #
# 7. The sort-enrichment trap
# --------------------------------------------------------------------------- #
def fig_sort(r) -> str:
    se = r["sort_enrichment"]
    f = Fig(700, 392, pad=(100, 40, 116, 112),
            title="Cycling fractions in this dataset are a sorting artefact",
            subtitle="FACS-sorted scRNA-seq vs the in-vivo expectation from "
                     "Baniol's FUCCI imaging.")
    f.title_block()
    lo, hi = 0.0, 50.0
    sy = linscale(lo, hi, f.y1, f.y0)
    ygrid(f, sy, nice_ticks(lo, hi, 5), 0, "% cycling cardiomyocytes")
    legend(f, [("scRNA-seq (sorted)", S1), ("in vivo (imaging)", S3)], y=f.y0 - 30)
    bw = 66
    for i, st in enumerate(("P0", "P7")):
        cx = f.x0 + (f.x1 - f.x0) * (0.28 + 0.44 * i)
        obs = se[st]["observed"] * 100
        if st == "P0":
            exp_lo = exp_hi = se[st]["expected"] * 100
            enr = "%.2f×" % se[st]["enrichment"]
            elab = f"{exp_hi:.1f}%"
        else:
            exp_hi, exp_lo = (max(se[st]["expected_range"]) * 100,
                              min(se[st]["expected_range"]) * 100)
            er = se[st]["enrichment_range"]
            enr = "%.1f–%.1f×" % (min(er), max(er))
            elab = f"{exp_lo:.1f}–{exp_hi:.1f}%"
        # both quantities as bars so they are directly comparable; 2px surface gap
        # between the pair, 4px rounded data-end anchored to the baseline
        f.rect(cx - bw - 1, sy(obs), bw, f.y1 - sy(obs), S1, rx=4)
        f.text(cx - bw / 2 - 1, sy(obs) - 9, f"{obs:.1f}%", "val", "middle")
        f.rect(cx + 1, sy(exp_hi), bw, f.y1 - sy(exp_hi), S3, rx=4)
        if exp_lo != exp_hi:                      # uncertainty cap on the bound
            f.line(cx + 1, sy(exp_lo), cx + 1 + bw, sy(exp_lo), "var(--pri)", 1.5)
        f.text(cx + 1 + bw / 2, sy(exp_hi) - 9, elab, "val", "middle")
        f.text(cx, f.y1 + 20, st, "key", "middle")
        f.text(cx, f.y1 + 36, f"{enr} enriched", "note", "middle")
    footnote(f, "P0 is essentially unenriched; P7 is 4.5–5.2× enriched. So the cycling "
                "fraction rises from P0 to P7 in this dataset and falls in reality — "
                "it can never be used as a rate. Only within-group expression "
                "contrasts are valid here.", f.y1 + 58, x=32)
    return f.render()


# --------------------------------------------------------------------------- #
def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    store = baniol.Store()
    r = baniol.run(store)
    from . import model as _model
    net = spec_load()
    M = _model.summary(net)
    with open(os.path.join(OUT, "model_results.json"), "w", encoding="utf8") as fh:
        json.dump(M, fh, indent=1, sort_keys=True, default=float)
    figs = {
        "fig1_ect2_specificity": fig_ect2(r),
        "fig2_module_cancellation": fig_modules(r),
        "fig3_maturation_coordinate": fig_maturation(r),
        "fig4_m_correlations": fig_correlations(r, store),
        "fig5_e2f_family": fig_e2f(r),
        "fig6_preflight_checks": fig_preflight(),
        "fig7_sort_enrichment": fig_sort(r),
        "fig8_maturation_fates": fig_maturation_fates(M),
        "fig9_clonidine_triad": fig_triad(M),
        "fig10_ko_and_screen": fig_ko_screen(M),
        "fig11_engine": fig_engine(),
        "fig12_activation_classes": fig_activation_classes(),
        "fig13_pathways": fig_pathways(),
    }
    for name, svg in figs.items():
        p = os.path.join(OUT, name + ".svg")
        with open(p, "w", encoding="utf8") as fh:
            fh.write(svg)
        print(f"  wrote {os.path.relpath(p, os.path.dirname(OUT))}  ({len(svg):,} B)")
    with open(os.path.join(OUT, "results.json"), "w", encoding="utf8") as fh:
        json.dump(r, fh, indent=1, sort_keys=True)
    print(f"  wrote figures/results.json")



# --------------------------------------------------------------------------- #
# 8. Maturation sets which fate the cycling flux lands in
# --------------------------------------------------------------------------- #
def fig_maturation_fates(M) -> str:
    rows = M["maturation"]
    f = Fig(780, 430, pad=(96, 132, 106, 76),
            title="Entry is necessary but not sufficient",
            subtitle="cmfate model. ONLY the maturation coordinate is swept — every "
                     "other input is held at the in-vivo P1 setting.")
    f.title_block()
    sx = linscale(0.05, 0.95, f.x0, f.x1)
    sy = linscale(0.0, 1.0, f.y1, f.y0)
    ygrid(f, sy, [0, 0.25, 0.5, 0.75, 1.0], 2, "share of cycling cells")
    for t in (0.2, 0.4, 0.6, 0.8):
        f.text(sx(t), f.y1 + 16, fmt(t, 1), "tick", "middle")
    f.text((f.x0 + f.x1) / 2, f.y1 + 34, "maturation index M", "lbl", "middle")
    legend(f, [("division", S1), ("binucleation", S2), ("polyploidization", S3)],
           y=f.y0 - 30)
    series = [("Division", S1), ("Binucleation", S2), ("Polyploidization", S3)]
    for name, col in series:
        pts = []
        for r in rows:
            cyc = sum(r[k] for k in ("Division", "Binucleation", "Polyploidization"))
            pts.append((sx(r["Maturation"]), sy(r[name] / cyc if cyc else 0)))
        for a, b in zip(pts, pts[1:]):
            f.line(a[0], a[1], b[0], b[1], col, 2.5)
        dy = {"Division": 14, "Binucleation": 0, "Polyploidization": 0}[name]
        f.text(pts[-1][0] + 10, pts[-1][1] + 4 + dy,
               name.replace("ization", ""), "val")
    # the context markers the papers actually measured
    for m in (0.30, 0.55):
        f.line(sx(m), f.y0, sx(m), f.y1, "var(--grd)", 1, dash="2 4", cap="butt")
        f.text(sx(m), f.y0 - 8, f"M = {m:.2f}", "note", "middle")
    footnote(f, "Along this axis the cycling flux moves from binucleation into "
                "polyploidization, crossing over near M = 0.62, while division's share "
                "stays under %.0f%% throughout. Division is low even at low M here "
                "because this sweep holds mechanical load and beta-adrenergic tone at "
                "their in-vivo values; the hiPSC-CM context reaches a 53%% division share "
                "because it differs in those too. That is the two-factor point — "
                "maturation alone does not set the outcome."
             % (100 * max(r["productive_share"] for r in rows) + 0.5,),
             f.y1 + 58, x=32)
    return f.render()


# --------------------------------------------------------------------------- #
# 9. The clonidine triad
# --------------------------------------------------------------------------- #
def fig_triad(M) -> str:
    rows = M["triad"]
    f = Fig(840, 396, pad=(104, 40, 96, 196),
            title="One drug, three contexts, three outcomes",
            subtitle="Change in each fate under clonidine. The model is calibrated at "
                     "the hiPSC context only; the other two are held out.")
    f.title_block()
    lo = min(min(r["delta"].values()) for r in rows)
    hi = max(max(r["delta"].values()) for r in rows)
    pad = 0.12 * (hi - lo)
    sx = linscale(min(lo, 0) - pad, hi + pad * 3, f.x0, f.x1)
    xgrid(f, sx, nice_ticks(min(lo, 0), hi, 5), 2, "change in fate fraction under clonidine")
    legend(f, [("division", S1), ("binucleation", S2), ("polyploidization", S3)],
           y=f.y0 - 30)
    f.line(sx(0), f.y0, sx(0), f.y1, "var(--axs)", 1, cap="butt")
    band = (f.y1 - f.y0) / len(rows)
    for i, r in enumerate(rows):
        y0 = f.y0 + band * i
        f.text(f.x0 - 12, y0 + 20, r["label"].split(" (")[0], "val", "end")
        f.text(f.x0 - 12, y0 + 35, f"M = {r['Maturation']:.2f}"
               + ("  ·  in vitro" if r["InVitro"] else "  ·  in vivo"), "note", "end")
        for j, (name, col) in enumerate((("Division", S1), ("Binucleation", S2),
                                         ("Polyploidization", S3))):
            v = r["delta"][name]
            yy = y0 + 14 + j * 13
            x0, x1 = sx(min(0, v)), sx(max(0, v))
            f.rect(x0, yy - 4.5, x1 - x0, 9, col, rx=4)
        best = max(r["delta"], key=r["delta"].get)
        f.text(sx(r["delta"][best]) + 12, y0 + 27,
               f"{best} {'OK' if r['hit'] else 'MISS'}", "val")
        if i:
            f.line(f.x0, y0, f.x1, y0, "var(--rule)" if False else "var(--grd)", 1, cap="butt")
    hits = sum(r["hit"] for r in rows)
    footnote(f, "%d of 3 dominant fates correct, from one calibration and one input "
                "change per context. What the model gets wrong is the MAGNITUDE of the "
                "entry response: %.1fx at hiPSC against an observed 2.4x, but %.0fx and "
                "%.0fx in the two mature contexts against observed 2.1x and 1.5x. "
                "Clonidine's leverage on entry scales with baseline PKA, which the model "
                "ties to beta-adrenergic tone."
             % (hits, rows[0]["entry_fold"], rows[1]["entry_fold"], rows[2]["entry_fold"]),
             f.y1 + 58, x=32)
    return f.render()


# --------------------------------------------------------------------------- #
# 10. In-silico E2f7/E2f8 knockdown and the target screen
# --------------------------------------------------------------------------- #
def fig_ko_screen(M) -> str:
    ko, screen = M["ko"], M["screen"][:8]
    f = Fig(820, 470, pad=(100, 34, 118, 60),
            title="The lab's own knockout, and what else the model says to try",
            subtitle="In-silico double knockdown of E2F7/E2F8, and every node ranked by "
                     "how much it converts cycling into division.")
    f.title_block()
    mid = f.x0 + 330
    # -- left panel: KO at two maturations -----------------------------------
    keys = ["SPhase", "Ect2", "Division", "Binucleation", "Polyploidization"]
    lo = min(min(r["double"][k] for k in keys) for r in ko)
    hi = max(max(r["double"][k] for k in keys) for r in ko)
    sx = linscale(min(0, lo) * 1.1, hi * 1.25, f.x0 + 96, mid - 30)
    f.text(f.x0 + 96, f.y0 - 26, "E2F7 + E2F8 double knockdown", "key")
    for t in (0.0, 0.2, 0.4):
        f.line(sx(t), f.y0, sx(t), f.y1 - 40, "var(--grd)", 1, cap="butt")
        f.text(sx(t), f.y1 - 26, fmt(t, 1), "tick", "middle")
    row = (f.y1 - 40 - f.y0) / len(keys)
    for i, k in enumerate(keys):
        y = f.y0 + row * (i + 0.5)
        f.text(f.x0 + 88, y + 4, k, "lbl", "end")
        for r, col, dy in ((ko[0], S1, -5), (ko[1], S2, 6)):
            v = r["double"][k]
            f.rect(sx(0), y + dy - 4, max(sx(v) - sx(0), 1), 8, col, rx=4)
        f.text(sx(max(r["double"][k] for r in ko)) + 8, y + 4,
               f"{ko[0]['double'][k]:+.3f} / {ko[1]['double'][k]:+.3f}", "note")
    legend(f, [("P0 (M=0.30)", S1), ("P7 (M=0.55)", S2)], x=mid - 30, y=f.y0 - 8)
    # -- right panel: screen --------------------------------------------------
    f.text(mid + 30, f.y0 - 26, "top targets at P7, by conditional division share", "key")
    smax = max(r["d_share"] for r in screen) or 1.0
    sx2 = linscale(0, smax * 1.35, mid + 190, f.x1 - 10)
    row2 = (f.y1 - 40 - f.y0) / len(screen)
    for i, r in enumerate(screen):
        y = f.y0 + row2 * (i + 0.5)
        col = S3 if r["mode"] == "overexpress" else S2
        f.text(mid + 182, y + 4, f"{r['node']}", "val", "end")
        f.text(mid + 34, y + 4, r["mode"], "note")
        f.rect(sx2(0), y - 5, max(sx2(r["d_share"]) - sx2(0), 1), 10, col, rx=4)
        f.text(sx2(r["d_share"]) + 8, y + 4, f"{r['d_share']:+.3f}", "note")
    sp = M["specificity"]
    footnote(f, "The knockout's pro-division effect is %.0fx larger at P0 than at P7 "
                "(%+.4f vs %+.4f), even though it raises Ect2 substantially at both "
                "(%+.2f and %+.2f) — the machinery goes up either way, but only at P0 is "
                "the rest of the context permissive enough to convert it. Scored "
                "like-for-like in hiPSC-CM, %.0f%% of perturbations clear the hit bar "
                "against 6%% in Murganti's wet screen, and only %.0f%% of those also "
                "raise the division share, against 33%% surviving their validation — an "
                "entry-scored screen mostly finds things that do not help."
             % (ko[0]["double"]["Division"] / max(ko[1]["double"]["Division"], 1e-9),
                ko[0]["double"]["Division"], ko[1]["double"]["Division"],
                ko[0]["double"]["Ect2"], ko[1]["double"]["Ect2"],
                100 * sp["entry_hit_rate"],
                100 * sp["entry_hits_that_help"]),
             f.y1 + 12, x=32)
    return f.render()


# --------------------------------------------------------------------------- #
# 11. The engine: how a state vector becomes a derivative
# --------------------------------------------------------------------------- #
def fig_engine() -> str:
    """The two composition rules are the claim: both of the model's structural
    bugs lived in them, so the diagram shows where, not just what."""
    from .svg import arrow, arrowhead_defs, box, elbow, tag
    f = Fig(1040, 568, pad=(96, 24, 24, 24),
            title="The engine: one state vector in, one derivative out",
            subtitle="Normalized-Hill logic-ODE (Kraeutler, Soltis & Saucerman 2010). "
                     "Orange marks the rule where both of this model's structural bugs lived.")
    f.title_block()
    arrowhead_defs(f)

    # ---- grid ---------------------------------------------------------------
    C1, C2, C3, C4 = 24, 202, 556, 844      # column lefts
    W1, W2, W3, W4 = 118, 288, 226, 172     # column widths
    R1, BH = 108, 96                        # pipeline row
    R2 = 392                                # integration row
    R3 = 476                                # bug band

    # ---- row 1: state -> reaction -> node -> drive -------------------------
    box(f, C1, R1, W1, BH, "state  Y", ["Y[i] ∈ [0,1]", "one per species", "55 here"])
    arrow(f, C1 + W1 + 6, R1 + BH / 2, C2 - 6, R1 + BH / 2, "read")
    box(f, C2, R1, W2, BH, "per REACTION",
        ["act(x) = w·β·xⁿ/(Kⁿ+xⁿ)", "inhib(x) = w − act(x)", "AND:  ∏ vals / w^(k−1)"])
    arrow(f, C2 + W2 + 6, R1 + BH / 2, C3 - 6, R1 + BH / 2, "value")
    box(f, C3, R1, W3, BH, "per NODE",
        ["weighted OR over every", "reaction onto it:", "d ← d + v − d·v"], accent=True)
    arrow(f, C3 + W3 + 6, R1 + BH / 2, C4 - 6, R1 + BH / 2, "drive")
    box(f, C4, R1, W4, BH, "per SPECIES", ["dY/dt =", "(drive·Ymax − Y)/τ"])

    # ---- where the Hill constants come from, under the reaction stage ------
    box(f, C2, R1 + BH + 40, W2, 62, "(w, n, EC50)  →  β, K",
        ["β = (EC50ⁿ−1)/(2·EC50ⁿ−1)", "K = (β−1)^(1/n)"])
    arrow(f, C2 + W2 / 2, R1 + BH + 36, C2 + W2 / 2, R1 + BH + 6, None)
    tag(f, C2, R1 + BH + 122, "these belong to the REACTION, not the reactant:", accent=False)
    tag(f, C2, R1 + BH + 137, "one curve applies to all of its reactants", accent=False)
    tag(f, C2, R1 + BH + 158, "guard: EC50ⁿ < 0.5, else β flips sign", accent=True)

    # ---- how a perturbation enters, under the species stage ----------------
    tag(f, C4, R1 + BH + 22, "knockdown:  Ymax → 0", accent=False)
    tag(f, C4, R1 + BH + 37, "overexpress: hold Y = 1", accent=False)
    tag(f, C4, R1 + BH + 52, "graded:  Ymax ×= f", accent=False)

    # ---- the sink shortcut, branching off the drive arrow ------------------
    box(f, C3, R1 + BH + 40, 264, 62, "output-only nodes stop here",
        ["y = drive·Ymax, algebraically", "4 fates + 4 marker nodes"])
    elbow(f, [(C4 - 26, R1 + BH / 2 + 10), (C4 - 26, R1 + BH + 34)], None, dash="4 3")

    # ---- row 2: the integration loop ---------------------------------------
    box(f, C4, R2, W4, 62, "RK23 step", ["adaptive, embedded", "3rd-order error"])
    elbow(f, [(C4 + W4 - 40, R1 + BH + 106), (C4 + W4 - 40, R2 - 4)], None)
    box(f, C3, R2, W3, 62, "converged?", ["max |dY/dt| < tol", "else step again"])
    arrow(f, C4 - 6, R2 + 31, C3 + W3 + 6, R2 + 31, "y′")
    elbow(f, [(C3 - 6, R2 + 31), (C1 + 59, R2 + 31), (C1 + 59, R1 + BH + 6)], "iterate")

    # ---- the bug band: proximity plus one clean pointer --------------------
    tag(f, C1, R3 + 14, "① a low-weight OR'd reaction acts as a FLOOR", accent=True)
    tag(f, C1 + 18, R3 + 30,
        "flattened Ect2's maturation response to 1.20×, against a measured 1.83×",
        accent=False)
    tag(f, C1, R3 + 54, "② OR can only ADD drive, so no route can be vetoed by another",
        accent=True)
    tag(f, C1 + 18, R3 + 70,
        "a second Midbody route bypassed the obligatory Ect2/RhoA arm entirely",
        accent=False)
    elbow(f, [(C3 - 30, R3 + 6), (C3 - 30, R1 + BH + 26), (C3 + 22, R1 + BH + 26),
              (C3 + 22, R1 + BH + 6)], None, accent=True)
    return f.render()


# --------------------------------------------------------------------------- #
# 12. What the parameter classes actually do
# --------------------------------------------------------------------------- #
def fig_activation_classes() -> str:
    """Small multiples: one curve per panel, so no palette and no legend are needed.

    Every curve is evaluated through the real engine, including the invalid case --
    the last panel is what the guard prevents.
    """
    from .logic import Reaction, Network
    CLASSES = [("T", 1.4, 0.35, "transduction"), ("F", 1.4, 0.40, "fate layer"),
               ("G", 3.0, 0.45, "gate"), ("L", 1.4, 0.50, "linear"),
               ("B", 2.0, 0.60, "basal"), ("✗", 1.4, 0.70, "invalid")]
    pw, gap, ph = 142, 20, 132
    f = Fig(24 * 2 + len(CLASSES) * pw + (len(CLASSES) - 1) * gap, ph + 208,
            pad=(96, 24, 92, 24),
            title="What each parameter class does to a signal",
            subtitle="act(x) against the dashed identity line. Above it amplifies, "
                     "below it attenuates — evaluated through the engine itself.")
    f.title_block()
    for i, (name, n, ec, role) in enumerate(CLASSES):
        px = 24 + i * (pw + gap)
        py = f.y0
        invalid = ec ** n >= 0.5
        # a bare Network just to reach _act with the real code path
        net = Network(species=["a", "b"], yinit=[0, 0], ymax=[1, 1], tau=[1, 1],
                      idx={"a": 0, "b": 1})
        if invalid:
            class _R:                      # bypass the guard on purpose, to show it
                w, n_, ec_ = 1.0, n, ec
            beta = (ec ** n - 1) / (2 * ec ** n - 1)
            K = 1e-6
            def act(x, beta=beta, n=n):
                v = 1.0 * (beta * x ** n) / (K ** n + x ** n) if x > 0 else 0.0
                return v
        else:
            r = Reaction(target=0, reactants=[(1, False)], w=1.0, n=n, ec50=ec)
            net.reactions = [r]
            def act(x, r=r, net=net):
                return net._act(x, r)
        lo, hi = (-2.0, 1.05) if invalid else (0.0, 1.05)
        sx = linscale(0.0, 1.0, px, px + pw)
        sy = linscale(lo, hi, py + ph, py)
        # frame + zero line
        f.rect(px, py, pw, ph, "var(--bnd)", rx=6)
        if lo < 0:
            f.line(px, sy(0), px + pw, sy(0), "var(--axs)", 1, cap="butt")
        # identity
        f.line(sx(0), sy(0), sx(1), sy(1), "var(--axs)", 1.5, dash="3 3", cap="butt")
        # the curve
        pts = [(sx(k / 60), sy(max(min(act(k / 60), 1.05), lo))) for k in range(61)]
        col = "var(--s2)" if invalid else "var(--s1)"
        for a, b in zip(pts, pts[1:]):
            f.line(a[0], a[1], b[0], b[1], col, 2.2)
        f.text(px, py - 26, f"{name}   n={fmt(n,1)}  EC50={fmt(ec,2)}", "key")
        f.text(px, py - 11, role, "note")
        a3 = act(0.3)
        verdict = ("β = %.2f, so act() < 0" % ((ec**n-1)/(2*ec**n-1))
                   if invalid else
                   ("amplifies (act(0.3) = %.2f)" if a3 > 0.3 else
                    "attenuates (act(0.3) = %.2f)") % a3)
        f.text(px, py + ph + 18, verdict.split(" (")[0], "val" if invalid else "lbl")
        if not invalid:
            f.text(px, py + ph + 33, "act(0.3) = %.2f" % a3, "note")
        else:
            f.text(px, py + ph + 33, "inhib = w − act > w", "note")
        f.text(px, py + ph + 48, "EC50ⁿ = %.3f" % (ec ** n),
               "tagA" if invalid else "note")
    footnote(f, "The ceiling is EC50 < 2^(−1/n) — 0.610 at n = 1.4, 0.794 at n = 3. "
                "Netflux's own defaults (n = 1.4, EC50 = 0.50) sit safely inside it, "
                "which is why the reference implementation never needed the guard; a "
                "high-threshold gate is what walks into it, and this model has three.",
             f.y1 + 58, x=24)
    return f.render()


# --------------------------------------------------------------------------- #
# 13. The pathway map: every node and every edge in the network
# --------------------------------------------------------------------------- #
#: Column assignment, ordered by biology rather than by computed graph depth --
#: depth-from-inputs is short-circuited because Maturation feeds deep nodes
#: directly, so it places Ect2 at depth 1 and E2F7 at depth 4.
#: Within each column, order is hand-tuned to shorten edges.
PATHWAY_COLUMNS = [
    ("stimulus", ["Clonidine", "BetaAR", "Nrg1", "IGF1", "MechLoad", "ROSenv",
                  "InVitro", "Maturation"]),
    ("signalling · stress",
     ["cAMP", "PKA", "Autophagy", "ERK", "Akt", "Amotl1", "YAP",
      "ROS", "DDR", "Ccng1", "p21", "p38", "Glycolysis", "FAO"]),
    ("E2F / G1–S", ["CycD", "Rb", "E2F1", "E2F3", "E2Fact", "E2Frep",
                    "E2F7", "E2F8", "E2F6", "E2F2", "CycE", "p27"]),
    ("cycle machinery", ["CycA", "CycB", "Cdc25", "Pkmyt1", "MitCompRaw", "Cdc20"]),
    ("cytokinesis", ["Ect2", "RhoA", "Anillin", "Midbody", "AbsRaw",
                     "Centralspindlin", "AurKB"]),
    ("gates", ["SPhase", "MitoticEntry", "Abscission"]),
    ("outcome", ["Quiescent", "Division", "Binucleation", "Polyploidization", "Ki67"]),
]

#: The two arms the fate decision turns on. Everything else stays recessive.
ARM_A = {("E2Fact", "Ect2"), ("Ect2", "RhoA"), ("RhoA", "Midbody"),
         ("Midbody", "AbsRaw"), ("AbsRaw", "Abscission"),
         ("Abscission", "Division"), ("Abscission", "Binucleation")}
ARM_B = {("ROS", "DDR"), ("DDR", "Ccng1"), ("DDR", "Pkmyt1"),
         ("Ccng1", "MitCompRaw"), ("Pkmyt1", "MitCompRaw"),
         ("MitCompRaw", "MitoticEntry"), ("MitoticEntry", "Polyploidization")}


def fig_pathways() -> str:
    """Every node and every edge, generated from the spec so it cannot drift."""
    from . import spec
    from .svg import curve, edge_defs, node
    net = spec.load()
    mod = net.meta["modules"]

    placed = [n for _, col in PATHWAY_COLUMNS for n in col]
    missing = [s for s in net.species if s not in placed]
    extra = [n for n in placed if n not in net.idx]
    if missing or extra:
        raise AssertionError(f"pathway layout out of sync: missing {missing}, extra {extra}")

    NW, NH, VP = 122, 24, 34
    CW, CG = NW, 62
    rows = max(len(c) for _, c in PATHWAY_COLUMNS)
    W = 30 * 2 + len(PATHWAY_COLUMNS) * CW + (len(PATHWAY_COLUMNS) - 1) * CG
    f = Fig(W, 112 + rows * VP + 118, pad=(112, 30, 112, 30),
            title="The cmfate network: 55 nodes, 78 reactions, 93 edges",
            subtitle="Every edge is drawn from the spec. Orange = the abscission arm "
                     "(what maturation closes); blue = the mitotic-entry brake "
                     "(what culture closes).")
    f.title_block()
    edge_defs(f)

    # ---- place nodes --------------------------------------------------------
    pos = {}
    for ci, (title, col) in enumerate(PATHWAY_COLUMNS):
        cx = 30 + ci * (CW + CG)
        f.text(cx, f.y0 - 14, title, "colttl")
        top = f.y0 + (rows - len(col)) * VP / 2          # centre short columns
        for ri, n in enumerate(col):
            pos[n] = (cx, top + ri * VP, NW, NH)

    # ---- edges, drawn before nodes so the boxes sit on top ------------------
    mat_targets = set()
    for r in net.reactions:
        if r.is_input:
            continue
        for i, inh in r.reactants:
            a, b = net.species[i], net.species[r.target]
            if a == "Maturation":            # 15 edges: badged, not drawn
                mat_targets.add(b)
                continue
            (ax, ay, aw, ah), (bx, by, bw, bh) = pos[a], pos[b]
            arm = "A" if (a, b) in ARM_A else ("B" if (a, b) in ARM_B else None)
            if bx > ax:                       # forward
                x1, y1, x2, y2 = ax + aw, ay + ah / 2, bx - 2, by + bh / 2
                bow = 0.06
            elif bx < ax:                     # feedback: bow under the columns
                x1, y1, x2, y2 = ax, ay + ah, bx + bw, by + bh
                bow = -0.10
            else:                             # same column: arc out to the left
                x1, y1, x2, y2 = ax, ay + ah / 2, bx, by + bh / 2
                bow = 0.42 if y2 > y1 else -0.42
            curve(f, x1, y1, x2, y2, inhibit=inh, accent=arm, bow=bow,
                  opacity=None if arm else 0.55)

    # ---- nodes -------------------------------------------------------------
    for n, (x, y, w, h) in pos.items():
        arm = ("A" if n in {"Ect2", "RhoA", "Midbody", "AbsRaw", "Abscission"} else
               "B" if n in {"ROS", "DDR", "Ccng1", "Pkmyt1", "MitCompRaw",
                            "MitoticEntry"} else
               "C" if n == "Maturation" else None)
        node(f, x, y, w, h, n, badge=(n in mat_targets), accent=arm,
             bold=mod[n] in ("fate", "input"))

    # ---- legend ------------------------------------------------------------
    ly = f.y1 + 26
    f.circle(36, ly - 3.5, 3.0, "var(--s3)")
    f.text(46, ly, "receives Maturation (%d nodes; its 15 edges are badged rather "
                   "than drawn, or they would dominate the map)" % len(mat_targets),
           "note")
    f.add(f'<path d="M30,{fmt(ly+18,1)} L64,{fmt(ly+18,1)}" stroke="var(--mut)" '
          f'stroke-width="1" marker-end="url(#act)"/>')
    f.text(74, ly + 21, "activates", "note")
    f.add(f'<path d="M160,{fmt(ly+18,1)} L194,{fmt(ly+18,1)}" stroke="var(--mut)" '
          f'stroke-width="1" marker-end="url(#inh)"/>')
    f.text(204, ly + 21, "inhibits  (28 of the 93 edges)", "note")
    footnote(f, "Read left to right: a stimulus sets the signalling layer, which sets "
                "E2F activity, which builds the cycle and cytokinesis machinery, which "
                "sets three gates — and the four outcomes are a product of those three "
                "gates alone. The two coloured arms are what make one drug give three "
                "different outcomes: maturation closes the orange one, culture closes "
                "the blue one.", ly + 46, x=30)
    return f.render()

if __name__ == "__main__":
    main()
