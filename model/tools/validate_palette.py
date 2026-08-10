"""Python twin of the dataviz skill's ``validate_palette.js``.

Same constants, same math, same verdicts -- there is no JS runtime on this
machine and the rule is to *compute* palette safety rather than eyeball it.

Checks (the ones measurable from hexes alone):
  2.  Lightness band     -- OKLCH L inside the mode's band
  3.  Chroma floor       -- OKLCH C >= floor, below which a hue reads as gray
  4.  CVD separation     -- OKLab dE x100 under simulated protanopia/deuteranopia
                            (Machado-Oliveira-Fernandes 2009, severity 1.0)
  4b. Normal-vision floor-- worst unsimulated OKLab dE x100 on the active pairlist
  5.  Contrast vs surface-- WCAG ratio of each mark against the chart surface

Checks 1 (fixed hue order) and 6 (values come from the documented palette) are
structural and enforced by review, not computable here.

Usage:
    python3 validate_palette.py "#2a78d6,#eb6834,#1baf7a" --mode light --pairs all
Exit code 1 on any hard FAIL; WARN bands exit 0.
"""
from __future__ import annotations

import math
import re
import sys

BAND = {"light": (0.43, 0.77), "dark": (0.48, 0.67)}   # OKLCH L
CHROMA_FLOOR = 0.10
CVD_TARGET, CVD_FLOOR = 8.0, 6.0                       # OKLab dE x100, min(protan, deutan)
NORMAL_FLOOR = 15.0
CONTRAST_MIN = 3.0
DEFAULT_SURFACE = {"light": "#fcfcfb", "dark": "#1a1a19"}

MACHADO = {
    "protan": ((0.152286, 1.052583, -0.204868),
               (0.114503, 0.786281, 0.099216),
               (-0.003882, -0.048116, 1.051998)),
    "deutan": ((0.367322, 0.860646, -0.227968),
               (0.280085, 0.672501, 0.047413),
               (-0.011820, 0.042940, 0.968881)),
    "tritan": ((1.255528, -0.076749, -0.178779),
               (-0.078411, 0.930809, 0.147602),
               (0.004733, 0.691367, 0.303900)),
}

# Kept in lockstep with the JS twin: the intersection of what JS trim() and
# Python str.strip() both remove, so a hex list pasted from a rendered page
# (NBSP / em-space padding) normalises identically in both.
_WS = "[ \t\n\v\f\r   -     　]+"
_WS_RE = re.compile(rf"^{_WS}|{_WS}$")
_HEX_RE = re.compile(r"^#?[0-9a-fA-F]{6}$")


def split_colors(raw: str) -> list[str]:
    return [c for c in (_WS_RE.sub("", p) for p in (raw or "").split(",")) if c]


def is_hex(v: str) -> bool:
    return bool(_HEX_RE.match(v))


def _srgb(h: str) -> tuple[float, float, float]:
    h = _WS_RE.sub("", h).lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


def _s2lin(c: float) -> float:
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def _lin(h: str) -> tuple[float, float, float]:
    return tuple(_s2lin(c) for c in _srgb(h))


def rel_lum(h: str) -> float:
    r, g, b = _lin(h)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast(a: str, b: str) -> float:
    hi, lo = sorted((rel_lum(a), rel_lum(b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


def _oklab_from_lin(rgb) -> tuple[float, float, float]:
    r, g, b = rgb
    l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b) ** (1 / 3)
    m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b) ** (1 / 3)
    s = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b) ** (1 / 3)
    return (0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
            1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
            0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s)


def oklch(h: str) -> tuple[float, float]:
    L, a, b = _oklab_from_lin(_lin(h))
    return L, math.hypot(a, b)


def _simulate(h: str, kind: str):
    r, g, b = _lin(h)
    M = MACHADO[kind]
    return tuple(min(1.0, max(0.0, M[i][0] * r + M[i][1] * g + M[i][2] * b))
                 for i in range(3))


def delta_e(h1: str, h2: str, kind: str | None = None) -> float:
    a = _oklab_from_lin(_simulate(h1, kind) if kind else _lin(h1))
    b = _oklab_from_lin(_simulate(h2, kind) if kind else _lin(h2))
    return 100 * math.dist(a, b)


def _pairlist(n: int, pairs: str):
    if pairs == "all":
        return [(i, j) for i in range(n) for j in range(i + 1, n)]
    return [(i, i + 1) for i in range(n - 1)]


def validate(palette, mode="light", surface=None, pairs="adjacent"):
    surface = surface or DEFAULT_SURFACE[mode]
    lo, hi = BAND[mode]
    report, ok = [], True

    off = [(c, round(oklch(c)[0], 3)) for c in palette
           if not (lo <= oklch(c)[0] <= hi)]
    ok &= not off
    report.append(("Lightness band", not off,
                   f"outside band: {off}" if off else
                   f"all {len(palette)} inside L {lo}-{hi}"))

    gray = [(c, round(oklch(c)[1], 3)) for c in palette if oklch(c)[1] < CHROMA_FLOOR]
    ok &= not gray
    report.append(("Chroma floor", not gray,
                   f"below {CHROMA_FLOOR}: {gray}" if gray else
                   f"all >= {CHROMA_FLOOR}"))

    pl = _pairlist(len(palette), pairs)
    worst_cvd, worst_cvd_pair, fails, warns = math.inf, None, [], []
    for i, j in pl:
        p = delta_e(palette[i], palette[j], "protan")
        d = delta_e(palette[i], palette[j], "deutan")
        v = min(p, d)
        if v < worst_cvd:
            worst_cvd, worst_cvd_pair = v, (palette[i], palette[j])
        if v < CVD_FLOOR:
            fails.append((palette[i], palette[j], round(v, 1)))
        elif v < CVD_TARGET:
            warns.append((palette[i], palette[j], round(v, 1)))
    cvd_ok = not fails
    ok &= cvd_ok
    msg = f"worst {pairs} CVD dE {worst_cvd:.1f} {worst_cvd_pair}"
    if fails:
        msg += f" | below floor {CVD_FLOOR}: {fails}"
    elif warns:
        msg += f" | WARN 6-8 band (needs secondary encoding): {warns}"
    report.append((f"CVD separation ({pairs})", cvd_ok, msg))

    worst_n, worst_n_pair = math.inf, None
    for i, j in pl:
        v = delta_e(palette[i], palette[j])
        if v < worst_n:
            worst_n, worst_n_pair = v, (palette[i], palette[j])
    n_ok = worst_n >= NORMAL_FLOOR
    ok &= n_ok
    report.append(("Normal-vision floor", n_ok,
                   f"worst {worst_n:.1f} {worst_n_pair} (floor {NORMAL_FLOOR})"))

    low = [(c, round(contrast(c, surface), 2)) for c in palette
           if contrast(c, surface) < CONTRAST_MIN]
    report.append(("Contrast vs surface", True,
                   f"WARN sub-{CONTRAST_MIN}:1 -> relief rule (visible labels or "
                   f"table view): {low}" if low else
                   f"all >= {CONTRAST_MIN}:1 vs {surface}"))
    return ok, report


def main(argv) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    raw = argv[1]
    mode = "light"
    surface = None
    pairs = "adjacent"
    for i, a in enumerate(argv):
        if a == "--mode" and i + 1 < len(argv):
            mode = argv[i + 1]
        if a == "--surface" and i + 1 < len(argv):
            surface = argv[i + 1]
        if a == "--pairs" and i + 1 < len(argv):
            pairs = argv[i + 1]
    palette = split_colors(raw)
    bad = [c for c in palette if not is_hex(c)]
    if not palette or bad:
        print(f"FAIL: not 6-digit hex: {bad or '(empty palette)'}")
        return 1
    palette = ["#" + c.lstrip("#") for c in palette]

    ok, report = validate(palette, mode=mode, surface=surface, pairs=pairs)
    print(f"palette ({len(palette)}) mode={mode} pairs={pairs} "
          f"surface={surface or DEFAULT_SURFACE[mode]}")
    print("-" * 78)
    for name, passed, detail in report:
        print(f"  {'PASS' if passed else 'FAIL'}  {name:26s} {detail}")
    print("-" * 78)
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
