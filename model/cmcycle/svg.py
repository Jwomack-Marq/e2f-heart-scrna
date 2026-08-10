"""Dependency-free SVG chart primitives.

No matplotlib, no numpy -- the figures must regenerate anywhere the analysis
runs. Output is an inline-able ``<svg>`` fragment carrying its own scoped
``<style>``, so it drops into an HTML page and follows the viewer's theme.

Theme handling: every colour is a CSS custom property defined three times --
on the bare scope (light), under ``prefers-color-scheme: dark`` guarded with
``:not([data-theme="light"])``, and under ``[data-theme="dark"]`` so an explicit
toggle wins in both directions. Dark steps are chosen for the dark surface, not
flipped.

Palette: the dataviz reference instance, validated by
``model/tools/validate_palette.py``. Categorical use is capped at the first three
slots because only those clear the all-pairs CVD and normal-vision floors in both
modes. Aqua is sub-3:1 on the light surface, so the relief rule applies -- every
series that uses it also carries a direct label.
"""
from __future__ import annotations

import html
import itertools
import math
from dataclasses import dataclass, field

_uid = itertools.count(1)

# Categorical slots 1-3 (light, dark). Validated all-pairs in both modes.
SERIES = [("#2a78d6", "#3987e5"),   # blue
          ("#eb6834", "#d95926"),   # orange
          ("#1baf7a", "#199e70")]   # aqua  (light contrast 2.74 -> label required)

# Sequential blue ramp, light->dark, for magnitude encoding.
SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
       "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

INK = {
    "surface": ("#fcfcfb", "#1a1a19"),
    "primary": ("#0b0b0b", "#ffffff"),
    "secondary": ("#52514e", "#c3c2b7"),
    "muted": ("#898781", "#898781"),
    "grid": ("#e1e0d9", "#2c2c2a"),
    "axis": ("#c3c2b7", "#383835"),
    "band": ("#f0efec", "#262624"),
}


def esc(s) -> str:
    return html.escape(str(s), quote=True)


def fmt(x: float, nd: int = 2) -> str:
    s = f"{x:.{nd}f}"
    return s.rstrip("0").rstrip(".") if "." in s else s


@dataclass
class Fig:
    """An SVG figure. Append shapes with the helpers, then ``render()``."""
    width: float
    height: float
    pad: tuple = (56, 18, 44, 132)      # top, right, bottom, left
    title: str = ""
    subtitle: str = ""
    body: list = field(default_factory=list)
    cls: str = field(default_factory=lambda: f"vz{next(_uid)}")

    # --- plot box ---------------------------------------------------------
    @property
    def x0(self): return self.pad[3]
    @property
    def x1(self): return self.width - self.pad[1]
    @property
    def y0(self): return self.pad[0]
    @property
    def y1(self): return self.height - self.pad[2]

    # --- primitives ------------------------------------------------------
    def add(self, s: str): self.body.append(s); return self

    def rect(self, x, y, w, h, fill, rx=0, opacity=None, stroke=None, sw=0):
        o = f' opacity="{opacity}"' if opacity is not None else ""
        st = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
        return self.add(f'<rect x="{fmt(x,1)}" y="{fmt(y,1)}" width="{fmt(max(w,0),1)}" '
                        f'height="{fmt(max(h,0),1)}" rx="{rx}" fill="{fill}"{o}{st}/>')

    def line(self, x1, y1, x2, y2, stroke, sw=2, dash=None, cap="round"):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        return self.add(f'<line x1="{fmt(x1,1)}" y1="{fmt(y1,1)}" x2="{fmt(x2,1)}" '
                        f'y2="{fmt(y2,1)}" stroke="{stroke}" stroke-width="{sw}" '
                        f'stroke-linecap="{cap}"{d}/>')

    def circle(self, cx, cy, r, fill, stroke=None, sw=2, opacity=None):
        st = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
        o = f' opacity="{opacity}"' if opacity is not None else ""
        return self.add(f'<circle cx="{fmt(cx,1)}" cy="{fmt(cy,1)}" r="{fmt(r,1)}" '
                        f'fill="{fill}"{st}{o}/>')

    def text(self, x, y, s, cls="lbl", anchor="start", dy=0):
        return self.add(f'<text x="{fmt(x,1)}" y="{fmt(y+dy,1)}" class="{cls}" '
                        f'text-anchor="{anchor}">{esc(s)}</text>')

    def title_block(self):
        self.text(self.pad[3] - 0 if False else 12, 22, self.title, "ttl")
        if self.subtitle:
            self.text(12, 40, self.subtitle, "sub")
        return self

    # --- render ----------------------------------------------------------
    def render(self) -> str:
        c = self.cls
        def css(mode: int) -> str:
            v = [f"--srf:{INK['surface'][mode]}", f"--pri:{INK['primary'][mode]}",
                 f"--sec:{INK['secondary'][mode]}", f"--mut:{INK['muted'][mode]}",
                 f"--grd:{INK['grid'][mode]}", f"--axs:{INK['axis'][mode]}",
                 f"--bnd:{INK['band'][mode]}"]
            v += [f"--s{i+1}:{SERIES[i][mode]}" for i in range(len(SERIES))]
            return ";".join(v)
        style = (
            f".{c}{{{css(0)}}}"
            f"@media (prefers-color-scheme:dark){{:root:where(:not([data-theme=\"light\"])) .{c}{{{css(1)}}}}}"
            f":root[data-theme=\"dark\"] .{c}{{{css(1)}}}"
            f".{c} text{{font-family:system-ui,-apple-system,'Segoe UI',sans-serif}}"
            f".{c} .ttl{{font-size:14.5px;font-weight:650;fill:var(--pri)}}"
            f".{c} .sub{{font-size:12px;fill:var(--sec)}}"
            f".{c} .lbl{{font-size:11.5px;fill:var(--sec)}}"
            f".{c} .tick{{font-size:11px;fill:var(--mut);font-variant-numeric:tabular-nums}}"
            f".{c} .val{{font-size:11px;font-weight:620;fill:var(--pri);font-variant-numeric:tabular-nums}}"
            f".{c} .note{{font-size:11px;fill:var(--mut)}}"
            f".{c} .key{{font-size:11.5px;font-weight:600;fill:var(--pri)}}"
            + DIAGRAM_CSS.replace("{c}", c)
            + PATHWAY_CSS.replace("{c}", c)
        )
        return (f'<svg class="{c}" viewBox="0 0 {fmt(self.width,0)} {fmt(self.height,0)}" '
                f'width="100%" role="img" aria-label="{esc(self.title)}" '
                f'xmlns="http://www.w3.org/2000/svg"><style>{style}</style>'
                f'<rect width="100%" height="100%" fill="var(--srf)"/>'
                + "".join(self.body) + "</svg>")


# --------------------------------------------------------------------------- #
# scales
# --------------------------------------------------------------------------- #
def linscale(d0, d1, r0, r1):
    span = (d1 - d0) or 1.0
    return lambda v: r0 + (v - d0) / span * (r1 - r0)


def nice_ticks(lo, hi, n=5):
    if hi <= lo:
        return [lo]
    raw = (hi - lo) / n
    mag = 10 ** math.floor(math.log10(raw))
    step = min((m * mag for m in (1, 2, 2.5, 5, 10)), key=lambda s: abs(s - raw))
    start = math.ceil(lo / step) * step
    out, v = [], start
    while v <= hi + step * 1e-9:
        out.append(round(v, 10))
        v += step
    return out


def xgrid(f: Fig, sx, ticks, nd=2, label=None):
    for t in ticks:
        x = sx(t)
        f.line(x, f.y0, x, f.y1, "var(--grd)", 1, cap="butt")
        f.text(x, f.y1 + 16, fmt(t, nd), "tick", "middle")
    f.line(f.x0, f.y1, f.x1, f.y1, "var(--axs)", 1, cap="butt")
    if label:
        f.text((f.x0 + f.x1) / 2, f.y1 + 34, label, "lbl", "middle")


def ygrid(f: Fig, sy, ticks, nd=2, label=None):
    for t in ticks:
        y = sy(t)
        f.line(f.x0, y, f.x1, y, "var(--grd)", 1, cap="butt")
        f.text(f.x0 - 8, y + 4, fmt(t, nd), "tick", "end")
    f.line(f.x0, f.y0, f.x0, f.y1, "var(--axs)", 1, cap="butt")
    if label:
        f.add(f'<text transform="translate(16,{fmt((f.y0+f.y1)/2,1)}) rotate(-90)" '
              f'class="lbl" text-anchor="middle">{esc(label)}</text>')


def legend(f: Fig, entries, x=None, y=22):
    """entries = [(label, colorvar)]. Always present for >= 2 series.

    Right-aligned to the plot edge. Place it on its own line -- overlapping the
    subtitle is the collision this signature makes easy to get wrong.
    """
    x = f.x1 if x is None else x
    cursor = x
    for lab, col in reversed(list(entries)):
        w = 7.2 * len(lab) + 22
        cursor -= w
        f.circle(cursor + 6, y - 4, 5, col)
        f.text(cursor + 17, y, lab, "key")
    return f


def wrap(text: str, width_px: float, px_per_char: float = 5.62) -> list[str]:
    """Greedy wrap by estimated advance width. 11px system sans ~5.6px/char."""
    limit = max(8, int(width_px / px_per_char))
    out, cur = [], ""
    for w in text.split():
        if cur and len(cur) + 1 + len(w) > limit:
            out.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        out.append(cur)
    return out


def footnote(f: Fig, text: str, y: float, x: float | None = None,
             width: float | None = None, cls: str = "note", lh: float = 14.5):
    """Wrapped caption. Returns the y after the last line so callers can stack."""
    x = 12 if x is None else x
    width = (f.width - x - 14) if width is None else width
    for i, line in enumerate(wrap(text, width)):
        f.text(x, y + i * lh, line, cls)
    return y + len(wrap(text, width)) * lh


def ordinal(n: float) -> str:
    """1 -> 1st, 3.1 -> 3rd, 11 -> 11th. Used for percentile labels."""
    i = int(round(n))
    if 10 <= i % 100 <= 20:
        suf = "th"
    else:
        suf = {1: "st", 2: "nd", 3: "rd"}.get(i % 10, "th")
    return f"{i}{suf}"


# --------------------------------------------------------------------------- #
# mechanism-diagram primitives
# --------------------------------------------------------------------------- #
def arrowhead_defs(f: Fig, ids=("arr", "arrA")) -> Fig:
    """Marker defs for arrowheads. Ids are fragment-internal, per the SVG rules."""
    parts = []
    for i, (idn, col) in enumerate(zip(ids, ("var(--sec)", "var(--s2)"))):
        parts.append(
            f'<marker id="{idn}" viewBox="0 0 10 10" refX="9" refY="5" '
            f'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
            f'<path d="M0,1 L9,5 L0,9 z" fill="{col}"/></marker>')
    return f.add("<defs>" + "".join(parts) + "</defs>")


def box(f: Fig, x, y, w, h, title, lines=(), accent=False, rx=8):
    """A labelled stage box. ``lines`` are short body lines (a formula, a phrase)."""
    stroke = "var(--s2)" if accent else "var(--axs)"
    f.rect(x, y, w, h, "var(--surface-2)" if False else "var(--bnd)", rx=rx,
           stroke=stroke, sw=1.5)
    f.text(x + w / 2, y + 21, title, "boxttl", "middle")
    for i, ln in enumerate(lines):
        f.text(x + w / 2, y + 42 + i * 15, ln, "boxln", "middle")
    return f


def arrow(f: Fig, x1, y1, x2, y2, label=None, accent=False, dash=None, above=True):
    """A labelled arrow. An unlabelled arrow only says 'related somehow'."""
    col = "var(--s2)" if accent else "var(--sec)"
    mk = "arrA" if accent else "arr"
    d = f' stroke-dasharray="{dash}"' if dash else ""
    f.add(f'<line x1="{fmt(x1,1)}" y1="{fmt(y1,1)}" x2="{fmt(x2,1)}" y2="{fmt(y2,1)}" '
          f'stroke="{col}" stroke-width="1.6" marker-end="url(#{mk})"{d}/>')
    if label:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        f.text(mx, my + (-7 if above else 15), label, "edge", "middle")
    return f


def elbow(f: Fig, pts, label=None, accent=False, dash=None):
    """A right-angled polyline arrow through ``pts`` = [(x, y), ...]."""
    col = "var(--s2)" if accent else "var(--sec)"
    mk = "arrA" if accent else "arr"
    d = f' stroke-dasharray="{dash}"' if dash else ""
    path = " ".join(f"{fmt(x,1)},{fmt(y,1)}" for x, y in pts)
    f.add(f'<polyline points="{path}" fill="none" stroke="{col}" stroke-width="1.6" '
          f'marker-end="url(#{mk})"{d}/>')
    if label:
        (ax, ay), (bx, by) = pts[0], pts[1]
        f.text((ax + bx) / 2, (ay + by) / 2 - 7, label, "edge", "middle")
    return f


def tag(f: Fig, x, y, text_, accent=True, anchor="start"):
    """A small annotation chip -- for the thing the argument hinges on."""
    cls = "tagA" if accent else "tagN"
    f.text(x, y, text_, cls, anchor)
    return f


#: Extra text classes the diagram primitives need, appended to Fig's scoped style.
#: Applied with ``.replace("{c}", cls)``, so braces here are SINGLE. They were
#: doubled once, for ``.format()``, which silently emitted invalid CSS -- every
#: class below then fell back to the SVG default 16px, which is how a monospace
#: formula ended up rendering as oversized sans.
DIAGRAM_CSS = (
    ".{c} .boxttl{font-size:12.5px;font-weight:650;fill:var(--pri)}"
    ".{c} .boxln{font-size:11px;fill:var(--sec);font-family:ui-monospace,"
    "'SF Mono',Menlo,Consolas,monospace}"
    ".{c} .edge{font-size:10.5px;fill:var(--mut)}"
    ".{c} .tagA{font-size:10.5px;font-weight:620;fill:var(--s2)}"
    ".{c} .tagN{font-size:10.5px;fill:var(--mut)}"
)


def edge_defs(f: Fig) -> Fig:
    """Markers for a pathway map: arrowhead = activates, bar = inhibits.

    The bar end is the standard biological convention for inhibition and reads
    without a legend lookup, unlike a dashed line.
    """
    parts = []
    for idn, col in (("act", "var(--mut)"), ("actA", "var(--s2)"), ("actB", "var(--s1)")):
        parts.append(
            f'<marker id="{idn}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5.5" '
            f'markerHeight="5.5" orient="auto-start-reverse">'
            f'<path d="M0,1.5 L9,5 L0,8.5 z" fill="{col}"/></marker>')
    for idn, col in (("inh", "var(--mut)"), ("inhA", "var(--s2)"), ("inhB", "var(--s1)")):
        parts.append(
            f'<marker id="{idn}" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="5.5" '
            f'markerHeight="5.5" orient="auto-start-reverse">'
            f'<rect x="6" y="0.5" width="2.2" height="9" fill="{col}"/></marker>')
    return f.add("<defs>" + "".join(parts) + "</defs>")


def curve(f: Fig, x1, y1, x2, y2, inhibit=False, accent=None, bow=0.28, opacity=None):
    """A quadratic Bezier edge with the right end-marker.

    ``accent`` is None (muted), "A" (orange) or "B" (blue) -- reserved for the two
    arms the model's mechanism turns on, so the rest of the graph stays recessive.
    """
    suffix = accent or ""
    mk = ("inh" if inhibit else "act") + suffix
    col = {"A": "var(--s2)", "B": "var(--s1)"}.get(accent, "var(--mut)")
    sw = 1.7 if accent else 1.0
    dx, dy = x2 - x1, y2 - y1
    # bow perpendicular to the chord, scaled by its length
    cx = (x1 + x2) / 2 - dy * bow
    cy = (y1 + y2) / 2 + dx * bow
    o = f' opacity="{opacity}"' if opacity is not None else ""
    return f.add(f'<path d="M{fmt(x1,1)},{fmt(y1,1)} Q{fmt(cx,1)},{fmt(cy,1)} '
                 f'{fmt(x2,1)},{fmt(y2,1)}" fill="none" stroke="{col}" '
                 f'stroke-width="{sw}" marker-end="url(#{mk})"{o}/>')


def node(f: Fig, x, y, w, h, label, badge=False, accent=None, bold=False):
    """A pathway node. ``badge`` marks a node the maturation rail feeds."""
    stroke = {"A": "var(--s2)", "B": "var(--s1)", "C": "var(--s3)"}.get(
        accent, "var(--axs)")
    f.rect(x, y, w, h, "var(--bnd)", rx=5, stroke=stroke, sw=1.6 if accent else 1)
    f.text(x + w / 2, y + h / 2 + 3.6, label,
           "nodeb" if (bold or accent) else "noden", "middle")
    if badge:
        f.circle(x + 4.5, y + 4.5, 3.0, "var(--s3)")
    return f


PATHWAY_CSS = (
    ".{c} .noden{font-size:10px;fill:var(--sec)}"
    ".{c} .nodeb{font-size:10px;font-weight:650;fill:var(--pri)}"
    ".{c} .colttl{font-size:11px;font-weight:650;letter-spacing:.05em;"
    "text-transform:uppercase;fill:var(--mut)}"
)
