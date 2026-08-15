#!/usr/bin/env python3
"""Assemble the lab-meeting deck: inline the repo's generated SVGs verbatim,
plus five hand-authored diagrams. Numbers here come from model_results.json /
results.json / cmfate_model.toml -- never from RESULTS.md sections 6-8, which
are stale against the regenerated JSON."""
import pathlib

REPO = pathlib.Path("/home/justin/Documents/GitHub/e2f-heart-scrna/model")
HERE = pathlib.Path(__file__).parent
FIGS = REPO / "figures"


# ---------------------------------------------------------------- tree
def tree(sfx):
    """Three sequential gates, four outcomes. The shape of the fate layer."""
    gates = [("S-phase entry", 274), ("Mitotic entry", 502), ("Abscission", 730)]
    outs = [
        (274, "Quiescent", "var(--mut)", "never enters S", ""),
        (502, "Polyploidization", "var(--s2)", "G2 arrest — APC/C never fires", "ploidy ↑, nuclei flat"),
        (730, "Binucleation", "var(--s1)", "mitosis completes, furrow fails", "nuclei ↑, ploidy flat"),
    ]
    p = []
    p.append('<defs><marker id="ar%s" viewBox="0 0 10 10" refX="9" refY="5" '
             'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
             '<path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>' % sfx)

    # start
    p.append('<rect x="6" y="34" width="142" height="58" rx="3" fill="none" stroke="var(--axs)"/>')
    p.append('<text x="77" y="68" text-anchor="middle" font-size="13" font-weight="600" '
             'fill="currentColor">Cardiomyocyte</text>')

    # gate boxes + horizontal chain
    for name, cx in gates:
        p.append('<rect x="%d" y="34" width="156" height="58" rx="3" fill="none" stroke="var(--axs)"/>'
                 % (cx - 78))
        p.append('<text x="%d" y="68" text-anchor="middle" font-size="13" font-weight="600" '
                 'fill="currentColor">%s</text>' % (cx, name))

    for x1, x2 in [(148, 190), (352, 418), (580, 646), (808, 856)]:
        p.append('<line x1="%d" y1="63" x2="%d" y2="63" stroke="currentColor" stroke-width="1.4" '
                 'marker-end="url(#ar%s)"/>' % (x1, x2, sfx))
    for xm in [169, 385, 613, 832]:
        p.append('<text x="%d" y="53" text-anchor="middle" font-size="10" fill="var(--mut)">yes</text>' % xm)

    # division terminal
    p.append('<rect x="856" y="34" width="140" height="58" rx="3" fill="none" '
             'stroke="var(--s3)" stroke-width="1.6"/>')
    p.append('<text x="926" y="66" text-anchor="middle" font-size="13.5" font-weight="700" '
             'fill="var(--s3)">Division</text>')
    p.append('<text x="926" y="112" text-anchor="middle" font-size="10.5" fill="var(--mut)">'
             'AurKB⁺ midbody</text>')
    p.append('<text x="926" y="127" text-anchor="middle" font-size="10.5" fill="var(--mut)">'
             'cell number rises</text>')

    # no-branches
    for cx, name, col, s1, s2 in outs:
        p.append('<line x1="%d" y1="92" x2="%d" y2="208" stroke="currentColor" stroke-width="1.4" '
                 'marker-end="url(#ar%s)"/>' % (cx, cx, sfx))
        p.append('<text x="%d" y="156" font-size="10" fill="var(--mut)">no</text>' % (cx + 10))
        w = 168
        p.append('<rect x="%d" y="214" width="%d" height="56" rx="3" fill="none" stroke="%s" '
                 'stroke-width="1.6"/>' % (cx - w // 2, w, col))
        p.append('<text x="%d" y="248" text-anchor="middle" font-size="13.5" font-weight="700" '
                 'fill="%s">%s</text>' % (cx, col, name))
        p.append('<text x="%d" y="290" text-anchor="middle" font-size="10.5" fill="var(--mut)">%s</text>'
                 % (cx, s1))
        if s2:
            p.append('<text x="%d" y="305" text-anchor="middle" font-size="10.5" '
                     'fill="var(--mut)">%s</text>' % (cx, s2))

    return ('<svg viewBox="0 0 1010 320" width="100%%" role="img" aria-label="Three sequential gates '
            '— S-phase entry, mitotic entry and abscission — produce four cardiomyocyte '
            'outcomes: quiescence, polyploidization, binucleation and division." '
            'xmlns="http://www.w3.org/2000/svg">%s</svg>') % "".join(p)


# ---------------------------------------------------------------- KO bars
def kobars():
    """CM %cycling, KO vs WT, from ko_export/tbl_cellcycle_fraction.csv."""
    rows = [("P0", 15.5, 16.5, "12,263", "13,019", 99, 239),
            ("P7", 31.6, 25.6, "10,597", "6,537", 0, 0)]
    top, bot, vmax = 44, 272, 35.0

    def y(v):
        return bot - (v / vmax) * (bot - top)

    p = []
    for g in (0, 10, 20, 30):
        p.append('<line x1="54" y1="%.1f" x2="400" y2="%.1f" stroke="var(--grd)" stroke-width="1"/>'
                 % (y(g), y(g)))
        p.append('<text x="46" y="%.1f" text-anchor="end" font-size="10.5" fill="var(--mut)" '
                 'font-family="ui-monospace,monospace">%d</text>' % (y(g) + 3.5, g))
    p.append('<text x="46" y="%.1f" text-anchor="end" font-size="10.5" fill="var(--mut)" '
             'font-family="ui-monospace,monospace">35</text>' % (y(35) + 3.5))

    groups = [("P0", 15.5, 16.5, "12,263", "13,019", 99), ("P7", 31.6, 25.6, "10,597", "6,537", 239)]
    for label, ko, wt, nko, nwt, x0 in groups:
        for v, x, col, fill in [(ko, x0, "var(--s1)", "var(--s1)"), (wt, x0 + 56, "var(--axs)", "none")]:
            p.append('<rect x="%d" y="%.1f" width="46" height="%.1f" fill="%s" stroke="%s" '
                     'stroke-width="1.4"/>' % (x, y(v), bot - y(v), fill, col))
            p.append('<text x="%d" y="%.1f" text-anchor="middle" font-size="11" font-weight="650" '
                     'fill="currentColor" font-family="ui-monospace,monospace">%.1f</text>'
                     % (x + 23, y(v) - 7, v))
        cx = x0 + 51
        p.append('<text x="%d" y="292" text-anchor="middle" font-size="12.5" font-weight="650" '
                 'fill="currentColor">%s</text>' % (cx, label))
        p.append('<text x="%d" y="309" text-anchor="middle" font-size="9.5" fill="var(--mut)" '
                 'font-family="ui-monospace,monospace">n %s / %s</text>' % (cx, nko, nwt))

    p.append('<line x1="54" y1="272" x2="400" y2="272" stroke="var(--axs)" stroke-width="1.2"/>')

    # the delta that matters -- bracketed clear of both bars and their value labels
    gx = 352
    p.append('<line x1="285" y1="%.1f" x2="%d" y2="%.1f" stroke="var(--grd)" stroke-width="1" '
             'stroke-dasharray="3 2"/>' % (y(31.6), gx, y(31.6)))
    p.append('<line x1="341" y1="%.1f" x2="%d" y2="%.1f" stroke="var(--grd)" stroke-width="1" '
             'stroke-dasharray="3 2"/>' % (y(25.6), gx, y(25.6)))
    p.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="var(--s1)" stroke-width="1.3"/>'
             % (gx, y(31.6), gx, y(25.6)))
    for yy in (y(31.6), y(25.6)):
        p.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="var(--s1)" stroke-width="1.3"/>'
                 % (gx - 4, yy, gx + 4, yy))
    p.append('<text x="%d" y="%.1f" font-size="11" font-weight="650" fill="var(--s1)" '
             'font-family="ui-monospace,monospace">+6.0</text>' % (gx + 8, (y(31.6) + y(25.6)) / 2 + 1))
    p.append('<text x="%d" y="%.1f" font-size="9.5" fill="var(--mut)">pp</text>'
             % (gx + 8, (y(31.6) + y(25.6)) / 2 + 13))

    # legend
    p.append('<rect x="54" y="14" width="12" height="12" fill="var(--s1)"/>')
    p.append('<text x="72" y="24" font-size="11" fill="var(--sec)">KO</text>')
    p.append('<rect x="108" y="14" width="12" height="12" fill="none" stroke="var(--axs)" stroke-width="1.4"/>')
    p.append('<text x="126" y="24" font-size="11" fill="var(--sec)">WT</text>')
    p.append('<text x="400" y="24" text-anchor="end" font-size="10.5" fill="var(--mut)">'
             '% cycling (S + G2M)</text>')

    return ('<svg viewBox="0 0 420 320" width="100%%" role="img" aria-label="Cardiomyocyte percent '
            'cycling, knockout versus wild type: 15.5 versus 16.5 at P0, and 31.6 versus 25.6 at P7, '
            'a 6 percentage point rise in the knockout at P7." '
            'xmlns="http://www.w3.org/2000/svg">%s</svg>') % "".join(p)


# ---------------------------------------------------------------- pipeline
def pipeline():
    lanes = [
        dict(y=36, col="var(--s1)", src="OUR LAB",
             lines=["10x · P0 + P7 · KO vs WT", "30,030 cells / 21,598 CM"],
             warn="⚠ 29 of 63 model genes in panel",
             proc=["shiny_app browser", "→ ko_export/*.csv"],
             give="MOTIVATION", sub="+ validation target"),
        dict(y=166, col="var(--s2)", src="BANIOL 2021",
             lines=["Smart-seq2, FUCCI-sorted", "285 cells · ENA PRJEB47622"],
             warn="⚠ P7 sort 4.5–5.2× enriched",
             proc=["cmcycle.baniol", "Ect2(M), E2f6(M), E2F roles"],
             give="TOPOLOGY", sub="which edges exist"),
        dict(y=296, col="var(--s3)", src="PUBLISHED MEASUREMENTS",
             lines=["cmcycle_targets.csv · 42 rows", "28 Murganti · 13 Baniol · 1 Hashimoto"],
             warn="✓ 4 of 42 fitted; 38 held out",
             proc=["cmcycle.preflight — 6 checks", "closed-form gate inversion"],
             give="PARAMETERS", sub="rates and fractions"),
    ]
    p = ['<defs><marker id="arp" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" '
         'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>']

    for x, t in [(8, "SOURCE"), (322, "PROCESSING"), (662, "CONTRIBUTES")]:
        p.append('<text x="%d" y="20" font-size="9.5" font-weight="650" letter-spacing="1.4" '
                 'fill="var(--mut)">%s</text>' % (x, t))

    for L in lanes:
        y, c = L["y"], L["col"]
        p.append('<rect x="8" y="%d" width="292" height="108" rx="3" fill="none" stroke="var(--grd)"/>' % y)
        p.append('<rect x="8" y="%d" width="3" height="108" fill="%s"/>' % (y, c))
        p.append('<text x="24" y="%d" font-size="12" font-weight="700" fill="%s">%s</text>'
                 % (y + 24, c, L["src"]))
        for k, ln in enumerate(L["lines"]):
            p.append('<text x="24" y="%d" font-size="10.5" fill="var(--sec)">%s</text>' % (y + 46 + k * 16, ln))
        wc = "var(--neg)" if L["warn"].startswith("⚠") else "var(--s3)"
        p.append('<text x="24" y="%d" font-size="10" fill="%s">%s</text>' % (y + 94, wc, L["warn"]))

        p.append('<line x1="300" y1="%d" x2="318" y2="%d" stroke="currentColor" stroke-width="1.3" '
                 'marker-end="url(#arp)"/>' % (y + 54, y + 54))
        p.append('<rect x="322" y="%d" width="272" height="108" rx="3" fill="none" stroke="var(--grd)"/>' % y)
        for k, ln in enumerate(L["proc"]):
            p.append('<text x="338" y="%d" font-size="10.5" fill="var(--sec)" '
                     'font-family="ui-monospace,monospace">%s</text>' % (y + 48 + k * 18, ln))

        p.append('<line x1="594" y1="%d" x2="654" y2="%d" stroke="currentColor" stroke-width="1.3" '
                 'marker-end="url(#arp)"/>' % (y + 54, y + 54))
        p.append('<text x="662" y="%d" font-size="12.5" font-weight="700" letter-spacing="0.6" '
                 'fill="%s">%s</text>' % (y + 50, c, L["give"]))
        p.append('<text x="662" y="%d" font-size="10" fill="var(--mut)">%s</text>' % (y + 66, L["sub"]))

    # the firewall
    p.append('<line x1="626" y1="10" x2="626" y2="424" stroke="var(--axs)" stroke-width="1.2" '
             'stroke-dasharray="4 4"/>')

    # converge
    for L, ty in zip(lanes, (206, 231, 256)):
        sy = L["y"] + 54
        p.append('<path d="M834,%d C862,%d 862,%d 884,%d" stroke="%s" stroke-width="1.5" fill="none" '
                 'marker-end="url(#arp)"/>' % (sy, sy, ty, ty, L["col"]))

    p.append('<rect x="888" y="180" width="200" height="104" rx="3" fill="none" stroke="currentColor" '
             'stroke-width="1.8"/>')
    p.append('<text x="988" y="214" text-anchor="middle" font-size="16" font-weight="700" '
             'fill="currentColor" font-family="ui-monospace,monospace">cmfate</text>')
    p.append('<text x="988" y="236" text-anchor="middle" font-size="10.5" fill="var(--sec)">'
             '55 nodes · 77 reactions</text>')
    p.append('<text x="988" y="252" text-anchor="middle" font-size="10.5" fill="var(--sec)">'
             '3 gates → 4 fates</text>')
    p.append('<text x="988" y="270" text-anchor="middle" font-size="10" fill="var(--mut)">'
             'stdlib only</text>')

    p.append('<line x1="988" y1="284" x2="988" y2="326" stroke="currentColor" stroke-width="1.5" '
             'marker-end="url(#arp)"/>')
    p.append('<rect x="888" y="330" width="200" height="58" rx="3" fill="none" stroke="var(--grd)"/>')
    p.append('<text x="988" y="352" text-anchor="middle" font-size="10.5" fill="var(--sec)">'
             'clonidine triad · in-silico KO</text>')
    p.append('<text x="988" y="370" text-anchor="middle" font-size="10.5" fill="var(--sec)">'
             'regeneration screen</text>')

    return ('<svg viewBox="0 0 1100 420" width="100%%" role="img" aria-label="Three data lanes converge '
            'on one model: our knockout single-cell data contributes motivation, the Baniol 285-cell '
            're-analysis contributes topology, and 42 curated published measurements contribute '
            'parameters. A dashed firewall separates sources from what each contributes." '
            'xmlns="http://www.w3.org/2000/svg">%s</svg>') % "".join(p)


# ---------------------------------------------------------------- two arms
def arms():
    p = ['<defs><marker id="arm" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5.5" '
         'markerHeight="5.5" orient="auto-start-reverse">'
         '<path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>']

    orange = ["E2Fact", "Ect2", "RhoA", "Midbody", "Abscission"]
    blue = ["ROS", "DDR", "Ccng1/Pkmyt1", "MitoticEntry"]

    for i, n in enumerate(orange):
        x = 12 + i * 88
        p.append('<rect x="%d" y="60" width="76" height="28" rx="14" fill="none" stroke="var(--s2)" '
                 'stroke-width="1.5"/>' % x)
        p.append('<text x="%d" y="78" text-anchor="middle" font-size="9.5" font-weight="600" '
                 'fill="var(--s2)">%s</text>' % (x + 38, n))
        if i < len(orange) - 1:
            p.append('<line x1="%d" y1="74" x2="%d" y2="74" stroke="var(--s2)" stroke-width="1.3" '
                     'marker-end="url(#arm)"/>' % (x + 76, x + 86))

    for i, n in enumerate(blue):
        x = 12 + i * 110
        p.append('<rect x="%d" y="168" width="98" height="28" rx="14" fill="none" stroke="var(--s1)" '
                 'stroke-width="1.5"/>' % x)
        p.append('<text x="%d" y="186" text-anchor="middle" font-size="9.5" font-weight="600" '
                 'fill="var(--s1)">%s</text>' % (x + 49, n))
        if i < len(blue) - 1:
            p.append('<line x1="%d" y1="182" x2="%d" y2="182" stroke="var(--s1)" stroke-width="1.3" '
                     'marker-end="url(#arm)"/>' % (x + 98, x + 108))

    p.append('<text x="12" y="34" font-size="10.5" font-weight="650" fill="var(--s2)">'
             'ABSCISSION ARM — closed by maturation</text>')
    p.append('<text x="12" y="228" font-size="10.5" font-weight="650" fill="var(--s1)">'
             'MITOTIC-ENTRY BRAKE — closed by culture</text>')
    p.append('<text x="440" y="132" text-anchor="end" font-size="10" fill="var(--mut)" '
             'font-style="italic">which shuts first is the whole 2×2</text>')

    return ('<svg viewBox="0 0 448 244" width="100%%" role="img" aria-label="Two arms: an abscission arm '
            'from E2F activity through Ect2, RhoA and midbody to abscission, closed by maturation; and a '
            'mitotic-entry brake from ROS through DDR and Ccng1/Pkmyt1, closed by culture." '
            'xmlns="http://www.w3.org/2000/svg">%s</svg>') % "".join(p)


# ---------------------------------------------------------------- unfitted
def unfitted():
    p = []
    # ---- panel 1: the level
    p.append('<text x="0" y="20" font-size="11" font-weight="650" letter-spacing="1.2" '
             'fill="var(--mut)">THE LEVEL — NOTHING FITTED</text>')
    x0, x1, vmax = 46, 430, 35.0

    def X(v):
        return x0 + (v / vmax) * (x1 - x0)

    p.append('<line x1="%d" y1="150" x2="%d" y2="150" stroke="var(--axs)" stroke-width="1.2"/>' % (x0, x1))
    for g in (0, 10, 20, 30):
        p.append('<line x1="%.1f" y1="150" x2="%.1f" y2="156" stroke="var(--axs)"/>' % (X(g), X(g)))
        p.append('<text x="%.1f" y="170" text-anchor="middle" font-size="10" fill="var(--mut)" '
                 'font-family="ui-monospace,monospace">%d%%</text>' % (X(g), g))

    # struck-through naive comparator
    p.append('<circle cx="%.1f" cy="150" r="4.5" fill="none" stroke="var(--mut)" stroke-width="1.4"/>'
             % X(32.5))
    p.append('<text x="%.1f" y="104" text-anchor="middle" font-size="11" fill="var(--mut)" '
             'font-family="ui-monospace,monospace">32.5%%</text>' % X(32.5))
    p.append('<line x1="%.1f" y1="100" x2="%.1f" y2="100" stroke="var(--mut)" stroke-width="1.3"/>'
             % (X(32.5) - 20, X(32.5) + 20))
    p.append('<text x="%.1f" y="122" text-anchor="middle" font-size="9.5" fill="var(--mut)">'
             'naive 1 − mKO2⁺</text>' % X(32.5))
    p.append('<line x1="%.1f" y1="128" x2="%.1f" y2="142" stroke="var(--mut)" stroke-width="1" '
             'stroke-dasharray="2 2"/>' % (X(32.5), X(32.5)))

    # the pair that nearly coincides
    p.append('<line x1="%.1f" y1="150" x2="%.1f" y2="150" stroke="var(--s3)" stroke-width="4"/>'
             % (X(17.17), X(17.7)))
    p.append('<circle cx="%.1f" cy="150" r="5.5" fill="var(--s3)"/>' % X(17.7))
    p.append('<circle cx="%.1f" cy="150" r="5.5" fill="var(--srf)" stroke="var(--s3)" stroke-width="2"/>'
             % X(17.17))
    p.append('<line x1="%.1f" y1="145" x2="%.1f" y2="72" stroke="var(--grd)" stroke-width="1"/>' % (X(17.7), X(17.7)))
    p.append('<text x="%.1f" y="64" text-anchor="middle" font-size="12.5" font-weight="700" '
             'fill="var(--s3)" font-family="ui-monospace,monospace">1.03×</text>' % X(17.7))
    p.append('<text x="%.1f" y="196" text-anchor="middle" font-size="10.5" fill="currentColor" '
             'font-family="ui-monospace,monospace">17.70%%</text>' % (X(17.7) - 30))
    p.append('<text x="%.1f" y="210" text-anchor="middle" font-size="9.5" fill="var(--mut)">observed</text>'
             % (X(17.7) - 30))
    p.append('<text x="%.1f" y="196" text-anchor="middle" font-size="10.5" fill="currentColor" '
             'font-family="ui-monospace,monospace">17.17%%</text>' % (X(17.17) + 46))
    p.append('<text x="%.1f" y="210" text-anchor="middle" font-size="9.5" fill="var(--mut)">model ceiling</text>'
             % (X(17.17) + 46))
    p.append('<text x="0" y="240" font-size="10" fill="var(--mut)" font-family="ui-monospace,monospace">'
             'corrected with Baniol’s own Ki-67 co-staining (Suppl 1G)</text>')

    # ---- panel 2: the slope
    o = 500
    p.append('<text x="%d" y="20" font-size="11" font-weight="650" letter-spacing="1.2" '
             'fill="var(--mut)">THE SLOPE — THE SAME CORRECTION EXPOSED IT</text>' % o)
    bx, bw, smax = o + 96, 268, 14.0
    for lab, v, col, y in [("observed", 3.23, "var(--s3)", 74), ("modelled", 13.1, "var(--neg)", 122)]:
        w = (v / smax) * bw
        p.append('<rect x="%d" y="%d" width="%.1f" height="30" fill="%s"/>' % (bx, y, w, col))
        p.append('<text x="%d" y="%d" text-anchor="end" font-size="11" fill="var(--sec)">%s</text>'
                 % (bx - 12, y + 20, lab))
        p.append('<text x="%.1f" y="%d" font-size="12" font-weight="700" fill="%s" '
                 'font-family="ui-monospace,monospace">%.1f×</text>' % (bx + w + 10, y + 20, col, v))
    p.append('<line x1="%d" y1="66" x2="%d" y2="160" stroke="var(--axs)" stroke-width="1.2"/>' % (bx, bx))
    p.append('<text x="%d" y="184" font-size="10.5" fill="var(--sec)">P0→P7 fall in cycling '
             'cardiomyocytes</text>' % (o + 96))
    p.append('<text x="%d" y="206" font-size="10.5" fill="var(--neg)">≈ 4× too steep — kept as a '
             'test that fails</text>' % (o + 96))
    p.append('<text x="%d" y="222" font-size="10.5" fill="var(--neg)">by design, not tuned away</text>' % (o + 96))
    p.append('<text x="%d" y="240" font-size="10" fill="var(--mut)" font-family="ui-monospace,monospace">'
             'root cause: CycE travels ~36,000×</text>' % o)

    return ('<svg viewBox="0 0 900 256" width="100%%" role="img" aria-label="Left: observed 17.70 percent '
            'against a model ceiling of 17.17 percent, a 1.03-fold match with nothing fitted, versus a '
            'struck-through naive comparator of 32.5 percent. Right: the maturation slope, 3.2-fold '
            'observed against 13.1-fold modelled." '
            'xmlns="http://www.w3.org/2000/svg">%s</svg>') % "".join(p)


# ---------------------------------------------------------------- assemble
html = (HERE / "deck_template.html").read_text()

# repo-generated SVGs, inlined verbatim so a slide can never disagree with its figure
for tok, name in [("FIG1", "fig1_ect2_specificity"), ("FIG2", "fig2_module_cancellation"),
                  ("FIG5", "fig5_e2f_family"), ("FIG6", "fig6_preflight_checks"),
                  ("FIG7", "fig7_sort_enrichment"), ("FIG9", "fig9_clonidine_triad"),
                  ("FIG10", "fig10_ko_and_screen"), ("FIG11", "fig11_engine")]:
    svg = (FIGS / (name + ".svg")).read_text().strip()
    assert svg.startswith("<svg"), name
    html = html.replace("<!--%s-->" % tok, svg)

# hand-authored diagrams
html = html.replace("<!--FIG_TREE-->", tree("w"), 1)   # slide 1 wash
html = html.replace("<!--FIG_TREE-->", tree("a"), 1)   # slide 2 figure
html = html.replace("<!--FIG_KOBARS-->", kobars())
html = html.replace("<!--FIG_PIPE-->", pipeline())
html = html.replace("<!--FIG_ARMS-->", arms())
html = html.replace("<!--FIG_UNFIT-->", unfitted())

assert "<!--" not in html.split("<script>")[0].replace("<!-- ", "@@"), "unsubstituted token remains"
left = [t for t in ["FIG_TREE", "FIG_KOBARS", "FIG_PIPE", "FIG_ARMS", "FIG_UNFIT",
                    "FIG1", "FIG2", "FIG5", "FIG6", "FIG7", "FIG9", "FIG10", "FIG11"]
        if "<!--%s-->" % t in html]
assert not left, "unsubstituted: %s" % left

out = REPO / "deck.html"   # generated; gitignored like report.html
out.write_text(html)
print("wrote %s  (%.0f KB)" % (out, len(html) / 1024))
print("slides:", html.count('<section class="slide'))
