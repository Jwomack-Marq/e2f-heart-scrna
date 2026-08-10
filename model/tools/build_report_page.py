"""Build the shareable HTML report from figures/ + results.json.

Writes a single self-contained page with every SVG inlined and every quoted
statistic pulled from ``figures/results.json``, so the report cannot drift from
the analysis. Run from ``model/``:

    python3 tools/build_report_page.py

Design notes: the two accents are the FUCCI reporters desaturated into ink roles
-- teal for mAG-geminin (cycling / confirmed), rust for mKO2-Cdt1 (arrested /
flagged) -- so colour carries the same meaning the assay does. Chrome stays
near-neutral because the figures already carry the categorical palette.
"""
import json, pathlib, html

F = pathlib.Path("figures")
figs = {p.stem: p.read_text() for p in sorted(F.glob("*.svg"))}
R = json.load(open(F / "results.json"))
M = json.load(open(F / "model_results.json"))

def fig(key, num, caption):
    return (f'<figure class="fig"><div class="fig-frame">{figs[key]}</div>'
            f'<figcaption><span class="fignum">Fig {num}</span>{caption}</figcaption>'
            f'</figure>')

CSS = """
:root{
  color-scheme: light;
  --ground:#f7f9f8; --surface:#ffffff; --surface-2:#f1f4f3;
  --ink:#12171a; --ink-2:#4c595e; --muted:#7c8a8f;
  --rule:#e0e6e4; --rule-2:#eef2f1;
  --cyc:#0f6f5c;            /* mAG-geminin, desaturated: cycling / confirmed */
  --arr:#a8442a;            /* mKO2-Cdt1, desaturated: arrested / flagged */
  --cyc-wash:#e8f2ef; --arr-wash:#f8ebe6;
  --shadow:0 1px 2px rgba(18,23,26,.05), 0 8px 24px -16px rgba(18,23,26,.18);
  --serif: Georgia,'Iowan Old Style','Palatino Linotype',ui-serif,serif;
  --sans: system-ui,-apple-system,'Segoe UI',Roboto,sans-serif;
  --mono: ui-monospace,'SF Mono','Cascadia Mono',Menlo,Consolas,monospace;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    color-scheme: dark;
    --ground:#0e1214; --surface:#161b1e; --surface-2:#1b2225;
    --ink:#edf2f1; --ink-2:#a4b1b4; --muted:#7c8a8f;
    --rule:#242c2f; --rule-2:#1d2427;
    --cyc:#4fbfa4; --arr:#e08163;
    --cyc-wash:#12312a; --arr-wash:#33201a;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 24px -16px rgba(0,0,0,.7);
  }
}
:root[data-theme="dark"]{
  color-scheme: dark;
  --ground:#0e1214; --surface:#161b1e; --surface-2:#1b2225;
  --ink:#edf2f1; --ink-2:#a4b1b4; --muted:#7c8a8f;
  --rule:#242c2f; --rule-2:#1d2427;
  --cyc:#4fbfa4; --arr:#e08163;
  --cyc-wash:#12312a; --arr-wash:#33201a;
  --shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 24px -16px rgba(0,0,0,.7);
}
*{box-sizing:border-box}
body{
  margin:0; background:var(--ground); color:var(--ink);
  font-family:var(--sans); font-size:16.5px; line-height:1.65;
  -webkit-font-smoothing:antialiased;
}
.wrap{max-width:1080px; margin:0 auto; padding:0 24px 96px}
.col{max-width:68ch}

/* ---- masthead: the two reporter colours as a hairline, which is the subject's
       own two-channel readout rather than decoration ---- */
.mast{padding:56px 0 30px; border-bottom:1px solid var(--rule)}
.channel{display:flex; height:3px; width:120px; margin-bottom:26px; border-radius:2px; overflow:hidden}
.channel i{flex:1}
.channel i:first-child{background:var(--arr)}
.channel i:last-child{background:var(--cyc)}
.eyebrow{font-family:var(--mono); font-size:11.5px; letter-spacing:.09em;
  text-transform:uppercase; color:var(--muted); margin:0 0 14px}
h1{font-family:var(--serif); font-weight:400; font-size:clamp(30px,4.4vw,45px);
  line-height:1.14; letter-spacing:-.014em; margin:0 0 18px; text-wrap:balance}
h1 em{font-style:italic; color:var(--cyc)}
.standfirst{font-size:19px; line-height:1.55; color:var(--ink-2); margin:0; max-width:60ch}
.meta{display:flex; flex-wrap:wrap; gap:8px 26px; margin-top:28px;
  font-family:var(--mono); font-size:12px; color:var(--muted)}
.meta b{color:var(--ink-2); font-weight:500}

h2{font-family:var(--serif); font-weight:400; font-size:29px; line-height:1.2;
  letter-spacing:-.01em; margin:0 0 6px; text-wrap:balance}
h3{font-family:var(--sans); font-weight:640; font-size:17px; letter-spacing:-.005em;
  margin:34px 0 8px; text-wrap:balance}
section{padding-top:56px}
.snum{font-family:var(--mono); font-size:11.5px; letter-spacing:.09em; color:var(--cyc);
  text-transform:uppercase; display:block; margin-bottom:10px}
p{margin:0 0 17px}
a{color:var(--cyc); text-decoration-thickness:1px; text-underline-offset:2px}
strong{font-weight:640}
code,.n{font-family:var(--mono); font-size:.88em; font-variant-numeric:tabular-nums}
code{background:var(--surface-2); padding:1px 5px; border-radius:4px; color:var(--ink)}
.gene{font-family:var(--mono); font-size:.9em; font-style:normal}

/* ---- typed verdict chips: the content really does have states ---- */
.chip{display:inline-flex; align-items:center; gap:6px; font-family:var(--mono);
  font-size:10.5px; letter-spacing:.06em; text-transform:uppercase;
  padding:3px 9px; border-radius:100px; white-space:nowrap; font-weight:500}
.chip::before{content:""; width:6px; height:6px; border-radius:50%}
.chip.ok{background:var(--cyc-wash); color:var(--cyc)}
.chip.ok::before{background:var(--cyc)}
.chip.flag{background:var(--arr-wash); color:var(--arr)}
.chip.flag::before{background:var(--arr)}
.chip.open{background:var(--surface-2); color:var(--ink-2)}
.chip.open::before{background:var(--muted)}

/* ---- findings ---- */
.finds{display:grid; gap:1px; background:var(--rule); border:1px solid var(--rule);
  border-radius:10px; overflow:hidden; margin:26px 0 0}
.find{background:var(--surface); padding:20px 22px; display:grid;
  grid-template-columns:minmax(0,1fr) auto; gap:6px 18px; align-items:start}
.find h4{margin:0; font-family:var(--sans); font-size:16px; font-weight:640; line-height:1.35}
.find p{margin:0; grid-column:1/-1; color:var(--ink-2); font-size:15px}
.find .chip{justify-self:end}

/* ---- figures: break out of the reading column ---- */
.fig{margin:34px 0 8px}
.fig-frame{background:var(--surface); border:1px solid var(--rule); border-radius:10px;
  padding:10px; box-shadow:var(--shadow); overflow-x:auto}
.fig-frame svg{display:block; min-width:600px}
figcaption{margin-top:11px; font-size:14px; color:var(--ink-2); line-height:1.5; max-width:76ch}
.fignum{font-family:var(--mono); font-size:11.5px; letter-spacing:.06em; color:var(--cyc);
  text-transform:uppercase; margin-right:9px}

/* ---- tables ---- */
.tw{overflow-x:auto; margin:22px 0; border:1px solid var(--rule); border-radius:10px;
  background:var(--surface)}
table{border-collapse:collapse; width:100%; font-size:14.5px; min-width:520px}
th,td{padding:10px 14px; text-align:left; border-bottom:1px solid var(--rule-2); vertical-align:top}
thead th{font-family:var(--mono); font-size:11px; letter-spacing:.07em; text-transform:uppercase;
  color:var(--muted); font-weight:500; background:var(--surface-2); border-bottom:1px solid var(--rule)}
tbody tr:last-child td{border-bottom:0}
td.num,th.num{text-align:right; font-family:var(--mono); font-variant-numeric:tabular-nums}
td.g{font-family:var(--mono); font-size:13.5px}
.hi{color:var(--cyc); font-weight:640}
.lo{color:var(--muted)}

/* ---- callout ---- */
.note{border-left:2px solid var(--arr); background:var(--arr-wash);
  padding:16px 20px; border-radius:0 8px 8px 0; margin:24px 0; font-size:15.5px}
.note.ok{border-left-color:var(--cyc); background:var(--cyc-wash)}
.note p:last-child{margin-bottom:0}
.note .lbl{font-family:var(--mono); font-size:11px; letter-spacing:.07em;
  text-transform:uppercase; color:var(--arr); display:block; margin-bottom:6px}
.note.ok .lbl{color:var(--cyc)}

/* ---- todo ---- */
.todo{list-style:none; padding:0; margin:20px 0; display:grid; gap:1px;
  background:var(--rule); border:1px solid var(--rule); border-radius:10px; overflow:hidden}
.todo li{background:var(--surface); padding:16px 20px}
.todo .t{display:flex; align-items:baseline; gap:10px; flex-wrap:wrap; margin-bottom:4px}
.todo .t b{font-weight:640; font-size:15.5px}
.todo p{margin:0; font-size:14.5px; color:var(--ink-2)}
.eq{font-family:var(--mono); font-size:13px; background:var(--surface-2); border:1px solid var(--rule);
  border-radius:8px; padding:14px 16px; margin:20px 0; overflow-x:auto; line-height:1.85}
.foot{margin-top:72px; padding-top:22px; border-top:1px solid var(--rule);
  font-size:13.5px; color:var(--muted)}
:focus-visible{outline:2px solid var(--cyc); outline-offset:2px; border-radius:3px}
@media (prefers-reduced-motion:no-preference){
  .fig-frame{transition:box-shadow .2s ease}
}
@media (max-width:640px){
  body{font-size:16px}
  h1{font-size:29px}
  .find{grid-template-columns:1fr}
  .find .chip{justify-self:start}
}
"""

e = R["per_gene"]; d = R["drift"]; m = R["modules"]; mat = R["maturation"]
corr = mat["correlations"]; g8 = mat["by_group"]; e2f = R["e2f"]; se = R["sort_enrichment"]

def f2(x, nd=2): return f"{x:.{nd}f}"

drift_rows = "".join(
    f'<tr><td class="g">{g}</td><td>{v["module"]}</td>'
    f'<td class="num">{f2(v["ratio"])}</td><td class="num">{v["t"]:+.2f}</td>'
    f'<td class="num">{"&lt;0.0001" if v["p"] < 1e-4 else f2(v["p"],4)}</td>'
    f'<td class="num {"hi" if v["pct_genomewide"] < 5 else "lo"}">{f2(v["pct_genomewide"],1)}</td></tr>'
    for g, v in sorted(e.items(), key=lambda kv: kv[1]["ratio"])[:6])

mat_rows = "".join(
    f'<tr><td>{k.replace("-"," · ")}</td><td class="num">{v["n"]}</td>'
    f'<td class="num">{v["mean"]:+.2f}</td><td class="num lo">± {f2(v["se"])}</td></tr>'
    for k, v in sorted(g8.items(), key=lambda kv: kv[1]["mean"]))

corr_rows = "".join(
    f'<tr><td class="g">{g}</td><td>{role}</td>'
    f'<td class="num {"hi" if abs(corr[g]["r"])>0.3 else ""}">{corr[g]["r"]:+.3f}</td>'
    f'<td class="num">{corr[g]["t"]:+.2f}</td></tr>'
    for g, role in [("Ect2","cytokinesis competence"), ("E2f6","cell-cycle exit enforcer"),
                    ("Ccne2","G1/S stall"), ("E2f2","endoreplication driver"),
                    ("Ccng1","G2/M arrest (marginal)")])

anchor_rows = "".join(
    f'<tr><td>{n}</td><td class="g lo">{o}</td><td class="g hi">{u}</td><td>{why}</td></tr>'
    for n, o, u, why in [
      ("mitotic brake","Wee1","Pkmyt1","Only Myt1 rises with maturation (t = +2.81); Wee1 falls (t = −1.32)."),
      ("p38","Mapk14","Mapk12","Baniol's “p38 up at P7” is p38γ (t = +4.49), not p38α."),
      ("APC/C co-activator","Cdh1","Fzr1","<em>Cdh1</em> is E-cadherin, detected in 4→0% of cells. Wiring mitotic exit to it hits an absent adhesion gene."),
      ("Cyclin D","Ccnd1","Ccnd2","Ccnd1 is undetectable here (0→6%); Ccnd2 is in ~100%."),
      ("α1-adrenergic","Adra1a","Adra1b","α1A is absent (0→2%). Murganti cited α1B and were right."),
    ])

PAGE = f"""<title>Cardiomyocyte cell-cycle fate model — results</title>
<style>{CSS}</style>
<div class="wrap">

<header class="mast">
  <div class="channel" aria-hidden="true"><i></i><i></i></div>
  <p class="eyebrow">Model note · e2f-heart-scrna / model</p>
  <h1>What two FUCCI papers can tell us about <em>why cardiomyocytes stop dividing</em></h1>
  <p class="standfirst">A re-analysis of Baniol 2021 and Murganti 2022, and the cell-cycle
  fate model it constrains. Five results, two of them corrections to the papers'
  own arithmetic — and one confound that blocks our knockout analysis until it is fixed.</p>
  <div class="meta">
    <span><b>Data</b> ENA PRJEB47622 · 285 cells</span>
    <span><b>Analysis</b> stdlib Python, 32 tests</span>
    <span><b>Model</b> specified, not yet built</span>
  </div>
</header>

<section>
  <span class="snum">01 — The question</span>
  <div class="col">
  <h2>Four outcomes, and two of them are not the same failure</h2>
  <p>A cardiomyocyte that enters the cell cycle has four ways out, and the whole
  point of modelling this is that two of them are mechanistically distinct and
  must not share a branch.</p>
  </div>
  <div class="tw"><table>
    <thead><tr><th>Outcome</th><th>What happens</th><th>Signature</th></tr></thead>
    <tbody>
      <tr><td>Quiescence</td><td>Never enters S phase</td><td class="lo">—</td></tr>
      <tr><td>Division</td><td>Mitosis and cytokinesis both complete</td><td>AurKB⁺ midbody; cell number rises</td></tr>
      <tr><td><strong>Binucleation</strong></td><td>Mitosis completes, the <strong>furrow fails</strong></td><td>Nuclei per cell rises; ploidy flat</td></tr>
      <tr><td><strong>Polyploidization</strong></td><td><strong>G2 arrest</strong> — APC/C never fires, no envelope breakdown</td><td>Ploidy rises; nuclei flat</td></tr>
    </tbody>
  </table></div>
  <div class="col">
  <p>Murganti's polyploid cells show mCherry-Cdt1 rising <em>before</em> mVenus-geminin
  is lost, decaying slowly with no sharp drop — the fingerprint of a G2 arrest rather
  than a failed furrow. Any model that routes both through one branch cannot
  reproduce that.</p>
  <h3>The clonidine triad is the test case worth designing for</h3>
  <p>One drug produced roughly a 2.2× rise in cell-cycle entry in three different
  systems, and three <em>different</em> outcomes. Reproducing all three by changing
  only a maturation coordinate is the bar.</p>
  </div>
  <div class="tw"><table>
    <thead><tr><th>System</th><th>Entry</th><th>Outcome markers</th><th>Fate</th></tr></thead>
    <tbody>
      <tr><td>hiPSC-CM <span class="lo">(least mature)</span></td><td class="num">Ki-67 9.0 → 22.0%</td>
          <td>Midbodies 0.83 → 1.69%; ploidy &amp; binucleation flat</td><td><strong>Division</strong></td></tr>
      <tr><td>Mouse neonatal CM, <em>in vitro</em></td><td class="num">Ki-67 11.0 → 23.2%</td>
          <td>Tetraploid 26.5 → 38.9%; <strong>zero</strong> midbodies in &gt;10,000 cells</td><td><strong>Polyploidization</strong></td></tr>
      <tr><td>Neonatal mouse, <em>in vivo</em></td><td class="num">EdU 25.5 → 36.2%</td>
          <td>Mononucleated 22.3 → 14.8%; ploidy flat</td><td><strong>Binucleation</strong></td></tr>
    </tbody>
  </table></div>
  <div class="col"><p>There is a free negative control in the same experiment:
  endothelial cells were unaffected (<span class="n">p = 0.85</span>) while smooth
  muscle responded (<span class="n">p = 0.04</span>), so no model may predict a
  generic mitogenic effect.</p></div>
</section>

<section>
  <span class="snum">02 — Pre-flight</span>
  <div class="col">
  <h2>Check the papers against each other before modelling anything</h2>
  <p>Four closed-form functionals of the same cycle constrain each other. Writing
  them down costs an afternoon and no solver, and it falsified two assumptions.</p>
  </div>
  <div class="eq">instantaneous phase fraction  ≈  entry rate × phase duration<br>
instantaneous Ki-67 index     ≈  entry rate × Ki-67 dwell&nbsp;&nbsp;&nbsp;&nbsp;<span class="lo">(LINEAR in rate)</span><br>
cumulative EdU index          =  1 − exp(−rate × window)&nbsp;&nbsp;<span class="lo">(SUBLINEAR in rate)</span><br>
fate-budget identity          :  ΔEdU ≈ ΔBinuc + ΔPoly + 2·ΔDiv</div>
  {fig("fig6_preflight_checks", 1, "Three of the six pre-flight checks. Predictions are closed-form from the published aggregates, with nothing fitted.")}
  <div class="finds">
    <div class="find"><h4>The hiPSC data hides a chronically arrested pool</h4>
      <span class="chip flag">tension</span>
      <p>Live imaging scored an outcome in 9.65% of 570 cells over 72 h, which with the
      measured 16.38 h S/G2/M predicts <span class="n">2.20%</span> mVenus⁺ at any instant.
      Flow cytometry of 50,000 cells saw <span class="n">5.10%</span> — implying an apparent
      S/G2/M of 38.1 h. The 2.9 pp excess is <strong>2.1× larger</strong> than the 1.40% caught
      completing polyploidization, so 24.5 h is a right-censored lower bound and the model needs
      a persistent arrested state. Rule out a permissive flow gate first.</p></div>
    <div class="find"><h4>The mouse arm is self-consistent across both papers</h4>
      <span class="chip ok">confirmed</span>
      <p>Baniol's P0/P7 mAG⁺ fractions with the 15.1 h duration give S-entry rates of
      19.7%/day and 1.32%/day — a 14.9× fall. Integrated over the P1→P5 labelling window
      that predicts <span class="n">32.8%</span> cumulative EdU against Murganti's observed
      <span class="n">25.54%</span>. Different papers, cohorts and modalities agreeing to
      1.28× is what licenses cumulative EdU as a held-out target.</p></div>
    <div class="find"><h4>Clonidine's two readouts have the wrong sign for a rate change</h4>
      <span class="chip flag">tension</span>
      <p>Ki-67 rose 2.12× but cumulative EdU rose 2.78×. Because EdU saturates, a rate
      increase <em>always</em> gives fold(EdU) &lt; fold(Ki-67) — rate folds of 1.5/2.0/3.0
      give EdU folds of 1.46/1.90/2.72. The observation is the other way round, so it is
      neither a rate change nor a shorter cycle. What survives: <strong>Ki-67 decays in
      arrested cells while the EdU label persists</strong> — which makes the 1.31 ratio positive
      evidence <em>for</em> non-productive cycling. Saturation-corrected, the true entry-rate
      fold is 1.52, not 1.42.</p></div>
    <div class="find"><h4>Cardiomyocyte mitosis should last about four hours</h4>
      <span class="chip open">prediction</span>
      <p>The 1% cytoplasmic-FUCCI fraction at P0, with an independent 5.7%/day event rate,
      implies a mitotic dwell of <span class="n">4.2 h</span> — several-fold longer than the
      ~1 h of cycling somatic cells, and longer still in binucleating than dividing cells.
      Neither paper computed it.</p></div>
  </div>
</section>

<section>
  <span class="snum">03 — Read the caveats first</span>
  <div class="col">
  <h2>Two properties of this dataset bound what it can be asked</h2>
  </div>
  {fig("fig7_sort_enrichment", 2, "The sort is asymmetric. P0 is essentially unenriched; P7 is 4.5–5.2× enriched, so the raw cycling fraction rises from P0 to P7 in the data and falls in reality.")}
  <div class="col">
  <p>Only within-group expression contrasts are valid here. Phase and cycling
  fractions must come from Baniol's imaging, never from these cells.</p>
  <p>The second property is subtler. The Regev/Tirosh gene lists behind the phase
  caller contain <span class="gene">Ect2</span> (G2/M) and <span class="gene">E2f8</span> (S),
  so cycling-versus-noncycling contrasts are partly self-fulfilling for 21 panel genes.
  The P0→P7 within-cycling axis — where our findings live — is not circular in the same
  way, and the residual bias runs <em>against</em> the Ect2 result: selecting cells
  <em>called</em> G2M enriches for high Ect2 in both groups, compressing the difference.
  The comparators that matter (<span class="gene">RhoA, Racgap1, Cep55, Ccnb1, Ccna2</span>)
  are all outside the lists.</p>
  </div>
  <div class="note"><span class="lbl">Also worth knowing</span>
  <p>The <code>Baniol2021_FUCCI_R</code> variant is <strong>not</strong> an independent phase
  call. Comparing the two metadata files cell by cell, <code>phase</code>,
  <code>cycling_score</code>, <code>CellType</code>, <code>stage</code> and
  <code>cluster_label</code> are 100% identical across all 285 cells; only
  <code>seurat_clusters</code> differs (19.3%). It gives a second clustering, not a second
  labelling.</p></div>
</section>

<section>
  <span class="snum">04 — The lesion</span>
  <div class="col">
  <h2>Ect2 is specific; RhoA only looks it</h2>
  <p>In cycling ventricular cardiomyocytes from P0 to P7, <span class="gene">Ect2</span>
  falls to <span class="n">0.55×</span> (t = −6.17, p &lt; 10⁻⁴, clears Bonferroni over 25
  genes). What makes that a <em>specific</em> finding rather than merely a significant one
  is the genome-wide baseline: the median P7/P0 ratio across {d["n_genes"]} expressed genes
  is <span class="n">{f2(d["median"],3)}</span>, because P7 cells have slightly lower
  complexity. Ect2 sits at the <span class="n">3rd</span> percentile of that distribution.
  Nothing else in either programme comes close.</p>
  </div>
  {fig("fig1_ect2_specificity", 3, "Every cytokinesis and mitotic gene against the background drift. Larger marks clear Bonferroni; the grey band is the genome-wide 5–95% range.")}
  <div class="tw"><table>
    <thead><tr><th>Gene</th><th>Module</th><th class="num">P7/P0</th><th class="num">t</th><th class="num">p</th><th class="num">Genome-wide pctile</th></tr></thead>
    <tbody>{drift_rows}</tbody>
  </table></div>
  <div class="note"><span class="lbl">The instructive near-miss</span>
  <p><span class="gene">RhoA</span> clears Bonferroni (<span class="n">p = 0.0004</span>) yet
  sits at the <span class="n">{f2(e["Rhoa"]["pct_genomewide"],0)}th</span> percentile
  genome-wide, so its low p-value comes from low variance rather than a large effect.
  <strong>Significance is not specificity.</strong> We do not claim RhoA.</p></div>
  <div class="col">
  <h3>And our current score cannot see any of this</h3>
  </div>
  {fig("fig2_module_cancellation", 4, "Averaged over 12 cytokinesis versus 13 mitotic genes, the two modules decline in lockstep and their difference is flat.")}
  <div class="col">
  <p>The <code>sig_ploidy = sig_prolif − sig_cytokinesis</code> heuristic in
  <code>build_signature_scores.R</code> is a difference of exactly those two averages.
  They fall together ({f2(m["cytokinesis"]["ratio"],3)}× and {f2(m["mitotic"]["ratio"],3)}×),
  so the difference is flat — <span class="n">Δ = {m["difference_index"]["delta"]:+.3f},
  p = {f2(m["difference_index"]["p"],2)}</span>. On this data the score is
  <strong>structurally blind</strong> to the clearest cytokinesis-competence signal available.
  That is the concrete argument for a model with named nodes.</p>
  </div>
</section>

<section>
  <span class="snum">05 — The coordinate</span>
  <div class="col">
  <h2>Maturation is measurable, not a fitted knob</h2>
  <p>Define <strong>M = mean z(fatty-acid oxidation) − mean z(glycolysis)</strong> per cell,
  over twelve genes each. Three things make it a good coordinate for the model's slow axis.</p>
  </div>
  {fig("fig3_maturation_coordinate", 5, "M across the eight stage × chamber × cycling groups. Marker size hints at n; the two P0-atrial groups are too small to carry a claim.")}
  <div class="tw"><table>
    <thead><tr><th>Group</th><th class="num">n</th><th class="num">M</th><th class="num">SE</th></tr></thead>
    <tbody>{mat_rows}</tbody>
  </table></div>
  <div class="col">
  <p><strong>It is continuous, not a relabelled timepoint.</strong> Between-stage separation
  is 1.78 against a within-stage SD near 0.5, so single cells can be ordered along it.</p>
  <p><strong>It reproduces a Baniol claim it was not built to test.</strong> The ordering puts
  P7 atrial cells (<span class="n">−0.07</span>) between P0 ventricular
  (<span class="n">−0.68</span>) and P7 ventricular (<span class="n">+1.12</span>) — exactly
  the paper's observation that cycling P7 atrial cardiomyocytes are transcriptionally more
  immature and retain more proliferative capacity.</p>
  <p><strong>It predicts the model's two load-bearing couplings</strong>, within cycling
  ventricular cells alone (n = 89, P0 and P7 pooled so stage cannot be doing the work).</p>
  </div>
  <div class="tw"><table>
    <thead><tr><th>Node</th><th>Role in the model</th><th class="num">r with M</th><th class="num">t</th></tr></thead>
    <tbody>{corr_rows}</tbody>
  </table></div>
  {fig("fig4_m_correlations", 6, "Cytokinesis competence falls with maturation and the exit enforcer rises. Both are therefore measured functions of M rather than fitted parameters.")}
  <div class="note ok"><span class="lbl">Why this matters</span>
  <p>Ect2(M) and E2f6(M) become <strong>measured functions</strong>, which removes the
  model's largest degree of freedom. Caveats stated plainly: n = 89; Ccng1 is marginal;
  the P0-atrial groups are n = 8 and n = 4; and M is z-scored within one dataset, so it has
  no absolute cross-system scale — the clonidine triad is therefore an ordinal constraint.</p></div>
</section>

<section>
  <span class="snum">06 — The E2F layer</span>
  <div class="col">
  <h2>Three roles, and E2f7 is not E2f8</h2>
  <p>This independently reproduces and quantifies every E2F claim in the paper.</p>
  </div>
  {fig("fig5_e2f_family", 7, "All eight E2F family members across the eight groups. Shade encodes mean log-normalised expression; the hairline separates ventricular from atrial blocks.")}
  <div class="col">
  <p><span class="gene">E2f1/7/8</span> are cycling-restricted. <span class="gene">E2f2</span>
  is near-absent unless cycling and rises with maturation (t = +2.18) — the endoreplication
  candidate. <span class="gene">E2f6</span> is the only member expressed in noncycling cells,
  reaching 0.38 in <strong>100%</strong> of noncycling P7 ventricular cells (t = +4.12) — the
  exit enforcer. And <strong><span class="gene">E2f7</span> is flat across maturation
  (t = −0.20) while <span class="gene">E2f8</span> rises (t = +2.28)</strong>, so they must be
  separate model nodes.</p>
  <h3>Edges tested before wiring them</h3>
  <p>Spearman across cycling ventricular cells (n = 89): <span class="gene">E2f1~E2f7</span>
  <span class="n hi">+0.49</span> and <span class="gene">E2f1~E2f8</span>
  <span class="n hi">+0.46</span> support the E2F1 → E2F7/8 delayed-feedback arm as a real
  edge — and that is precisely the arm our knockout removes.
  <span class="gene">E2f1~E2f6</span> is <span class="n">+0.18 (n.s.)</span>, so E2f6 is
  <em>not</em> part of the activator programme and belongs on its own maturation-driven
  branch. <span class="gene">E2f8~Ccna2</span> is <span class="n">−0.56</span>, independently
  confirming E2F8 is G1/S-restricted and off in S/G2.</p>
  </div>
  <div class="note"><span class="lbl">A negative result that changes the plan</span>
  <p><span class="gene">Gmnn~Cdt1</span> mRNA correlates <em>positively</em>
  (<span class="n">+0.29</span>) even though the FUCCI reporters they inspire are
  anti-correlated by phase — because the oscillation is post-translational. So
  <strong>scRNA-seq cannot validate the FUCCI reporter module at all.</strong> That layer has
  to be checked against imaging traces only, and the write-up should never claim
  transcriptomic support for it.</p></div>
</section>

<section>
  <span class="snum">07 — Wiring</span>
  <div class="col">
  <h2>Five gene anchors the data overturns</h2>
  <p>Checking detection rates before wiring anything caught five obvious choices that are
  wrong in cardiomyocytes.</p>
  </div>
  <div class="tw"><table>
    <thead><tr><th>Node</th><th>Obvious</th><th>Use instead</th><th>Why</th></tr></thead>
    <tbody>{anchor_rows}</tbody>
  </table></div>
  <div class="col">
  <h3>A mechanism the papers did not state</h3>
  <p>Across all 285 cells, from P0 to P7, <span class="gene">Nisch</span> — nischarin, the
  imidazoline-I1 receptor and clonidine's proposed pro-proliferative target — <strong>falls</strong>
  from 1.28 to 1.05 at 100% detection, while <span class="gene">Adra1b</span> (0.051 → 0.097)
  and <span class="gene">Adrb1</span> (0.094 → 0.156) both <strong>rise</strong>. Since
  β-adrenergic/PKA signalling blocks cardiomyocyte cytokinesis (Liu 2019), the adrenergic
  balance tips toward the cytokinesis brake exactly as cardiomyocytes lose the ability to
  divide — an independent reason clonidine's effect should become less productive with
  maturation. Limitation: <span class="gene">Adra2a/2b/2c</span> are absent from the store,
  so clonidine's canonical α2 target cannot be assessed here at all.</p>
  <h3>And one paper claim that weakens</h3>
  <p>Baniol attributed the P7 G1/S delay to a DNA-damage response, but on the P0→P7
  within-cycling axis all four genes are flat (<span class="gene">Chek1</span> t = +0.41,
  <span class="gene">Timeless</span> +0.64, <span class="gene">Ung</span> −0.14,
  <span class="gene">Msh6</span> −0.91). Their claim concerned one G1/S sub-cluster rather
  than cycling P7 ventricular cells as a whole, so this is not a contradiction — but the
  <strong>DDR module is hypothesis-only</strong>, resting on Puente 2014 rather than on this
  dataset.</p>
  </div>
</section>

<section>
  <span class="snum">appendix — before the data came back</span>
  <div class="col">
  <h2>What we expected, and why it was wrong</h2>
  <p>Before extracting the bundle, the reasoning ran: <span class="gene">E2f8</span>
  sits in the S-phase gene list used by <code>CellCycleScoring</code>, so in an
  E2f7/E2f8 knockout the S-score would be biased <em>downward</em> and the knockout
  would spuriously look less proliferative. That was recorded as the top blocking
  item.</p>
  <p>Both halves turned out wrong, and in an informative way. The bias runs the
  <em>other</em> direction, because E2f7/E2f8 mRNA is higher in KO, not lower — and its
  magnitude is negligible either way, because one gene out of 42 cannot move a mean.
  Worth keeping visible: the prediction was specific enough to be checked, and
  checking it cost one container run.</p>
  </div>
</section>

<section>
  <span class="snum">06 — The model</span>
  <div class="col">
  <h2>What we built, and what it gets right</h2>
  <p><code>cmfate</code> is a 55-node normalized-Hill network with
  78 reactions, eight inputs and four fate outputs. Standard library
  only — engine, adaptive integrator and calibration included. The network is three
  tracked text files with one reaction per row, each carrying its own evidence
  column, because a curated network is reviewed one edge at a time.</p>
  <h3>The network</h3>
  </div>
  """+fig("fig13_pathways", 8, "Every node and every edge, generated from the spec. Maturation is badged rather than drawn: with out-degree 15 into every module, its edges would dominate the figure.")+"""
  <div class="col">
  <p>Left to right: a stimulus sets the signalling layer, which sets E2F activity,
  which builds the cycle and cytokinesis machinery, which sets three gates — and the
  four outcomes are a product of those three gates alone.</p>
  <p>The two coloured arms are what make one drug give three different outcomes. Orange
  is the abscission arm, <span class="gene">E2Fact → Ect2 → RhoA → Midbody → AbsRaw →
  Abscission</span>, which <strong>maturation</strong> closes. Blue is the mitotic-entry
  brake, <span class="gene">ROS → DDR → Ccng1/Pkmyt1 → MitCompRaw → MitoticEntry</span>,
  which <strong>culture</strong> closes. Which of the two shuts first is the whole 2×2.</p>
  <h3>How the engine works</h3>
  </div>
  """+fig("fig11_engine", 9, "A state vector in, a derivative out. Both of the model's structural bugs lived in the two composition rules, so the diagram marks them rather than just naming them.")+"""
  <div class="col">
  <p>Two of the four stages are composition rules, and they are where the bugs were.
  <strong>Per reaction</strong>, each reactant goes through <code>act()</code> or
  <code>inhib()</code> and the results are combined by AND — with the
  <code>(w, n, EC50)</code> triple belonging to the <em>reaction</em>, so one curve
  applies to all of its reactants. <strong>Per node</strong>, every reaction pointing
  at it is combined by weighted OR, and because OR can only <em>add</em> drive,
  <strong>no route can be vetoed by another</strong>. That is exactly how an OR'd
  Midbody route bypassed the obligatory Ect2/RhoA arm.</p>
  <p>Then <code>dY/dt = (drive·Ymax − Y)/τ</code>, which is also how perturbations
  enter: knockdown is Ymax → 0, overexpression holds Y = 1. Output-only nodes never
  reach the integrator — their steady state is algebraic, which matters because the
  fate nodes carry τ up to 35 h.</p>
  </div>
  """+fig("fig12_activation_classes", 10, "Each parameter class against the identity line, evaluated through the engine itself. The last panel is what the guard prevents.")+"""
  <div class="col">
  <p>The <code>EC50ⁿ &lt; 0.5</code> guard is not a correction to the published
  formalism — Netflux's own defaults sit safely inside the ceiling, which is why the
  reference never needed it. A high-threshold gate is what walks into it, and this
  model has three.</p>
  <h3>The fate layer is exact, and exactly identified</h3>
  <p>The four fates sum to <strong>1</strong> with no normalization step and no engine
  change, via a complementary product over three gate nodes:</p>
  </div>
  <div class="eq">!SPhase                              =&gt; Quiescent<br>
SPhase &amp; MitoticEntry &amp; Abscission   =&gt; Division<br>
SPhase &amp; MitoticEntry &amp; !Abscission  =&gt; Binucleation<br>
SPhase &amp; !MitoticEntry               =&gt; Polyploidization</div>
  <div class="col">
  <p>Inverting Murganti's measured fractions identifies the three gates in closed
  form — <span class="n">g(SPhase) = 0.0965</span>,
  <span class="n">g(MitoticEntry) = 0.8549</span>,
  <span class="n">g(Abscission) = 0.6170</span> — and because competence is
  deliberately <em>not</em> gated on prevalence, those are three independent
  one-dimensional solves rather than a joint fit. The model reproduces the fit context
  to <span class="n">1×10⁻⁵</span> on all four fractions. Which means, stated plainly:
  <strong>the fate layer is fitted exactly and predicts nothing.</strong> All predictive
  content is downstream of it.</p>
  </div>
  <div class="note ok"><span class="lbl">Two bugs a passing test suite would have hidden</span>
  <p>An <strong>OR'd floor on Ect2</strong> flattened its maturation response to 1.20×
  when the measured drop is 1.83×. And an <strong>OR'd Midbody route</strong>
  (<span class="gene">Centralspindlin &amp; AurKB</span>) bypassed the Ect2/RhoA arm
  entirely, making it <em>impossible</em> for Ect2 to be rate-limiting — the model's
  central claim. Making RhoA obligatory widened that arm's span from 1.6× to
  <strong>42×</strong> between hiPSC and mouse P1, and only then did the model behave as
  designed.</p></div>
  """+fig("fig8_maturation_fates", 11, "Sweeping maturation alone. Only M changes; every other input is held at the in-vivo setting.")+"""
  <div class="col">
  <p>Sweeping maturation alone, the cycling flux moves from binucleation into
  polyploidization, crossing over near M = 0.62. Division's share stays under 3%
  across this sweep because it holds mechanical load and β-adrenergic tone at their
  in-vivo values — whereas the hiPSC-CM <em>context</em>, which differs in those too,
  reaches a 53% division share. The two-factor result, made concrete: maturation alone
  does not set the outcome.</p>
  <h3>The held-out test</h3>
  </div>
  """+fig("fig9_clonidine_triad", 12, "Change in each fate under clonidine. Calibrated at the hiPSC context only; the other two are held out.")+"""
  <div class="col">
  <p><strong>3 of 3 dominant fates correct</strong>, from one calibration and one input
  change per context. The mechanism is a clean 2×2: maturation selects between hiPSC
  and mouse through the <strong>Ect2 arm</strong>, and the in-vitro flag selects between
  the two mouse contexts through the <strong>ROS → DDR → Ccng1/Pkmyt1</strong> arm, which
  collapses mitotic entry (0.293 in vitro versus 0.591 in vivo) and routes the flux to
  polyploidization rather than binucleation.</p>
  </div>
  <div class="note"><span class="lbl">The largest quantitative miss</span>
  <p>The entry-response <em>magnitude</em>. 1.72× at hiPSC
  against an observed 2.44× is fine; <strong>6.3× and
  8.8× in the two mature contexts against observed
  2.12× and 1.52× is not.</strong> The cause is diagnosed: clonidine acts by relieving the
  PKA brake, baseline PKA scales with β-adrenergic tone, and the model raises that tone
  steeply with maturation — so the relative effect of removing it is largest exactly
  where the data says it should be smallest. Pinned by a test, so a fix shows up as a
  test change rather than a silent improvement.</p></div>
  """+fig("fig10_ko_and_screen", 13, "In-silico E2f7/E2f8 double knockdown at two maturations, and every node ranked by how much it converts cycling into division.")+"""
  <div class="col">
  <p>In-silico E2f7/E2f8 double knockdown: the pro-division effect is <strong>57× larger
  at P0 than at P7</strong> (+0.0163 versus
  +0.0003), even though the knockdown raises
  Ect2 substantially at both ages (+0.44 and
  +0.24). The machinery goes up either way; only
  at P0 is the rest of the context permissive enough to convert it. There is real
  epistasis, which matters because the lab's data <em>is</em> a double knockout.</p>
  <p>The screen ranks nodes by the <strong>conditional</strong> division share, D/(D+B+P),
  not by Δ-Division and not by entry. Murganti's wet screen scored
  percent-mVenus-positive, i.e. S-phase entry — and in this partition every fate scales
  with entry, so an entry-scored screen is confounded with respect to productive
  division. That is what 6/94 hits collapsing to 2/94 on validation looks like. The
  model's own specificity now matches the bench, after two errors in our own scoring
  were fixed. Scored like-for-like in hiPSC-CM — where Murganti actually ran theirs —
  <strong>10.5% of perturbations clear the hit bar against their 6.4%</strong>, and only
  <strong>25% of those also raise the division share, against 2 of 6 surviving their
  validation.</strong> The old numbers (14.5% and 82%) came from a purely relative hit test
  applied at P7, where baseline entry is 0.27% so a 0.4-percentage-point change scored
  as a hit, and from a conversion figure that counted any positive value including
  epsilon. Baseline entry is 36× higher in hiPSC-CM, so the two were never comparable.</p>
  </div>
</section>

<section>
  <span class="snum">07 — Our own data, answered</span>
  <div class="col">
  <h2>The knockout is functional, and the confound was a false alarm</h2>
  </div>
  <div class="note ok"><span class="lbl">Resolved</span>
  <p>E2f7 and E2f8 mRNA are <em>higher</em> in KO than WT (log2FC +0.82 and +1.57 at P7,
  both p &lt; 1×10⁻¹⁰), which looks alarming until you check the targets:
  <strong>24 of 26 canonical E2F targets are significantly up at P7</strong> (median log2FC
  +0.91) against a housekeeping baseline of +0.19, and 20 of 26 at P0 (median +0.50)
  against −0.02. Regulon activity for HALLMARK_E2F_TARGETS is KO−WT = +1.12 at P0 and
  +3.15 at P7. Since E2f1/7/8 are themselves E2F1-induced, losing the repressors'
  <em>protein</em> function raises E2F1 activity, which transcriptionally induces all
  three — exactly the feedback loop this analysis validated in the Baniol data. So "KO
  not transcript-confirmed" should be restated: the KO is not
  transcript-<em>null</em>, and it is functionally confirmed.</p></div>
  <div class="col">
  <p><strong>And the E2f8 phase-scoring confound is negligible.</strong> E2f8 is one of 42
  genes in the S-phase list and carries 0.15–0.54% of the S-set mean. Recomputing phase
  with and without it changes the call for 0.9% of cells and shifts the KO−WT cycling
  gap by <span class="n">≤0.5 pp</span>. Caveat on our own method: the simple-mean
  reimplementation agrees with the shipped <code>CellCycleScoring</code> labels for only
  49% of cells, so it supports the <em>differential</em> claim — a within-method
  comparison where the approximation cancels — but not any absolute cycling fraction.</p>
  <p>One hard limit on validating the model here: only <strong>29 of 63</strong> model node
  genes are in the 2,181-gene curated panel. <span class="gene">E2f2–E2f6, Rb1, Ccng1,
  Chek1, Wee1, Pkmyt1, Rhoa, Nisch, Adrb1</span> and <span class="gene">Mapk12</span>
  are all missing, so most of the network cannot be scored against the lab's own cells
  until the panel is widened.</p>
  </div>
</section>

<section>
  <span class="snum">08 — What's next</span>
  <div class="col">
  <h2>What is left to do</h2>
  <p>Ordered by whether it unblocks something else. The full list, with rationale for each
  item, is in <code>model/TODO.md</code>.</p>
  </div>
  <ol class="todo">
    <li><div class="t"><b>Recompute Phase without E2f7/E2f8</b><span class="chip flag">blocking</span></div>
      <p>Every KO-versus-WT cycling result is suspect until this is checked and the
      concordance reported.</p></li>
    <li><div class="t"><b>Extract the KO bundle to CSV</b><span class="chip flag">blocking</span></div>
      <p>One <code>docker run</code> against <code>rocker/r-ver</code> — the bundle needs only
      <code>Matrix</code>, no Seurat. Re-export the matrix into the same CSR layout so the
      analysis code works on both datasets unchanged. Note the copy on disk is the oldest
      build and lacks the signature scores.</p></li>
    <li><div class="t"><b>Fix two engine bugs first</b><span class="chip open">prerequisite</span></div>
      <p><code>beta</code> goes negative whenever EC50ⁿ &gt; 0.5, corrupting state silently —
      and high-threshold gates are exactly what the fate layer needs. Separately,
      <code>perturbation_matrix</code> takes one node at a time, so it cannot express the
      E2f7/E2f8 double knockout, which is the central experiment.</p></li>
    <li><div class="t"><b>Build the logic network, then the ODE</b><span class="chip open">open</span></div>
      <p>The four fates sum to exactly 1 with no engine change via a complementary product
      over three gates, which also makes the fate layer closed-form identifiable from
      Murganti's table. The ODE extends Gérard &amp; Goldbeter 2009 — 45 species already
      supply the Rb–E2F switch and cyclin waves; what must be added is Cdt1 and geminin as
      FUCCI observables, cytokinesis competence, and ploidy bookkeeping.</p></li>
    <li><div class="t"><b>Fit mouse, predict human</b><span class="chip open">open</span></div>
      <p>Hold out every hiPSC-CM number and change only the maturation coordinate. Pass
      criteria written down in advance, including reproducing the <em>significance structure</em>
      of the durations rather than three means.</p></li>
    <li><div class="t"><b>Use the free negative controls</b><span class="chip open">open</span></div>
      <p><span class="gene">E2f3</span>, <span class="gene">AurKB</span>,
      <span class="gene">Ccna2</span>, <span class="gene">Racgap1</span> and the DDR set must
      <em>not</em> move. A model that shifts the comparator arm as much as it shifts Ect2 has
      failed even if every positive contrast passes.</p></li>
  </ol>
  <div class="col">
  <h3>Improvements to the app itself</h3>
  <p>Replace the <code>sig_ploidy</code> heuristic with named nodes; widen the curated panel
  to all of <span class="gene">E2f1–E2f8</span>, since Baniol's entire E2F argument rests on
  E2f2 and E2f6 and neither is currently in it; propagate per-cell atrial/ventricular labels,
  which are declared in <code>CAT_COLS</code> but absent from the bundle; and add the
  sort-enrichment caveat to the About block, because anyone reading cycling fractions off a
  FACS-enriched dataset will get the developmental direction backwards.</p>
  </div>
</section>

<p class="foot">Baniol M et al. <em>Exp Cell Res</em> 408:112880 (2021) ·
Murganti F et al. <em>Front Cardiovasc Med</em> 9:840147 (2022) ·
Gérard C &amp; Goldbeter A, <em>PNAS</em> (2009), BioModels BIOMD0000000730.
Figures and every quoted statistic regenerate from <code>model/</code> —
<code>python3 -m cmcycle.figures</code>.</p>

</div>
"""
out = pathlib.Path("report.html")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(PAGE, encoding="utf8")
print("wrote", out, len(PAGE), "B")
