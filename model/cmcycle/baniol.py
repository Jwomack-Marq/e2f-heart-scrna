"""Re-analysis of the Baniol et al. 2021 FUCCI scRNA-seq (ENA PRJEB47622).

Every number quoted in ``model/RESULTS.md`` is produced here, so the write-up is
reproducible rather than transcribed. Stdlib only: the expression store is a raw
gene-major CSR dump, read with ``array`` + ``seek``.

Reference: Baniol M, Murganti F, Smialowska A, Panula J, Lazar E, Brockman V,
Giatrellis S, Derks W, Bergmann O. Identification and characterization of distinct
cell cycle stages in cardiomyocytes using the FUCCI transgenic system.
Exp Cell Res 408:112880 (2021). doi:10.1016/j.yexcr.2021.112880

The dataset is 285 FACS-sorted cardiomyocytes from P0 and P7 mouse hearts,
Smart-seq2, with per-cell atrial/ventricular and cycling/noncycling labels.

IMPORTANT -- two properties of this dataset constrain what may be asked of it:

1. The sort is *asymmetrically* enriched. P0 is essentially unenriched but P7 is
   4.5-5.2x enriched for cycling cells, so the raw cycling fraction RISES from P0
   to P7 in the data and FALLS in reality. Cycling fractions here are never a
   valid rate; only within-group expression contrasts are.
2. The ``phase`` call is not independent of several genes of interest. The
   Regev/Tirosh lists used by ``scanpy.tl.score_genes_cell_cycle`` contain Ect2
   (G2/M) and E2f8 (S), among others -- see :data:`PHASE_S` / :data:`PHASE_G2M`
   and :func:`circularity_audit`.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import os
from array import array
from dataclasses import dataclass, field

# --------------------------------------------------------------------------- #
# Where the data lives. Not in this repo -- it is a 40 GB sibling store.
# --------------------------------------------------------------------------- #
DEFAULT_ROOT = "/home/justin/Projects/lab-server/apps/cardiac-rnaseq-explorer"
DATASET = "Baniol2021_FUCCI"


def _root() -> str:
    return os.environ.get("CARDIAC_RNASEQ_ROOT", DEFAULT_ROOT)


# --------------------------------------------------------------------------- #
# The phase-calling gene lists, verbatim from the explorer's preprocess/ingest.py.
# Vendored here because the circularity audit must not silently drift if that
# pipeline is edited -- a copy that disagrees is itself the finding.
# --------------------------------------------------------------------------- #
PHASE_S = set("""MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1
UHRF1 HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2
RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1
E2F8""".split())

PHASE_G2M = set("""HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B
MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E
TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5
CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3
GAS2L3 CBX5 CENPA""".split())


# --------------------------------------------------------------------------- #
# Gene panels
# --------------------------------------------------------------------------- #
CYTOKINESIS = ["Ect2", "Anln", "Racgap1", "Kif23", "Cep55", "Aurkb", "Cit",
               "Kif20b", "Prc1", "Cdca8", "Incenp", "Rhoa"]
MITOTIC = ["Ccnb1", "Cdk1", "Cdc20", "Plk1", "Bub1", "Aurka", "Cenpa", "Cenpe",
           "Cenpf", "Ube2c", "Top2a", "Tpx2", "Nusap1"]
E2F_FAMILY = [f"E2f{i}" for i in range(1, 9)]

#: Fatty-acid oxidation, high in mature cardiomyocytes.
FAO = ["Hmgcs2", "Fabp3", "Pdk4", "Ucp2", "Cpt1a", "Acadm", "Acadvl", "Hadha",
       "Hadhb", "Cd36", "Ech1", "Decr1"]
#: Glycolysis, high in immature cardiomyocytes.
GLYCOLYSIS = ["Hk1", "Hk2", "Pfkl", "Pfkm", "Pkm", "Ldha", "Slc2a1", "Aldoa",
              "Gapdh", "Eno1", "Pgk1", "Tpi1"]

#: Anchors the data overturns: (node, obvious choice, what to use instead).
ANCHOR_CHECKS = [
    ("mitotic brake", "Wee1", "Pkmyt1"),
    ("p38", "Mapk14", "Mapk12"),
    ("APC/C co-activator", "Cdh1", "Fzr1"),
    ("Cyclin D", "Ccnd1", "Ccnd2"),
    ("alpha1-adrenergic", "Adra1a", "Adra1b"),
]

ADRENERGIC = ["Adra1a", "Adra1b", "Adra1d", "Adra2a", "Adra2b", "Adra2c",
              "Adrb1", "Adrb2", "Adrb3", "Nisch"]

DDR = ["Chek1", "Timeless", "Ung", "Msh6"]


# --------------------------------------------------------------------------- #
# The store
# --------------------------------------------------------------------------- #
@dataclass
class Store:
    """Gene-major CSR expression store plus per-cell metadata.

    Layout (all little-endian, no headers): ``indptr.bin`` int32 of length
    n_genes+1, ``indices.bin`` int32 cell indices, ``data.bin`` float32
    log-normalised values. One gene's vector across all cells is a single
    seek+read, so nothing scales with the full matrix.
    """
    root: str = field(default_factory=_root)
    dataset: str = DATASET

    def __post_init__(self):
        e = os.path.join(self.root, "expr", self.dataset)
        d = os.path.join(self.root, "data", self.dataset)
        self._edir = e
        with open(os.path.join(e, "meta.json")) as fh:
            self.info = json.load(fh)
        self.genes = [l.strip() for l in open(os.path.join(e, "genes.txt"))]
        self.cells = [l.strip() for l in open(os.path.join(e, "cells.txt"))]
        self._gidx = {g: i for i, g in enumerate(self.genes)}
        self._cpos = {c: i for i, c in enumerate(self.cells)}
        self._indptr = array("i")
        self._indptr.frombytes(open(os.path.join(e, "indptr.bin"), "rb").read())
        with gzip.open(os.path.join(d, "meta.csv.gz"), "rt") as fh:
            self.meta = {r["ID"]: r for r in csv.DictReader(fh)}
        self._fi = open(os.path.join(e, "indices.bin"), "rb")
        self._fd = open(os.path.join(e, "data.bin"), "rb")
        self._check()

    def _check(self) -> None:
        """Five integrity assertions. These files are raw dumps with no header,
        so a truncated or byte-swapped copy would otherwise read as plausible."""
        n_g, n_c, nnz = self.info["n_genes"], self.info["n_cells"], self.info["nnz"]
        assert len(self.genes) == n_g, (len(self.genes), n_g)
        assert len(self.cells) == n_c, (len(self.cells), n_c)
        assert len(self._indptr) == n_g + 1
        assert self._indptr[0] == 0 and self._indptr[-1] == nnz
        for f in ("indices.bin", "data.bin"):
            assert os.path.getsize(os.path.join(self._edir, f)) == 4 * nnz, f
        import sys
        assert sys.byteorder == "little", "store is little-endian"

    def has(self, gene: str) -> bool:
        return gene in self._gidx

    def gene(self, gene: str) -> list[float]:
        """Dense vector for one gene across all cells (285 floats -- cheap)."""
        i = self._gidx[gene]
        s, e = self._indptr[i], self._indptr[i + 1]
        idx, dat = array("i"), array("f")
        self._fi.seek(s * 4); idx.frombytes(self._fi.read((e - s) * 4))
        self._fd.seek(s * 4); dat.frombytes(self._fd.read((e - s) * 4))
        v = [0.0] * len(self.cells)
        for j, val in zip(idx, dat):
            v[j] = val
        return v

    def group(self, **crit) -> list[int]:
        """Cell indices matching all metadata criteria, e.g. stage='P0'."""
        return [self._cpos[c] for c in self.cells
                if all(self.meta[c][k] == v for k, v in crit.items())]

    def groups8(self) -> dict[str, list[int]]:
        """The eight stage x chamber x cycling groups, in maturation order."""
        out = {}
        for st in ("P0", "P7"):
            for ct in ("aCM", "vCM"):
                for cy in ("cycling", "noncycling"):
                    out[f"{st}-{ct}-{cy}"] = self.group(
                        stage=st, CellType=ct, cycling_score=cy)
        return out


# --------------------------------------------------------------------------- #
# Statistics -- stdlib, so the module has no install step
# --------------------------------------------------------------------------- #
def mean_sd_se(v):
    n = len(v)
    if n == 0:
        return 0.0, 0.0, 0.0
    m = sum(v) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1)) if n > 1 else 0.0
    return m, sd, sd / math.sqrt(n)


def _betai(a, b, x):
    """Regularised incomplete beta, for t-distribution tail probabilities."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0

    def cf(a, b, x):
        qab, qap, qam = a + b, a + 1, a - 1
        c, d = 1.0, 1 - qab * x / qap
        d = 1 / d if abs(d) > 1e-30 else 1e30
        h = d
        for m in range(1, 300):
            m2 = 2 * m
            aa = m * (b - m) * x / ((qam + m2) * (a + m2))
            d = 1 + aa * d; d = 1 / d if abs(d) > 1e-30 else 1e30
            c = 1 + aa / c; c = c if abs(c) > 1e-30 else 1e30
            h *= d * c
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
            d = 1 + aa * d; d = 1 / d if abs(d) > 1e-30 else 1e30
            c = 1 + aa / c; c = c if abs(c) > 1e-30 else 1e30
            de = d * c; h *= de
            if abs(de - 1) < 3e-12:
                break
        return h

    lb = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
          + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return math.exp(lb) * cf(a, b, x) / a
    return 1 - math.exp(lb) * cf(b, a, 1 - x) / b


def welch(a, b):
    """Welch's unequal-variance t test. Returns (mean_a, mean_b, t, p)."""
    n1, (m1, s1, _) = len(a), mean_sd_se(a)
    n2, (m2, s2, _) = len(b), mean_sd_se(b)
    se = math.sqrt(s1 * s1 / n1 + s2 * s2 / n2)
    if se == 0:
        return m1, m2, 0.0, 1.0
    t = (m2 - m1) / se
    num = (s1 * s1 / n1 + s2 * s2 / n2) ** 2
    den = (s1 ** 2 / n1) ** 2 / (n1 - 1) + (s2 ** 2 / n2) ** 2 / (n2 - 1)
    df = num / den if den else n1 + n2 - 2
    return m1, m2, t, _betai(df / 2, 0.5, df / (df + t * t))


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = math.sqrt(sum((x - ma) ** 2 for x in a))
    db = math.sqrt(sum((y - mb) ** 2 for y in b))
    if da * db == 0:
        return 0.0, 0.0
    r = num / (da * db)
    t = r * math.sqrt((n - 2) / max(1e-12, 1 - r * r))
    return r, t


def spearman(a, b):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(rank(a), rank(b))


def zscore(v):
    m, sd, _ = mean_sd_se(v)
    sd = sd or 1.0
    return [(x - m) / sd for x in v]


# --------------------------------------------------------------------------- #
# Analyses
# --------------------------------------------------------------------------- #
def global_drift(store: Store, p0, p7, n_probe=4000, seed=0):
    """Genome-wide P7/P0 ratio distribution -- the baseline any single gene must
    beat before it counts as specific.

    P7 cells have modestly lower library complexity, so almost everything drifts
    down. Without this control a 0.9x change looks meaningful when it is not.
    """
    import hashlib
    ratios = []
    step = max(1, len(store.genes) // n_probe)
    # deterministic stride + hash offset: reproducible without Math.random
    off = hashlib.sha256(str(seed).encode()).digest()[0] % step
    for gi in range(off, len(store.genes), step):
        s, e = store._indptr[gi], store._indptr[gi + 1]
        if e - s < 40:
            continue
        v = store.gene(store.genes[gi])
        ma = sum(v[i] for i in p0) / len(p0)
        mb = sum(v[i] for i in p7) / len(p7)
        if ma > 0.5 and mb > 0.0:
            ratios.append(mb / ma)
    ratios.sort()
    return ratios


def quantile(sorted_vals, q):
    if not sorted_vals:
        return float("nan")
    return sorted_vals[int(q * (len(sorted_vals) - 1))]


def percentile_of(sorted_vals, x):
    return 100.0 * sum(1 for r in sorted_vals if r < x) / len(sorted_vals)


def maturation_index(store: Store):
    """M = mean z(fatty-acid oxidation) - mean z(glycolysis), per cell.

    A metabolic surrogate for cardiomyocyte maturation. z-scored across all 285
    cells, so M is a relative index *within this dataset* and has no absolute
    cross-system scale.
    """
    fao = [g for g in FAO if store.has(g)]
    gly = [g for g in GLYCOLYSIS if store.has(g)]
    Z = {g: zscore(store.gene(g)) for g in fao + gly}
    n = len(store.cells)
    return [sum(Z[g][i] for g in fao) / len(fao)
            - sum(Z[g][i] for g in gly) / len(gly) for i in range(n)], fao, gly


def circularity_audit(panel):
    """Split a gene panel by whether it helped define the phase labels."""
    up = {g: g.upper() for g in panel}
    in_s = [g for g in panel if up[g] in PHASE_S]
    in_g2m = [g for g in panel if up[g] in PHASE_G2M]
    clean = [g for g in panel if up[g] not in PHASE_S | PHASE_G2M]
    return {"in_S": in_s, "in_G2M": in_g2m, "clean": clean}


def sort_enrichment(store: Store):
    """The sort is asymmetric: P0 unenriched, P7 heavily enriched.

    In-vivo expectation is the complement of Baniol's mKO2+ (G0/G1) fraction:
    67.5% at P0, >90% at P7. The four FUCCI states are exhaustive (their P0
    numbers sum to 99.96%), which is what licenses taking the complement.
    """
    out = {}
    invivo = {"P0": 1 - 0.675, "P7": (0.100, 0.0873)}  # P7: bound, component-sum
    for st in ("P0", "P7"):
        cyc = len(store.group(stage=st, cycling_score="cycling"))
        tot = len(store.group(stage=st))
        obs = cyc / tot
        exp = invivo[st] if st == "P0" else invivo[st]
        if isinstance(exp, tuple):
            out[st] = {"observed": obs, "expected_range": exp,
                       "enrichment_range": (obs / exp[0], obs / exp[1]),
                       "n_cycling": cyc, "n_total": tot}
        else:
            out[st] = {"observed": obs, "expected": exp,
                       "enrichment": obs / exp, "n_cycling": cyc, "n_total": tot}
    return out


def run(store: Store | None = None) -> dict:
    """Every quoted result, as one dict. Consumed by figures.py and RESULTS.md."""
    store = store or Store()
    g8 = store.groups8()
    p0c, p7c = g8["P0-vCM-cycling"], g8["P7-vCM-cycling"]
    p0n, p7n = g8["P0-vCM-noncycling"], g8["P7-vCM-noncycling"]

    res: dict = {
        "dataset": {**store.info, "dataset": store.dataset},
        "group_sizes": {k: len(v) for k, v in g8.items()},
    }

    # --- 1. Ect2 specificity, against the genome-wide drift baseline ---------
    drift = global_drift(store, p0c, p7c)
    panel = [g for g in CYTOKINESIS + MITOTIC if store.has(g)]
    per_gene = {}
    for g in panel:
        v = store.gene(g)
        m0, m7, t, p = welch([v[i] for i in p0c], [v[i] for i in p7c])
        ratio = (m7 / m0) if m0 > 0 else float("nan")
        per_gene[g] = {
            "P0_cyc": m0, "P7_cyc": m7, "ratio": ratio, "t": t, "p": p,
            "module": "cytokinesis" if g in CYTOKINESIS else "mitotic",
            "pct_genomewide": percentile_of(drift, ratio),
        }
    res["drift"] = {
        "n_genes": len(drift),
        "q05": quantile(drift, .05), "q25": quantile(drift, .25),
        "median": quantile(drift, .50),
        "q75": quantile(drift, .75), "q95": quantile(drift, .95),
    }
    res["per_gene"] = per_gene
    res["bonferroni_alpha"] = 0.05 / len(panel)

    # --- 2. Module scores cancel the signal ---------------------------------
    def module_vec(genes, cells):
        gs = [g for g in genes if store.has(g)]
        V = {g: store.gene(g) for g in gs}
        return [sum(V[g][i] for g in gs) / len(gs) for i in cells]

    mods = {}
    for lab, genes in (("cytokinesis", CYTOKINESIS), ("mitotic", MITOTIC)):
        a, b = module_vec(genes, p0c), module_vec(genes, p7c)
        m0, m7, t, p = welch(a, b)
        mods[lab] = {"P0": m0, "P7": m7, "ratio": m7 / m0, "t": t, "p": p,
                     "P0_se": mean_sd_se(a)[2], "P7_se": mean_sd_se(b)[2]}
    ca = [x - y for x, y in zip(module_vec(CYTOKINESIS, p0c), module_vec(MITOTIC, p0c))]
    cb = [x - y for x, y in zip(module_vec(CYTOKINESIS, p7c), module_vec(MITOTIC, p7c))]
    m0, m7, t, p = welch(ca, cb)
    mods["difference_index"] = {"P0": m0, "P7": m7, "delta": m7 - m0, "t": t, "p": p,
                                "P0_se": mean_sd_se(ca)[2], "P7_se": mean_sd_se(cb)[2]}
    res["modules"] = mods

    # --- 3. Maturation coordinate ------------------------------------------
    M, fao, gly = maturation_index(store)
    res["maturation"] = {
        "fao_genes": fao, "glycolysis_genes": gly,
        "by_group": {k: dict(zip(("mean", "sd", "se"), mean_sd_se([M[i] for i in v])),
                             n=len(v)) for k, v in g8.items()},
    }
    vcm_cyc = p0c + p7c
    corr = {}
    for g in ["Ect2", "E2f6", "Ccne2", "E2f2", "Ccng1"]:
        if not store.has(g):
            continue
        v = store.gene(g)
        r, t = pearson([M[i] for i in vcm_cyc], [v[i] for i in vcm_cyc])
        corr[g] = {"r": r, "t": t, "n": len(vcm_cyc)}
    res["maturation"]["correlations"] = corr
    for st, cells in (("P0", p0c + p0n), ("P7", p7c + p7n)):
        res["maturation"].setdefault("within_stage_vCM", {})[st] = \
            dict(zip(("mean", "sd", "se"), mean_sd_se([M[i] for i in cells])))

    # --- 4. E2F family -----------------------------------------------------
    e2f = {}
    for g in E2F_FAMILY:
        if not store.has(g):
            continue
        v = store.gene(g)
        e2f[g] = {k: {"mean": sum(v[i] for i in cells) / len(cells),
                      "detected_pct": 100 * sum(1 for i in cells if v[i] > 0) / len(cells)}
                  for k, cells in g8.items()}
        m0, m7, t, p = welch([v[i] for i in p0c], [v[i] for i in p7c])
        e2f[g]["maturation_t"] = t
        e2f[g]["maturation_p"] = p
    res["e2f"] = e2f

    # --- 5. Edges tested before wiring them --------------------------------
    edges = [("E2f1", "E2f7"), ("E2f1", "E2f8"), ("E2f7", "E2f8"),
             ("E2f1", "Ccne1"), ("E2f1", "Ccne2"), ("E2f1", "E2f2"),
             ("E2f1", "E2f6"), ("E2f6", "Ect2"), ("E2f8", "Ccna2"),
             ("Gmnn", "Cdt1")]
    res["edges"] = {}
    for a, b in edges:
        if not (store.has(a) and store.has(b)):
            continue
        va, vb = store.gene(a), store.gene(b)
        rho, t = spearman([va[i] for i in vcm_cyc], [vb[i] for i in vcm_cyc])
        res["edges"][f"{a}~{b}"] = {"rho": rho, "t": t, "n": len(vcm_cyc)}

    # --- 6. Circularity, anchors, sort enrichment ---------------------------
    res["circularity"] = circularity_audit(
        E2F_FAMILY + CYTOKINESIS + MITOTIC
        + ["Ccna2", "Ccne1", "Ccne2", "Cdt1", "Gmnn", "Cdkn1a", "Rb1", "Ccng1"]
        + DDR)
    res["sort_enrichment"] = sort_enrichment(store)

    anchors = []
    for node, obvious, better in ANCHOR_CHECKS:
        row = {"node": node}
        for role, g in (("obvious", obvious), ("use", better)):
            if not store.has(g):
                row[role] = {"gene": g, "absent": True}
                continue
            v = store.gene(g)
            m0, m7, t, p = welch([v[i] for i in p0c], [v[i] for i in p7c])
            row[role] = {"gene": g, "P0_cyc": m0, "P7_cyc": m7, "t": t,
                         "det_P0": 100 * sum(1 for i in p0c if v[i] > 0) / len(p0c),
                         "det_P7": 100 * sum(1 for i in p7c if v[i] > 0) / len(p7c)}
        anchors.append(row)
    res["anchors"] = anchors

    allp0, allp7 = store.group(stage="P0"), store.group(stage="P7")
    res["adrenergic"] = {}
    for g in ADRENERGIC:
        if not store.has(g):
            res["adrenergic"][g] = {"absent": True}
            continue
        v = store.gene(g)
        res["adrenergic"][g] = {
            "P0": sum(v[i] for i in allp0) / len(allp0),
            "P7": sum(v[i] for i in allp7) / len(allp7),
            "det_P0": 100 * sum(1 for i in allp0 if v[i] > 0) / len(allp0),
            "det_P7": 100 * sum(1 for i in allp7 if v[i] > 0) / len(allp7),
        }

    res["ddr"] = {}
    for g in DDR:
        if not store.has(g):
            continue
        v = store.gene(g)
        m0, m7, t, p = welch([v[i] for i in p0c], [v[i] for i in p7c])
        res["ddr"][g] = {"P0_cyc": m0, "P7_cyc": m7, "t": t, "p": p}

    return res


def main() -> None:
    r = run()
    d = r["dataset"]
    print(f"{d['dataset']}: {d['n_genes']} genes x {d['n_cells']} cells, nnz={d['nnz']}")
    print(f"groups: {r['group_sizes']}\n")

    dr = r["drift"]
    print(f"genome-wide P7/P0 drift (n={dr['n_genes']}): "
          f"median {dr['median']:.3f}  IQR {dr['q25']:.3f}-{dr['q75']:.3f}")
    print(f"Bonferroni alpha = {r['bonferroni_alpha']:.4f}\n")
    print(f"{'gene':10s}{'module':13s}{'ratio':>7s}{'t':>7s}{'p':>10s}{'pctile':>8s}")
    print("-" * 55)
    for g, v in sorted(r["per_gene"].items(), key=lambda kv: kv[1]["ratio"]):
        print(f"{g:10s}{v['module']:13s}{v['ratio']:7.2f}{v['t']:7.2f}"
              f"{v['p']:10.4f}{v['pct_genomewide']:7.1f}%")

    m = r["modules"]
    print(f"\nmodule means (cycling vCM):")
    for k in ("cytokinesis", "mitotic"):
        print(f"  {k:16s} P0 {m[k]['P0']:.3f}  P7 {m[k]['P7']:.3f}  "
              f"ratio {m[k]['ratio']:.3f}  p={m[k]['p']:.4f}")
    di = m["difference_index"]
    print(f"  {'cyto - mito':16s} P0 {di['P0']:+.3f}  P7 {di['P7']:+.3f}  "
          f"delta {di['delta']:+.3f}  p={di['p']:.4f}  <-- the signal cancels")

    print("\nmaturation index M by group:")
    for k, v in sorted(r["maturation"]["by_group"].items(), key=lambda kv: kv[1]["mean"]):
        print(f"  {k:22s} n={v['n']:3d}  M={v['mean']:+.3f} +/- {v['se']:.3f}")
    print("\nM correlations within cycling vCM (n=%d):" %
          next(iter(r["maturation"]["correlations"].values()))["n"])
    for g, v in r["maturation"]["correlations"].items():
        print(f"  {g:8s} r={v['r']:+.3f}  t={v['t']:+.2f}")

    c = r["circularity"]
    print(f"\ncircularity audit: {len(c['in_S'])} genes in _S_GENES, "
          f"{len(c['in_G2M'])} in _G2M_GENES")
    print(f"  in _S_GENES  : {c['in_S']}")
    print(f"  in _G2M_GENES: {c['in_G2M']}")

    se = r["sort_enrichment"]
    print(f"\nsort enrichment: P0 {se['P0']['observed']*100:.1f}% cycling "
          f"(expected {se['P0']['expected']*100:.1f}%, {se['P0']['enrichment']:.2f}x)")
    er = se["P7"]["enrichment_range"]
    print(f"                 P7 {se['P7']['observed']*100:.1f}% cycling "
          f"({min(er):.1f}-{max(er):.1f}x enriched) <-- fractions unusable as rates")


if __name__ == "__main__":
    main()
