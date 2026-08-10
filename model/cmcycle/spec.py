"""Load, lint and calibrate the ``cmfate`` network from its text spec.

The spec is three tracked, diffable text files rather than a binary workbook.
That is a deliberate choice with a concrete reason: in the sibling `cardiac-models`
repo a bare ``data/`` line in .gitignore silently swallowed every Excel network
spec, and a missing `.xls` is indistinguishable from a present one in a diff. A
curated network is reviewed one edge at a time, so the diff unit must be one edge
-- hence one reaction per row, each carrying its own ``evidence`` column.

The fate layer
--------------
The four fates are mutually exclusive and sum to exactly 1, with no normalization
step and no engine change, via a complementary product over three gate nodes:

    !SPhase                              => Quiescent
    SPhase & MitoticEntry & Abscission   => Division
    SPhase & MitoticEntry & !Abscission  => Binucleation
    SPhase & !MitoticEntry               => Polyploidization

Writing g() for the shared activation function, at steady state

    Q = 1 - g(S)
    D = g(S)*g(M)*g(A)
    B = g(S)*g(M)*(1 - g(A))
    P = g(S)*(1 - g(M))          =>  Q + D + B + P = 1  identically.

Three structural invariants make that hold, and :func:`lint` enforces all three:
exactly one reaction per fate node (a second OR'd reaction breaks the identity),
``Weight = 1`` on all four, and identical ``(n, EC50)`` across all four -- because
those are per-reaction rather than per-reactant, so the Division and Binucleation
rules must apply the *same* g to MitoticEntry.

The identity is exact only at steady state; during a transient the sum deviates
by the differing tau. That deviation is informative, so normalization belongs in a
reporting helper and never in the model.

Calibration
-----------
Because the three gates are independent by construction -- competence is not
gated on prevalence -- inverting Murganti's measured fate fractions gives three
*separate* one-dimensional solves rather than a joint fit:

    g(SPhase) = 1 - Q = 0.0965      -> SPhase activity 0.0992
    g(MitoticEntry) = (D+B)/g(S)    -> 0.7698
    g(Abscission) = D/(D+B)         -> 0.5028

So the fate layer is fitted *exactly* and predicts nothing. All predictive content
sits in whether the upstream network reproduces the OTHER contexts from that one
calibration -- which is what :func:`cmcycle.model.clonidine_triad` tests.
"""
from __future__ import annotations

import csv
import math
import re
from importlib.resources import files

from .logic import Network, Reaction, ec50_ceiling

_DATA = files(__package__).joinpath("data")

_ID_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
MAX_DEPTH = 7          # hops from any input to any fate node


def parse_rule(rule: str):
    """``'C & !D => E'`` -> ``('E', [('C', False), ('D', True)], False)``.

    An empty left side marks an input reaction, whose weight *is* the stimulus.
    """
    if "=>" not in rule:
        raise ValueError(f"rule has no '=>': {rule!r}")
    left, right = rule.split("=>")
    target = right.strip()
    left = left.strip()
    if not target:
        raise ValueError(f"rule has no target: {rule!r}")
    if not left:
        return target, [], True
    reactants = []
    for tok in left.split("&"):
        tok = tok.strip()
        if not tok:
            raise ValueError(f"empty reactant in {rule!r}")
        reactants.append((tok.lstrip("!").strip(), tok.startswith("!")))
    return target, reactants, False


def _read(name):
    with _DATA.joinpath(name).open("r", encoding="utf8") as fh:
        return list(csv.DictReader(fh))


def _toml(name):
    import tomllib
    with _DATA.joinpath(name).open("rb") as fh:
        return tomllib.load(fh)


def load(species_csv="cmfate_species.csv", reactions_csv="cmfate_reactions.csv",
         manifest="cmfate_model.toml") -> Network:
    """Parse the spec into a :class:`~cmcycle.logic.Network`.

    Strict on purpose, unlike the lenient Excel loader it replaces: an unknown
    node id in a hand-authored spec is a typo, and silently dropping the reaction
    produces a plausible-looking but wrong network.
    """
    sp = [r for r in _read(species_csv) if (r.get("ID") or "").strip()]
    rx = [r for r in _read(reactions_csv) if (r.get("Rule") or "").strip()]
    man = _toml(manifest)

    species, names, yinit, ymax, tau, modules, genes, orient, notes = ([] for _ in range(9))
    for r in sp:
        sid = r["ID"].strip()
        if not _ID_RE.match(sid):
            raise ValueError(f"bad species id {sid!r}")
        if sid in species:
            raise ValueError(f"duplicate species id {sid!r}")
        species.append(sid)
        names.append((r.get("name") or sid).strip())
        yinit.append(float(r.get("Yinit") or 0))
        ymax.append(float(r.get("Ymax") or 1))
        tau.append(float(r.get("tau") or 1))
        modules.append((r.get("module") or "").strip())
        genes.append([g for g in (r.get("genes") or "").split(";") if g])
        orient.append((r.get("data_orient") or "").strip())
        notes.append((r.get("notes") or "").strip())
    idx = {s: i for i, s in enumerate(species)}

    reactions = []
    seen = set()
    for r in rx:
        rid = (r.get("ID") or "").strip()
        if rid in seen:
            raise ValueError(f"duplicate reaction id {rid!r}")
        seen.add(rid)
        target, reactants, is_input = parse_rule(r["Rule"])
        if target not in idx:
            raise ValueError(f"{rid}: unknown target {target!r}")
        for s, _ in reactants:
            if s not in idx:
                raise ValueError(f"{rid}: unknown reactant {s!r}")
        for col in ("Weight", "n", "EC50"):
            if (r.get(col) or "").strip() == "":
                raise ValueError(f"{rid}: {col} is required (a blank is not 0)")
        reactions.append(Reaction(
            target=idx[target],
            reactants=[(idx[s], inh) for s, inh in reactants],
            w=float(r["Weight"]), n=float(r["n"]), ec50=float(r["EC50"]),
            is_input=is_input, rid=rid))

    net = Network(species=species, names=dict(zip(species, names)), yinit=yinit,
                  ymax=ymax, tau=tau, reactions=reactions, idx=idx)
    net.meta = dict(man)
    net.meta["modules"] = dict(zip(species, modules))
    net.meta["genes"] = dict(zip(species, genes))
    net.meta["data_orient"] = dict(zip(species, orient))
    net.meta["notes"] = dict(zip(species, notes))
    net.meta["rid"] = {r.rid: j for j, r in enumerate(reactions)}
    return net


# --------------------------------------------------------------------------- #
# lint
# --------------------------------------------------------------------------- #
def depth_from_inputs(net: Network) -> dict:
    """Shortest hop count from any input node to each node."""
    inputs = set(net.input_nodes())
    dist = {s: (0 if s in inputs else math.inf) for s in net.species}
    for _ in range(len(net.species)):
        changed = False
        for r in net.reactions:
            if r.is_input or not r.reactants:
                continue
            tgt = net.species[r.target]
            src = min(dist[net.species[i]] for i, _ in r.reactants)
            if src + 1 < dist[tgt]:
                dist[tgt] = src + 1
                changed = True
        if not changed:
            break
    return dist


def lint(net: Network) -> list:
    """Structural checks. Returns a list of problem strings; empty means clean."""
    p = []
    fates = net.meta.get("fate_nodes", [])
    gates = net.meta.get("gate_nodes", [])

    # -- the three invariants that make the fate partition exact --------------
    fate_rx = {f: [r for r in net.reactions if net.species[r.target] == f] for f in fates}
    for f, rs in fate_rx.items():
        if len(rs) != 1:
            p.append(f"fate {f}: {len(rs)} reactions, must be exactly 1 "
                     f"(a second OR'd reaction breaks Q+D+B+P=1)")
    ws = {f: rs[0].w for f, rs in fate_rx.items() if len(rs) == 1}
    if any(abs(w - 1.0) > 1e-12 for w in ws.values()):
        p.append(f"fate weights must all be 1.0, got {ws}")
    shapes = {(rs[0].n, rs[0].ec50) for rs in fate_rx.values() if len(rs) == 1}
    if len(shapes) > 1:
        p.append(f"fate reactions must share one (n, EC50); got {sorted(shapes)}")
    for f in fates:
        if net.ymax[net.idx[f]] != 1.0:
            p.append(f"fate {f}: Ymax must be 1 so activity == probability")

    # -- Hill validity (the engine raises, so this is a spec-level restatement)
    for r in net.reactions:
        if r.ec50 ** r.n >= 0.5:
            p.append(f"{r.rid}: EC50**n >= 0.5 (ceiling {ec50_ceiling(r.n):.4f})")

    # -- every non-input, non-integrator node should be anchorable to a gene ---
    integrators = set(gates) | set(fates) | {
        "E2Frep", "E2Fact", "ROS", "MitCompRaw", "AbsRaw"}
    for s in net.species:
        mod = net.meta["modules"].get(s)
        if mod in ("input", "fate") or s in integrators:
            continue
        if not net.meta["genes"].get(s):
            p.append(f"{s}: no gene anchor and not a declared integrator")

    # -- reachability and depth budget ---------------------------------------
    dist = depth_from_inputs(net)
    for f in fates:
        if math.isinf(dist[f]):
            p.append(f"fate {f} is unreachable from any input")
        elif dist[f] > MAX_DEPTH:
            p.append(f"fate {f} is {dist[f]} hops from an input (budget {MAX_DEPTH})")
    for s in net.species:
        if math.isinf(dist[s]) and net.meta["modules"].get(s) != "input":
            p.append(f"{s} is unreachable from any input")

    # -- every node should be read by something, or be a declared output ------
    # A node that drives nothing and is not reported is dead weight: it inflates
    # the node count without touching a prediction. Declared outputs are exempt
    # because being read by the analyst is the job (Ki67 and Cdc20 are the
    # instantaneous reporters the two papers actually stain for).
    used = {i for r in net.reactions for i, _ in r.reactants}
    reported = set(fates) | set(net.meta.get("outputs", []))
    for i, s in enumerate(net.species):
        if i not in used and s not in reported:
            p.append(f"{s} drives nothing and is not a declared output")
    return p


# --------------------------------------------------------------------------- #
# calibration
# --------------------------------------------------------------------------- #
def calibrate(net: Network, verbose=False) -> dict:
    """Solve each gate knob so the gates hit their target at the fit context.

    One bisection per gate. The gates are independent by construction, so no
    joint fit and no iteration are needed -- which is the practical payoff of
    keeping competence ungated on prevalence.
    """
    cal = net.meta["calibration"]
    ctx = net.meta["contexts"][cal["context"]]["inputs"]
    out = {}
    for node, target in cal["targets"].items():
        rid = cal["knobs"][node]
        j = net.meta["rid"][rid]
        r = net.reactions[j]
        lo, hi = 0.0, 1.0

        def at(w):
            saved, r.w = r.w, w
            try:
                return net.ss(inputs=ctx)[node]
            finally:
                r.w = saved

        if at(hi) < target:
            raise RuntimeError(
                f"{node}: even w=1 on {rid} only reaches {at(hi):.4f} < {target}. "
                f"The upstream drive is too weak - raise an activator, not this knob.")
        for _ in range(32):
            mid = 0.5 * (lo + hi)
            if at(mid) < target:
                lo = mid
            else:
                hi = mid
        r.w = 0.5 * (lo + hi)
        out[node] = {"reaction": rid, "weight": r.w, "target": target,
                     "achieved": net.ss(inputs=ctx)[node]}
        if verbose:
            print(f"  {node:14s} {rid}  w={r.w:.6f}  "
                  f"target {target:.4f} -> {out[node]['achieved']:.4f}")
    return out


CACHE = "cmfate_calibration.json"


def _spec_fingerprint() -> str:
    """Hash of the three spec files, so a cached fit can never outlive its spec."""
    import hashlib
    h = hashlib.sha256()
    for name in ("cmfate_species.csv", "cmfate_reactions.csv", "cmfate_model.toml"):
        h.update(_DATA.joinpath(name).read_bytes())
    return h.hexdigest()[:16]


def load_calibrated(verbose=False, use_cache=True):
    """The network as everything downstream should use it.

    The three bisections take a while in pure Python, so the solved weights are
    cached beside the spec and keyed on a hash of the spec files. Edit any spec
    file and the cache is ignored automatically -- a stale fit silently applied to
    a changed network is exactly the failure worth designing out.
    """
    import json

    net = load()
    problems = lint(net)
    if problems:
        raise ValueError("spec lint failed:\n  " + "\n  ".join(problems))

    fp = _spec_fingerprint()
    cached = None
    if use_cache:
        try:
            blob = json.loads(_DATA.joinpath(CACHE).read_text())
            if blob.get("fingerprint") == fp:
                cached = blob["calibrated"]
        except (FileNotFoundError, ValueError, KeyError):
            cached = None

    if cached:
        for node, rec in cached.items():
            net.reactions[net.meta["rid"][rec["reaction"]]].w = rec["weight"]
        if verbose:
            print(f"  using cached calibration (spec {fp})")
        net.meta["calibrated"] = cached
    else:
        cal = calibrate(net, verbose=verbose)
        net.meta["calibrated"] = cal
        try:
            import os
            path = os.path.join(os.path.dirname(str(_DATA.joinpath(CACHE))), CACHE)
            with open(path, "w", encoding="utf8") as fh:
                json.dump({"fingerprint": fp, "calibrated": cal}, fh, indent=1)
            if verbose:
                print(f"  wrote {CACHE} (spec {fp})")
        except OSError:
            pass          # read-only install: recalibrating each time is correct
    return net
