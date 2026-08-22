"""Normalized-Hill logic-ODE engine, dependency-free.

Implements the Netflux formalism (Kraeutler, Soltis & Saucerman 2010, BMC Syst
Biol 4:157) in pure Python, with an embedded adaptive Runge-Kutta integrator so
the package keeps its no-dependency property.

Each species activity Y is a normalized fractional activation in [0, 1]:

    act(x)      = w * beta * x^n / (K^n + x^n),   clamped to [0, w]
    beta        = (EC50^n - 1) / (2*EC50^n - 1)
    K           = (beta - 1)^(1/n)
    inhib(x)    = w - act(x)
    AND(a, b)   = a * b / w^(k-1)      (k reactants; renormalized to [0, w])
    OR(a, b)    = a + b - a*b          (reactions onto the same node)
    dY/dt       = (drive * Ymax - Y) / tau

Input reactions ("=> A") contribute their weight directly, so the weight of an
input reaction *is* the stimulus level.

Two corrections relative to the reference implementation, both load-bearing here
because this model deliberately uses high-threshold gates:

1. ``EC50 ** n < 0.5`` is a hard validity constraint. Above it ``beta`` turns
   NEGATIVE -- at EC50=0.7, n=1.4 it is -1.838, so ``act`` returns a negative
   multiple of w for *every* input. The reference clamps only the upper bound, so
   the sign error propagates through ``w - act(x) > w`` and the weighted-OR to
   push activities above 1, silently. :class:`Reaction` now raises instead, and
   ``_act`` clamps at zero as well.
2. ``perturbation_matrix`` accepts a tuple of node ids per column, so a *double*
   knockdown is expressible. The reference takes one node at a time, which cannot
   represent the E2f7/E2f8 double knockout this project is built around.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

#: Hard ceiling on EC50 for a given Hill coefficient: EC50 < 2**(-1/n).
def ec50_ceiling(n: float) -> float:
    return 2.0 ** (-1.0 / n)


@dataclass
class Reaction:
    """One reaction. ``reactants`` is a list of (species_index, is_inhibitor)."""
    target: int
    reactants: list
    w: float
    n: float
    ec50: float
    is_input: bool = False
    rid: str = ""
    beta: float = field(default=0.0, init=False)
    K: float = field(default=0.0, init=False)

    def __post_init__(self):
        n, ec = self.n, self.ec50
        if not (0.0 < ec < 1.0):
            raise ValueError(f"{self.rid or 'reaction'}: EC50 must be in (0,1), got {ec}")
        if n <= 0:
            raise ValueError(f"{self.rid or 'reaction'}: n must be > 0, got {n}")
        if ec ** n >= 0.5:
            raise ValueError(
                f"{self.rid or 'reaction'}: EC50**n = {ec ** n:.4f} >= 0.5 makes beta "
                f"negative, which inverts act() for every input. For n={n}, EC50 must "
                f"be < {ec50_ceiling(n):.4f}. A high-threshold gate needs a larger n."
            )
        denom = 2 * ec ** n - 1
        self.beta = (ec ** n - 1) / denom
        base = self.beta - 1
        self.K = base ** (1.0 / n) if base > 0 else 1e-6


@dataclass
class Network:
    species: list
    names: dict = field(default_factory=dict)
    yinit: list = field(default_factory=list)
    ymax: list = field(default_factory=list)
    tau: list = field(default_factory=list)
    reactions: list = field(default_factory=list)
    idx: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    # ---- helpers ---------------------------------------------------------
    def input_nodes(self) -> list:
        return [self.species[r.target] for r in self.reactions if r.is_input]

    def _act(self, x: float, r: Reaction) -> float:
        if r.w == 0:
            return 0.0
        x = 0.0 if x < 0 else x
        f = r.w * (r.beta * x ** r.n) / (r.K ** r.n + x ** r.n)
        if f < 0.0:            # unreachable for a validated reaction; belt and braces
            return 0.0
        return r.w if f > r.w else f

    def _reaction_value(self, r: Reaction, y: list, w: float) -> float:
        if r.is_input or not r.reactants:
            return w
        vals = []
        for sidx, is_inhib in r.reactants:
            a = self._act(y[sidx], r)
            vals.append((r.w - a) if is_inhib else a)
        if len(vals) == 1:
            return vals[0]
        prod = 1.0
        for v in vals:
            prod *= v
        return 0.0 if r.w == 0 else prod / (r.w ** (len(vals) - 1))

    def _compiled(self):
        """Flatten reactions to plain tuples once.

        The RHS is evaluated tens of thousands of times per screen, and in pure
        Python the dataclass attribute lookups in the inner loop dominated the
        run time, so the flattened form is cached on the instance.

        The cache key covers **both** the weights and the identity of every
        reaction object. Keying on weights alone was a correctness bug: replacing a
        reaction to change its reactants, ``n`` or ``EC50`` left the key unchanged,
        so the edit was silently ignored and the model kept running the old
        structure. That is exactly what programmatic spec edits do, and it produced
        a diagnosis that had to be thrown away.

        Note the contract this implies: a reaction's ``w`` may be mutated in place
        (calibration and the perturbation sweeps do), but anything else requires
        **replacing** the :class:`Reaction` — which is necessary anyway, since
        ``beta`` and ``K`` are derived in ``__post_init__`` and would otherwise go
        stale.
        """
        key = (tuple(r.w for r in self.reactions),
               tuple(id(r) for r in self.reactions))
        cache = getattr(self, "_cc", None)
        if cache is not None and cache[0] == key:
            return cache[1]
        flat = [(r.target, tuple(r.reactants), r.w, r.n, r.K ** r.n, r.beta, r.is_input)
                for r in self.reactions]
        self._cc = (key, flat)
        return flat

    def _drive(self, y: list, input_w: list) -> list:
        drive = [0.0] * len(self.species)
        for j, (tgt, reacts, w, n, Kn, beta, is_input) in enumerate(self._compiled()):
            if is_input:
                val = input_w[j]
            elif not reacts:
                val = w
            else:
                prod, k = 1.0, 0
                for sidx, inh in reacts:
                    x = y[sidx]
                    if x < 0.0:
                        x = 0.0
                    if w == 0.0:
                        a = 0.0
                    else:
                        xn = x ** n
                        a = w * (beta * xn) / (Kn + xn)
                        if a > w:
                            a = w
                        elif a < 0.0:
                            a = 0.0
                    prod *= (w - a) if inh else a
                    k += 1
                val = 0.0 if w == 0.0 else (prod if k == 1 else prod / (w ** (k - 1)))
            d = drive[tgt]
            drive[tgt] = d + val - d * val                  # weighted OR
        return drive

    def _rhs(self, y: list, ymax: list, input_w: list) -> list:
        drive = self._drive(y, input_w)
        return [(drive[i] * ymax[i] - y[i]) / self.tau[i]
                for i in range(len(self.species))]

    def sink_indices(self) -> set:
        """Nodes that drive nothing -- pure outputs (the four fates, Ki67, Cdc20).

        Their steady state is algebraic (``y = drive * ymax``), so they need not be
        integrated. That matters a lot here: the fate nodes carry tau up to 35.3 h,
        set from the measured S/G2/M durations, and integrating them to convergence
        dominated the run time while changing no steady-state value. Their tau
        therefore shapes transients only, never a steady-state prediction.
        """
        read = {i for r in self.reactions for i, _ in r.reactants}
        return set(range(len(self.species))) - read

    # ---- integration -----------------------------------------------------
    def _setup(self, inputs, knockdowns, overexpress, ymax_scale):
        ymax = list(self.ymax)
        for sid in (knockdowns or []):
            if sid in self.idx:
                ymax[self.idx[sid]] = 0.0
        for sid, f in (ymax_scale or {}).items():
            if sid in self.idx:
                ymax[self.idx[sid]] *= float(f)
        input_w = [r.w for r in self.reactions]
        for j, r in enumerate(self.reactions):
            if r.is_input:
                sid = self.species[r.target]
                if inputs and sid in inputs:
                    input_w[j] = float(inputs[sid])
        y0 = list(self.yinit)
        over = set(overexpress or [])
        for sid in over:
            if sid in self.idx:
                y0[self.idx[sid]] = 1.0
        return ymax, input_w, y0, {self.idx[s] for s in over if s in self.idx}

    def simulate(self, inputs=None, knockdowns=None, overexpress=None,
                 ymax_scale=None, t_end=60.0, n_points=200, rtol=1e-6):
        """Integrate to ``t_end``. Returns (times, Y) with Y[k] a state list.

        Adaptive Bogacki-Shampine (RK23) with a 3rd-order error estimate. The
        network is non-stiff by construction -- every node relaxes on its own tau
        with no fast/slow chemistry -- so an explicit method is appropriate and
        avoids a scipy dependency.
        """
        ymax, input_w, y, held = self._setup(inputs, knockdowns, overexpress, ymax_scale)
        ts = [t_end * k / (n_points - 1) for k in range(n_points)] if n_points > 1 else [t_end]
        out, t, h = [], 0.0, min(0.05, t_end / 50.0)
        nxt = 0

        def f(yy):
            d = self._rhs(yy, ymax, input_w)
            for i in held:
                d[i] = 0.0
            return d

        while nxt < len(ts):
            while nxt < len(ts) and ts[nxt] <= t + 1e-12:
                out.append(list(y)); nxt += 1
            if nxt >= len(ts):
                break
            h = min(h, ts[nxt] - t)
            k1 = f(y)
            y2 = [y[i] + 0.5 * h * k1[i] for i in range(len(y))]
            k2 = f(y2)
            y3 = [y[i] + 0.75 * h * k2[i] for i in range(len(y))]
            k3 = f(y3)
            yn = [y[i] + h * (2 * k1[i] + 3 * k2[i] + 4 * k3[i]) / 9 for i in range(len(y))]
            k4 = f(yn)
            err = max(abs(h * (-5 * k1[i] / 72 + k2[i] / 12 + k3[i] / 9 - k4[i] / 8))
                      for i in range(len(y)))
            tol = rtol * max(1.0, max(abs(v) for v in yn))
            if err <= tol or h <= 1e-9:
                t += h
                y = [0.0 if v < 0.0 else v for v in yn]
                h = min(h * (2.0 if err < tol / 8 else 1.2), t_end / 10.0)
            else:
                h *= 0.5
        return ts, out

    def steady_state(self, inputs=None, knockdowns=None, overexpress=None,
                     ymax_scale=None, t_max=None, tol=1e-7):
        """Integrate until every derivative is below ``tol``.

        ``t_max`` defaults to 40x the slowest tau in the network. That matters
        here: the fate nodes carry tau up to 35.3 h (set from the measured S/G2/M
        durations), so a fixed horizon that suits the signalling layer leaves them
        still relaxing and the run looks like a non-convergence.

        Returns the state list. Raises rather than returning a half-settled state,
        because a silent non-convergence would be read as a fate prediction.
        """
        ymax, input_w, y, held = self._setup(inputs, knockdowns, overexpress, ymax_scale)
        sinks = self.sink_indices() - held
        driver_tau = [self.tau[i] for i in range(len(self.species)) if i not in sinks]
        if t_max is None:
            t_max = max(200.0, 40.0 * max(driver_tau or self.tau))

        def finish(yy):
            """Assign the algebraic steady value to every output-only node."""
            drive = self._drive(yy, input_w)
            for i in sinks:
                yy[i] = drive[i] * ymax[i]
            return yy

        def f(yy):
            d = self._rhs(yy, ymax, input_w)
            for i in held:
                d[i] = 0.0
            for i in sinks:          # excluded from the convergence test
                d[i] = 0.0
            return d

        step_tol = 1e-8
        t, h = 0.0, 0.05
        while t < t_max:
            k1 = f(y)
            if max(abs(v) for v in k1) < tol:
                return finish(y)
            y2 = [y[i] + 0.5 * h * k1[i] for i in range(len(y))]
            k2 = f(y2)
            y3 = [y[i] + 0.75 * h * k2[i] for i in range(len(y))]
            k3 = f(y3)
            yn = [y[i] + h * (2 * k1[i] + 3 * k2[i] + 4 * k3[i]) / 9 for i in range(len(y))]
            k4 = f(yn)
            err = max(abs(h * (-5 * k1[i] / 72 + k2[i] / 12 + k3[i] / 9 - k4[i] / 8))
                      for i in range(len(y)))
            # Step acceptance stays tight (1e-8) even when the caller's convergence
            # tolerance is loose: with large steps the derivative estimate jitters
            # around the fixed point and the convergence test never fires. The max
            # step size, not the error tolerance, is the safe speed knob here.
            if err <= step_tol or h <= 1e-9:
                t += h
                y = [0.0 if v < 0.0 else v for v in yn]
                h = min(h * 1.7, 4.0)
            else:
                h *= 0.5
        d = f(y)
        worst = max(range(len(d)), key=lambda i: abs(d[i]))
        raise RuntimeError(
            f"steady state not reached by t={t_max:.0f}; worst node "
            f"{self.species[worst]} has |dy/dt| = {abs(d[worst]):.2e} (tol {tol:.0e}). "
            f"Either a loop is oscillating, or tau[{self.species[worst]}] = "
            f"{self.tau[worst]:g} needs a longer horizon."
        )

    def ss(self, **kw) -> dict:
        """Steady state as {species_id: activity} -- the convenient form."""
        y = self.steady_state(**kw)
        return dict(zip(self.species, y))

    def perturbation_matrix(self, perturb_nodes=None, output_nodes=None, inputs=None,
                           mode="knockdown", background_knockdowns=None):
        """Delta steady state of each output under each perturbation.

        ``perturb_nodes`` elements may be a node id **or a tuple of ids**, so a
        double knockdown is one column. ``background_knockdowns`` applies to the
        baseline and every perturbed run, for "what adds to X" screens.

        Returns (matrix[len(output)][len(perturb)], output_ids, column_labels).
        """
        perturb_nodes = list(perturb_nodes or self.species)
        output_nodes = list(output_nodes or self.species)
        bg = list(background_knockdowns or [])
        base = self.ss(inputs=inputs, knockdowns=bg)
        mat = [[0.0] * len(perturb_nodes) for _ in output_nodes]
        labels = []
        for c, p in enumerate(perturb_nodes):
            nodes = [p] if isinstance(p, str) else list(p)
            labels.append("+".join(nodes))
            kw = ({"knockdowns": bg + nodes} if mode == "knockdown"
                  else {"knockdowns": bg, "overexpress": nodes})
            pert = self.ss(inputs=inputs, **kw)
            for r_, o in enumerate(output_nodes):
                mat[r_][c] = pert[o] - base[o]
        return mat, output_nodes, labels

    # ---- fate layer ------------------------------------------------------
    def fates(self, inputs=None, **kw) -> dict:
        """The four fate activities, which sum to 1 at steady state.

        See :mod:`cmcycle.spec` for why the complementary-product form makes this
        exact without any normalization step.
        """
        s = self.ss(inputs=inputs, **kw)
        return {k: s[k] for k in self.meta.get("fate_nodes", []) if k in s}
