"""Rigorous equilibrium-stage columns: the MESH equations solved simultaneously by Newton's
method (Naphtali-Sandholm).

Stages are numbered 0..N-1 from the top. A condenser, if present, is stage 0 and a reboiler
is stage N-1. The unknowns of stage j are the liquid and vapor component flows leaving it
(l_ij, v_ij, mol/s) and its temperature T_j. With a total condenser no vapor leaves stage 0;
the distillate is a liquid draw s * l_0 (s = D / reflux). Three global unknowns follow the
stages: s, the condenser duty Qc and the reboiler duty Qr.

Equations of stage j (x = l / L, y = v / V):

    material      (1 + s_j) l_ij + v_ij - l_i,j-1 - v_i,j+1 - f_ij = 0
    equilibrium   y_ij - eta_j K_ij x_ij - (1 - eta_j) y_i,j+1 = 0     (Murphree, eta = 1: ideal)
    energy        (1 + s_j) L_j h_j + V_j H_j - L_j-1 h_j-1 - V_j+1 H_j+1 - HF_j - Q_j = 0

where s_j = s on a total condenser and 0 elsewhere. A total condenser replaces its equilibrium
equations by v_i0 = 0 and adds its bubble-point condition sum(K x) = 1 as the first global
equation (otherwise s = 0).

Side draws, pumparounds, stage duties and reactions extend these equations without new
unknowns:

- A side draw or a pumparound withdraws a ratio r of the flow continuing from its stage
  (w_ij = r l_ij, or r v_ij for a vapor draw). r is constant for a draw given as a fraction of
  the stage's outflow and W / L_j for one given as a flow, so it depends on stage j alone.
- A pumparound returns its liquid r l_a to stage b as a feed with enthalpy r L_a h_a + Q, or
  with the enthalpy of x_a at its return temperature (a TP flash), and so couples stage a to a
  stage that may be far away; its source stage gets its own finite-difference evaluation.
- A stage duty adds to Q_j.
- A kinetic reactive stage with liquid hold-up V_j forms g_ij = V_j sum_r nu_ri r_r(T_j, x_j)
  of component i, and its energy balance gains sum_i g_ij hf_i: the enthalpies exclude the
  formation enthalpies hf, so this is the heat of reaction at the stage (as in the reactors). Each existing duty pairs with one specification (reflux ratio,
distillate rate, bottoms rate or boilup ratio); a missing condenser or reboiler fixes its
duty at zero. The Jacobian is block-tridiagonal apart from the global rows and columns; it is
built by finite differences with column coloring, so each Newton step costs about four
property evaluations per stage and unknown.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.linalg import solve_banded
from scipy.sparse.linalg import spsolve

from ..errors import SimulationError
from ..thermo_package import PropertyPackage
from .saturation import bubble_point, dew_point

SPECS = {
    "reflux_ratio": "reflux ratio",
    "distillate_rate": "distillate rate",
    "bottoms_rate": "bottoms rate",
    "boilup_ratio": "boilup ratio",
}
_H_SCALE = 1.0e4  # J/mol, typical enthalpy difference; scales the energy equations
_MAX_DT = 25.0  # K, largest temperature change per Newton step
_MAX_GROWTH = 10.0  # residual growth accepted for a full Newton step


@dataclass
class StageFeed:
    stage: int  # 0-based, from the top
    flows: np.ndarray  # mol/s per component
    H: float  # enthalpy flow, W
    vapor_flows: np.ndarray  # the part of `flows` that is vapor
    T: float  # K, used for the initial temperature profile of columns without condenser/reboiler


@dataclass
class SideDraw:
    stage: int  # 0-based
    phase: str  # "liquid" or "vapor"
    rate: float | None = None  # mol/s
    fraction: float | None = None  # of the stage's total liquid (vapor) outflow


@dataclass
class Pumparound:
    """Liquid drawn from stage `source`, heated or cooled, returned to stage `target`."""

    source: int  # 0-based
    target: int
    rate: float | None = None  # mol/s
    fraction: float | None = None  # of the source stage's total liquid outflow
    duty: float | None = None  # W, heat added; zero if neither duty nor T is given
    T: float | None = None  # K, return temperature


@dataclass
class ColumnSetup:
    stages: int
    P: np.ndarray  # Pa, per stage
    feeds: list[StageFeed]
    condenser: str = "total"  # "total", "partial" or "none"
    reboiler: bool = True
    specs: dict[str, float] = field(default_factory=dict)
    efficiency: float = 1.0  # Murphree vapor efficiency of stages other than condenser/reboiler
    draws: list[SideDraw] = field(default_factory=list)
    pumparounds: list[Pumparound] = field(default_factory=list)
    stage_duties: dict[int, float] = field(default_factory=dict)  # stage -> W added
    # Kinetic liquid-phase reactions (reactions.Reaction) on stages with a liquid hold-up.
    reactions: list = field(default_factory=list)
    holdup: dict[int, float] = field(default_factory=dict)  # stage -> m3 of liquid
    formation_enthalpies: np.ndarray | None = None  # J/mol per component (0 if unknown)
    max_iterations: int = 60
    tolerance: float = 1e-10  # on the scaled residuals


@dataclass
class ColumnSolution:
    T: np.ndarray
    P: np.ndarray
    L: np.ndarray  # liquid flowing down from each stage (reflux on a condenser), mol/s
    V: np.ndarray  # vapor flowing up from each stage, mol/s
    x: np.ndarray  # liquid mole fractions, stages x components
    y: np.ndarray  # vapor mole fractions (equilibrium vapor on a total condenser)
    h: np.ndarray  # liquid molar enthalpy per stage, J/mol
    Hv: np.ndarray  # vapor molar enthalpy per stage, J/mol
    distillate: np.ndarray  # liquid distillate component flows (total condenser)
    vapor_distillate: np.ndarray  # vapor leaving the top stage
    bottoms: np.ndarray
    condenser_duty: float  # W, heat added (negative: removed)
    reboiler_duty: float
    iterations: int
    side_draws: list[np.ndarray] = field(default_factory=list)  # component flows per draw
    pumparound_flows: list[np.ndarray] = field(default_factory=list)
    pumparound_duties: list[float] = field(default_factory=list)  # W
    pumparound_return_H: list[float] = field(default_factory=list)  # W, enthalpy flow returned
    generation: np.ndarray | None = None  # stages x components, mol/s formed by reactions
    # The converged unknowns, for a warm start of the same column (see solve_column).
    restart: tuple | None = field(default=None, repr=False)


class _Model:
    """Unknown layout, stage properties and residuals for one column setup."""

    def __init__(self, pkg: PropertyPackage, setup: ColumnSetup):
        self.pkg, self.setup = pkg, setup
        self.N, self.C = N, C = setup.stages, pkg.n
        self.nv = 2 * C + 1
        self.size = N * self.nv + 3
        self.P = np.asarray(setup.P, dtype=float)
        self.total = setup.condenser == "total"
        self.has_cond = setup.condenser != "none"
        self.has_reb = setup.reboiler
        self.f = np.zeros((N, C))
        self.HF = np.zeros(N)
        self.fv = np.zeros((N, C))
        for feed in setup.feeds:
            self.f[feed.stage] += feed.flows
            self.HF[feed.stage] += feed.H
            self.fv[feed.stage] += feed.vapor_flows
        self.F = self.f.sum()
        self.E = self.F * _H_SCALE
        self.eta = np.full(N, setup.efficiency)
        if self.has_cond:
            self.eta[0] = 1.0
        if self.has_reb:
            self.eta[-1] = 1.0
        # Vapor entering the bottom stage from below: only a vapor feed (Murphree needs it).
        fv_bottom = self.fv[-1].sum()
        self.y_bottom_in = self.fv[-1] / fv_bottom if fv_bottom > 0 else None
        if self.y_bottom_in is None:
            self.eta[-1] = 1.0
        # Withdrawals per stage: (key, rate, fraction), key ("draw", k) or ("pa", k).
        self.draws, self.pas = list(setup.draws), list(setup.pumparounds)
        self.liq_w: list[list] = [[] for _ in range(N)]
        self.vap_w: list[list] = [[] for _ in range(N)]
        for k, d in enumerate(self.draws):
            (self.liq_w if d.phase == "liquid" else self.vap_w)[d.stage].append(
                (("draw", k), d.rate, d.fraction))
        for k, pa in enumerate(self.pas):
            self.liq_w[pa.source].append((("pa", k), pa.rate, pa.fraction))
        self.withdrawing = [j for j in range(N) if self.liq_w[j] or self.vap_w[j]]
        self.pa_sources = sorted({pa.source for pa in self.pas})
        self.Qs = np.zeros(N)
        for j, q in setup.stage_duties.items():
            self.Qs[j] += q
        self.holdup = np.zeros(N)
        for j, vol in setup.holdup.items():
            self.holdup[j] += vol
        self.reactions = list(setup.reactions) if self.holdup.any() else []
        self.nu = np.array([r.nu for r in self.reactions]) if self.reactions else None
        if not self.reactions:
            self.holdup[:] = 0.0
        self.hf = (np.zeros(C) if setup.formation_enthalpies is None
                   else np.asarray(setup.formation_enthalpies, dtype=float))
        self.spec_rows = []  # (global row index, spec name, value)
        duty_rows = [1] * self.has_cond + [2] * self.has_reb
        for row, (name, value) in zip(duty_rows, setup.specs.items()):
            self.spec_rows.append((row, name, value))

    # ---- layout --------------------------------------------------------------------------

    def unpack(self, x):
        X = x[: self.N * self.nv].reshape(self.N, self.nv)
        C = self.C
        return X[:, :C], X[:, C:2 * C], X[:, 2 * C], x[-3], x[-2], x[-1]

    # ---- properties ----------------------------------------------------------------------

    def stage_props(self, j: int, l, v, T: float):
        """K-values, liquid and vapor molar enthalpy of stage j, the components its reactions
        form (mol/s) and the molar enthalpy of each pumparound that returns its liquid at a
        fixed temperature ({pumparound: J/mol})."""
        pkg, P = self.pkg, self.P[j]
        lp = np.maximum(l, 0.0)
        x = lp / lp.sum() if lp.sum() > 0 else np.full(self.C, 1.0 / self.C)
        liq = pkg.phase_props("liquid", T, P, list(x))
        gen = self.generation(j, x, T, liq.V)
        ret = {k: pkg.flash(list(x), T=pa.T, P=self.P[pa.target]).H
               for k, pa in enumerate(self.pas) if pa.source == j and pa.T is not None}
        if j == 0 and self.total:
            # No vapor leaves a total condenser; use the vapor in equilibrium with its liquid.
            y = x
            for _ in range(2):
                vap = pkg.phase_props("vapor", T, P, list(y))
                K = np.exp(np.asarray(liq.lnphi) - np.asarray(vap.lnphi))
                y = K * x / np.sum(K * x)
            return K, liq.H, 0.0, gen, ret
        vp = np.maximum(v, 0.0)
        y = vp / vp.sum() if vp.sum() > 0 else x
        vap = pkg.phase_props("vapor", T, P, list(y))
        K = np.exp(np.asarray(liq.lnphi) - np.asarray(vap.lnphi))
        return K, liq.H, vap.H, gen, ret

    def generation(self, j: int, x, T: float, Vm: float) -> np.ndarray:
        """Components formed by the reactions on stage j, mol/s, at liquid mole fractions x
        and molar volume Vm."""
        if self.holdup[j] <= 0.0:
            return np.zeros(self.C)
        rates = np.array([r.rate(self.pkg, T, self.P[j], x, Vm, "liquid")
                          for r in self.reactions])
        return self.holdup[j] * (rates @ self.nu)

    def withdrawal_ratios(self, j: int, L: float, V: float) -> dict:
        """Each withdrawal of stage j over the flow continuing from it (L down, V up). A
        fraction is of the stage's whole outflow in that phase, every withdrawal included."""
        out = {}
        for items, flow in ((self.liq_w[j], L), (self.vap_w[j], V)):
            if not items:
                continue
            flow = max(flow, 1e-12 * self.F)
            rates = sum(r for _, r, f in items if r is not None)
            fsum = sum(f for _, r, f in items if f is not None)
            outflow = (flow + rates) / (1.0 - fsum)
            for key, r, f in items:
                out[key] = (r if r is not None else f * outflow) / flow
        return out

    def all_props(self, x, stages=None, base=None):
        l, v, T, *_ = self.unpack(x)
        props = list(base) if base is not None else [None] * self.N
        for j in range(self.N) if stages is None else stages:
            props[j] = self.stage_props(j, l[j], v[j], T[j])
        return props

    # ---- equations -----------------------------------------------------------------------

    def residual(self, x, props, raw_energy: bool = False):
        N, C = self.N, self.C
        l, v, _, s, Qc, Qr = self.unpack(x)
        L, V = l.sum(axis=1), v.sum(axis=1)
        K = np.array([p[0] for p in props])
        h = np.array([p[1] for p in props])
        Hv = np.array([p[2] for p in props])
        gen = np.array([p[3] for p in props])
        xm = l / np.where(L > 0, L, 1.0)[:, None]
        ym = v / np.where(V > 0, V, 1.0)[:, None]
        draw = np.zeros(N)
        if self.total:
            draw[0] = s
        draw_v = np.zeros(N)
        f_ret = np.zeros((N, C))  # pumparound returns
        H_ret = np.zeros(N)
        for j in self.withdrawing:
            for (kind, k), r in self.withdrawal_ratios(j, L[j], V[j]).items():
                if kind == "pa" or self.draws[k].phase == "liquid":
                    draw[j] += r
                else:
                    draw_v[j] += r
                if kind == "pa":
                    pa = self.pas[k]
                    f_ret[pa.target] += r * l[j]
                    H_ret[pa.target] += r * L[j] * (
                        props[j][4][k] if pa.T is not None else h[j]) + (pa.duty or 0.0)
        zero = np.zeros((1, C))
        l_in = np.vstack([zero, l[:-1]])
        v_in = np.vstack([v[1:], zero])

        R = np.empty((N, self.nv))
        R[:, :C] = ((1.0 + draw)[:, None] * l + (1.0 + draw_v)[:, None] * v - l_in - v_in
                    - self.f - f_ret - gen) / self.F

        y_in = np.vstack([ym[1:], self.y_bottom_in if self.y_bottom_in is not None else zero])
        eta = self.eta[:, None]
        R[:, C:2 * C] = ym - eta * K * xm - (1.0 - eta) * y_in
        if self.total:
            R[0, C:2 * C] = v[0] / self.F

        Q = np.zeros(N)
        if self.has_cond:
            Q[0] = Qc
        if self.has_reb:
            Q[-1] = Qr
        out = (1.0 + draw) * L * h + (1.0 + draw_v) * V * Hv
        energy = out - np.concatenate([[0.0], (L * h)[:-1]]) \
            - np.concatenate([(V * Hv)[1:], [0.0]]) - self.HF - Q - self.Qs - H_ret \
            + gen @ self.hf
        if raw_energy:
            return energy
        R[:, 2 * C] = energy / self.E

        g = np.zeros(3)
        g[0] = np.sum(K[0] * xm[0]) - 1.0 if self.total else s
        g[1] = 0.0 if self.has_cond else Qc / self.E
        g[2] = 0.0 if self.has_reb else Qr / self.E
        D, B = (s * L[0] if self.total else V[0]), L[-1]
        for row, name, value in self.spec_rows:
            if name == "reflux_ratio":
                g[row] = (L[0] - value * D) / self.F
            elif name == "distillate_rate":
                g[row] = (D - value) / self.F
            elif name == "bottoms_rate":
                g[row] = (B - value) / self.F
            else:  # boilup_ratio
                g[row] = (V[-1] - value * B) / self.F
        return np.concatenate([R.ravel(), g])

    # ---- Jacobian ------------------------------------------------------------------------

    def jacobian(self, x, props, r0):
        """Finite-difference Jacobian. Stages three apart share an evaluation; the bottom stage
        and the global unknowns get their own, because the global rows see stage 0 and N-1,
        and so does the source stage of each pumparound, which reaches its return stage."""
        N, nv = self.N, self.nv
        rows, cols, vals = [], [], []
        steps = self._steps(x)
        src = set(self.pa_sources)
        groups = [[j for j in range(m, N - 1, 3) if j not in src] for m in range(3)]
        groups += [[j] for j in sorted(src | {N - 1})]
        n_stage_rows = N * nv
        for group in groups:
            if not group:
                continue
            for k in range(nv):
                xp = x.copy()
                idx = [j * nv + k for j in group]
                xp[idx] += steps[idx]
                props_p = self.all_props(xp, stages=group, base=props)
                dr = self.residual(xp, props_p) - r0
                for j in group:
                    col = j * nv + k
                    lo, hi = max(j - 1, 0) * nv, min(j + 2, N) * nv
                    touched = list(range(lo, hi))
                    for b in sorted({pa.target for pa in self.pas if pa.source == j}):
                        if not lo <= b * nv < hi:
                            touched += list(range(b * nv, (b + 1) * nv))
                    if j in (0, N - 1):
                        touched += list(range(n_stage_rows, self.size))
                    d = dr[touched] / steps[col]
                    nz = d != 0.0
                    rows.extend(np.asarray(touched)[nz])
                    cols.extend([col] * int(nz.sum()))
                    vals.extend(d[nz])
        for col in range(n_stage_rows, self.size):  # s, Qc, Qr: properties are unaffected
            xp = x.copy()
            xp[col] += steps[col]
            d = (self.residual(xp, props) - r0) / steps[col]
            nz = np.nonzero(d)[0]
            rows.extend(nz)
            cols.extend([col] * len(nz))
            vals.extend(d[nz])
        return sp.csc_matrix((vals, (rows, cols)), shape=(self.size, self.size))

    def _steps(self, x):
        steps = np.empty_like(x)
        X = steps[: self.N * self.nv].reshape(self.N, self.nv)
        Xv = x[: self.N * self.nv].reshape(self.N, self.nv)
        C = self.C
        X[:, :2 * C] = 1e-7 * np.maximum(np.abs(Xv[:, :2 * C]), 1e-6 * self.F)
        X[:, 2 * C] = 1e-6 * Xv[:, 2 * C]
        steps[-3] = 1e-7 * max(abs(x[-3]), 1e-3)
        steps[-2:] = 1e-7 * np.maximum(np.abs(x[-2:]), self.E)
        return steps

    def describe_row(self, i: int) -> str:
        if i >= self.N * self.nv:
            g = i - self.N * self.nv
            if g == 0:
                return "condenser bubble point"
            spec = next((name for row, name, _ in self.spec_rows if row == g), None)
            return SPECS[spec] + " specification" if spec else "duty"
        j, k = divmod(i, self.nv)
        kind = ("material balance" if k < self.C else
                "phase equilibrium" if k < 2 * self.C else "energy balance")
        return f"{kind} on stage {j + 1}"


# ---- initial estimate ------------------------------------------------------------------------

def _withdrawals(m: _Model, L, V):
    """Side draw flows per stage (liquid, vapor) and pumparound flows at total flows L, V."""
    wl, wv, pa = np.zeros(m.N), np.zeros(m.N), np.zeros(len(m.pas))
    L, V = np.maximum(L, 1e-12 * m.F), np.maximum(V, 1e-12 * m.F)  # as withdrawal_ratios
    for j in m.withdrawing:
        for (kind, k), r in m.withdrawal_ratios(j, L[j], V[j]).items():
            if kind == "pa":
                pa[k] = r * L[j]
            elif m.draws[k].phase == "liquid":
                wl[j] += r * L[j]
            else:
                wv[j] += r * V[j]
    return wl, wv, pa


def _cmo_flows(m: _Model, liq, vap, wl, wv, pa):
    """Constant-molar-overflow flows from the specifications with the given withdrawals.
    Returns L, V, reflux Lr, boilup Vb, top and bottom product."""
    Fv = (vap[1:].sum() if m.has_cond else vap.sum()) - wv.sum()
    Fl = liq.sum() - wl.sum()
    # Unknowns: reflux Lr (with condenser) and boilup Vb (with reboiler). D = Vb + Fv - Lr,
    # B = Lr + Fl - Vb. Each specification is linear in them.
    names = ["Lr"] * m.has_cond + ["Vb"] * m.has_reb
    A, rhs = [], []
    for _, spec, val in m.spec_rows:
        coef = {"reflux_ratio": ({"Lr": 1.0 + val, "Vb": -val}, val * Fv),
                "distillate_rate": ({"Lr": -1.0, "Vb": 1.0}, val - Fv),
                "bottoms_rate": ({"Lr": 1.0, "Vb": -1.0}, val - Fl),
                "boilup_ratio": ({"Lr": -val, "Vb": 1.0 + val}, val * Fl)}[spec]
        A.append([coef[0][n] for n in names])
        rhs.append(coef[1])
    u = dict(zip(names, np.linalg.solve(np.array(A), np.array(rhs)) if names else []))
    Lr, Vb = u.get("Lr", 0.0), u.get("Vb", 0.0)

    L = Lr + np.cumsum(liq - wl)
    L[-1] -= Vb
    V = Vb + np.cumsum((vap - wv)[::-1])[::-1]
    for k, p in enumerate(m.pas):  # liquid carried past the stages between source and return
        a, b = p.source, p.target
        lo = max(min(a, b), 1 if m.has_cond else 0)  # the reflux stays at its specification
        if b < a:
            L[lo:a] += pa[k]
        else:
            L[lo:b] -= pa[k]
    D = V[1] - Lr if m.has_cond else V[0]
    if m.has_cond:
        V[0] = 0.0 if m.total else D
    return L, V, Lr, Vb, D, L[-1]


def _initial(m: _Model) -> np.ndarray:
    """Constant-molar-overflow flows from the specifications, then sweeps of tridiagonal
    component balances: for distillation with bubble-point temperatures and flows from enthalpy
    balances (Wang-Henke), for absorbers at a fixed linear temperature profile (sum-rates),
    because their gas may have no bubble point at all. Side draws given as fractions and
    pumparounds follow the flows; reactions enter the sweeps with the rates of the previous
    one."""
    N, C, pkg = m.N, m.C, m.pkg
    vap = m.fv.sum(axis=1)
    liq = m.f.sum(axis=1) - vap
    wl, wv, pa = np.zeros(N), np.zeros(N), np.zeros(len(m.pas))
    for _ in range(4):  # withdrawals given as fractions depend on the flows
        L, V, Lr, Vb, D, B = _cmo_flows(m, liq, vap, wl, wv, pa)
        wl, wv, pa = _withdrawals(m, L, V)
    L, V, Lr, Vb, D, B = _cmo_flows(m, liq, vap, wl, wv, pa)
    problems = []
    if D <= 1e-9 * m.F:
        problems.append(f"top product {D * 3.6:.4g} kmol/h")
    if B <= 1e-9 * m.F:
        problems.append(f"bottom product {B * 3.6:.4g} kmol/h")
    if m.has_cond and Lr <= 0.0:
        problems.append(f"reflux {Lr * 3.6:.4g} kmol/h")
    if m.has_reb and Vb <= 0.0:
        problems.append(f"boilup {Vb * 3.6:.4g} kmol/h")
    if problems:
        drawn = wl.sum() + wv.sum()
        raise SimulationError(
            "the specifications cannot be met with this feed: they imply "
            + ", ".join(problems),
            hint=f"the total feed is {m.F * 3.6:.4g} kmol/h"
            + (f" and the side draws take {drawn * 3.6:.4g} kmol/h" if drawn > 0 else "")
            + "; choose product rates between zero and the feed, and positive reflux and "
            "boilup ratios",
        )
    inner = range(1 if m.has_cond else 0, N - 1 if m.has_reb else N)
    dry = [j for j in inner if L[j] <= 0.0 or (V[j] <= 0.0 and j > 0)]
    if dry and (m.draws or m.pas):
        j = dry[0]
        raise SimulationError(
            f"the side draws and pumparounds leave no {'liquid' if L[j] <= 0 else 'vapor'} "
            f"flowing on stage {j + 1} (estimate {min(L[j], V[j]) * 3.6:.4g} kmol/h)",
            hint="lower the side draw or pumparound flows, or raise the reflux or boilup",
        )

    # Component balance sweeps with fixed flows.
    z = m.f.sum(axis=0) / m.F
    bubble = m.has_cond or m.has_reb
    if bubble:
        P_mean = float(np.mean(m.P))
        try:
            T0, K0 = bubble_point(pkg, z, P_mean)
        except SimulationError:
            # A feed carrying dissolved light gases (traces of hydrogen in a benzene column)
            # has no bubble point at column pressure; its dew point is a usable start.
            T0, K0 = dew_point(pkg, z, P_mean)
        T = np.full(N, T0)
        K = np.tile(K0, (N, 1))
    else:
        top = min(m.setup.feeds, key=lambda fd: fd.stage)
        bottom = max(m.setup.feeds, key=lambda fd: fd.stage)
        T = np.linspace(top.T, bottom.T, N)
        # K between the liquid and the vapor fed in; the overall composition may have no
        # liquid root at all, which would give the trivial K = 1.
        fl, fvap = (m.f - m.fv).sum(axis=0), m.fv.sum(axis=0)
        x_ref = fl / fl.sum() if fl.sum() > 0 else z
        y_ref = fvap / fvap.sum() if fvap.sum() > 0 else z
        K = np.array([pkg.k_values(T[j], m.P[j], list(x_ref), list(y_ref)) for j in range(N)])
    draw = D / Lr if m.total else 0.0
    xs = np.tile(z, (N, 1))
    gen = np.zeros((N, C))
    reactive = np.flatnonzero(m.holdup)
    for _ in range(20 if bubble else 40):
        # Reactions at the previous composition: production on the right-hand side,
        # consumption (proportional to x) on the diagonal, which keeps x positive.
        for j in reactive:
            Vm = pkg.phase_props("liquid", T[j], m.P[j], list(xs[j])).V
            gen[j] = 0.5 * gen[j] + 0.5 * m.generation(j, xs[j], T[j], Vm)
        consume = np.maximum(-gen, 0.0) / np.maximum(xs, 1e-14)
        rhs = m.f + np.maximum(gen, 0.0)
        for k, p in enumerate(m.pas):
            rhs[p.target] += pa[k] * xs[p.source]
        out_l = L + wl
        for k, p in enumerate(m.pas):
            out_l[p.source] += pa[k]
        ab = np.zeros((3, N))
        x_new = np.empty((N, C))
        for i in range(C):
            ab[0, 1:] = -(V[1:] * K[1:, i])  # super-diagonal: vapor from below
            ab[1, :] = out_l + (V + wv) * K[:, i] + consume[:, i]
            ab[1, 0] += draw * L[0]
            ab[2, :-1] = -L[:-1]  # sub-diagonal: liquid from above
            x_new[:, i] = solve_banded((1, 1), ab, rhs[:, i])
        x_new = np.maximum(x_new, 1e-14)
        if bubble:
            xs = x_new / x_new.sum(axis=1, keepdims=True)
            T_old, V_old = T.copy(), V.copy()
            for j in range(N):
                try:
                    T[j], K[j] = bubble_point(pkg, xs[j], m.P[j], T_guess=T[j])
                except SimulationError:
                    # No bubble point, e.g. a top stage estimated as almost pure hydrogen
                    # and methane: keep this stage's temperature, Newton corrects it.
                    K[j] = pkg.k_values(T[j], m.P[j], list(xs[j]), list(xs[j]))
            if m.has_cond:
                L, V = _enthalpy_flows(m, T, xs, K, L, V, Lr, D, wl, wv, pa, gen)
                wl, wv, pa = _withdrawals(m, L, V)
            change = np.max(np.abs(V - V_old)) / m.F
            if np.max(np.abs(T - T_old)) < 0.05 and change < 1e-4:
                break
        else:
            # Sum-rates: the component flows give new total flows; K follows at fixed T.
            l, v = x_new * L[:, None], K * x_new * V[:, None]
            L_new, V_new = l.sum(axis=1), v.sum(axis=1)
            change = max(np.max(np.abs(L_new / L - 1.0)), np.max(np.abs(V_new / V - 1.0)))
            L, V = L_new, V_new
            xs, ys = l / L[:, None], v / V[:, None]
            K = np.array([pkg.k_values(T[j], m.P[j], list(xs[j]), list(ys[j]))
                          for j in range(N)])
            if change < 1e-6:
                break

    ys = K * xs
    ys /= ys.sum(axis=1, keepdims=True)
    x0 = np.zeros(m.size)
    X = x0[: N * m.nv].reshape(N, m.nv)
    X[:, :C] = xs * L[:, None]
    X[:, C:2 * C] = ys * V[:, None]
    X[:, 2 * C] = T
    x0[-3] = draw
    # Duties that close the condenser and reboiler energy balances of the estimate.
    energy = m.residual(x0, m.all_props(x0), raw_energy=True)
    x0[-2] = energy[0] if m.has_cond else 0.0
    x0[-1] = energy[-1] if m.has_reb else 0.0
    return x0


def _enthalpy_flows(m: _Model, T, xs, K, L, V, Lr: float, D: float, wl, wv, pa, gen):
    """Vapor and liquid flows from enthalpy balances around the top of the column, with
    reflux and distillate held. Keeps the previous flows if the balances give nonsense.

    The envelope over stages 0..j: V_j+1 H_j+1 + HF + Qc + E_in = L_j h_j + D h_D + E_out
    with L_j = V_j+1 + F + M_in - D - M_out; side draws, pumparounds that cross the envelope,
    stage duties and reactions make up the M and E terms."""
    N, pkg = m.N, m.pkg
    ys = K * xs
    ys /= ys.sum(axis=1, keepdims=True)
    h = np.array([pkg.phase_props("liquid", T[j], m.P[j], list(xs[j])).H for j in range(N)])
    H = np.array([pkg.phase_props("vapor", T[j], m.P[j], list(ys[j])).H for j in range(N)])
    h_d = h[0] if m.total else H[0]
    # Per stage: material and enthalpy entering besides feeds and stage-to-stage flows.
    m_in = gen.sum(axis=1) - wl - wv
    e_in = m.Qs - gen @ m.hf - wl * h - wv * H
    Fs = np.cumsum(m.f.sum(axis=1) + m_in)
    HFs = np.cumsum(m.HF + e_in)
    pa_m, pa_e = np.zeros(N), np.zeros(N)  # pumparound terms per envelope 0..j
    for k, p in enumerate(m.pas):
        a, b = p.source, p.target
        q = p.duty or 0.0
        if p.T is not None:
            q = pa[k] * (pkg.flash(list(xs[a]), T=p.T, P=m.P[b]).H - h[a])
        for j in range(N):
            if b <= j < a:  # returns into the envelope
                pa_m[j] += pa[k]
                pa_e[j] += pa[k] * h[a] + q
            elif a <= j < b:  # leaves it
                pa_m[j] -= pa[k]
                pa_e[j] -= pa[k] * h[a]
            elif a <= j and b <= j:
                pa_e[j] += q
    Fs, HFs = Fs + pa_m, HFs + pa_e
    # The condenser balance (envelope j = 0, L_0 = Lr) gives its duty.
    V1 = Lr + D - Fs[0]
    q_cond = Lr * h[0] + D * h_d - V1 * H[1] - HFs[0]
    V_new, L_new = V.copy(), L.copy()
    V_new[1] = V1
    for j in range(1, N - 1):
        V_new[j + 1] = ((Fs[j] - D) * h[j] + D * h_d - HFs[j] - q_cond) / (H[j + 1] - h[j])
        L_new[j] = V_new[j + 1] - D + Fs[j]
    if np.all(V_new[1:] > 0) and np.all(L_new > 0):
        return L_new, V_new
    return L, V


# ---- Newton --------------------------------------------------------------------------------

def _step(m: _Model, x, dx, alpha: float):
    """Damped step that keeps flows (and the distillate draw ratio) positive: a flow that
    would turn negative drops to a tenth instead."""
    xn = x + alpha * dx
    positive = np.zeros(m.size, dtype=bool)
    positive[: m.N * m.nv] = True
    positive[2 * m.C:m.N * m.nv:m.nv] = False  # temperatures
    positive[-3] = m.total
    bad = positive & (xn <= 0.0)
    xn[bad] = 0.1 * x[bad]
    return xn


def _max_alpha(m: _Model, x, dx) -> float:
    """Largest fraction of the Newton step within the trust limits: no stage temperature moves
    more than _MAX_DT and no total stage flow more than doubles (or grows by more than the
    feed). The step is shortened as a whole, so it keeps its direction; clipping single
    variables would mix the Newton direction with arbitrary ones."""
    N, C, nv = m.N, m.C, m.nv
    X, D = x[: N * nv].reshape(N, nv), dx[: N * nv].reshape(N, nv)
    alpha = min(1.0, _MAX_DT / max(np.max(np.abs(D[:, 2 * C])), 1e-30))
    for part in (slice(0, C), slice(C, 2 * C)):  # liquid, vapor
        total, change = X[:, part].sum(axis=1), D[:, part].sum(axis=1)
        room = np.maximum(total, m.F)
        grow = change > room
        if np.any(grow):
            alpha = min(alpha, float(np.min(room[grow] / change[grow])))
    return alpha


def solve_column(pkg: PropertyPackage, setup: ColumnSetup,
                 start: ColumnSolution | None = None) -> ColumnSolution:
    """Solve the column. ``start``, a solution of the same column for a nearby feed (the
    previous pass of a recycle loop), replaces the start estimate: Newton then usually needs
    two or three iterations instead of building the estimate from bubble points. If Newton
    gets nowhere from there, the column is solved again from its own estimate."""
    m = _Model(pkg, setup)
    if start is not None and start.restart is not None:
        key, x_old, F_old = start.restart
        if key == _restart_key(m):
            x = x_old.copy()
            flows = np.ones(m.size, dtype=bool)
            flows[2 * m.C:m.N * m.nv:m.nv] = False  # temperatures
            flows[-3] = False  # the distillate draw ratio
            x[flows] *= m.F / F_old  # flows and duties follow the feed
            try:
                return _newton(m, x, max(setup.max_iterations // 3, 1))
            except SimulationError:
                pass
    return _newton(m, _initial(m), setup.max_iterations)


def _restart_key(m: _Model) -> tuple:
    return (m.N, m.C, m.setup.condenser, m.has_reb,
            tuple((d.stage, d.phase) for d in m.draws),
            tuple((p.source, p.target) for p in m.pas), tuple(np.flatnonzero(m.holdup)))


def _newton(m: _Model, x, max_iterations: int) -> ColumnSolution:
    props = m.all_props(x)
    r = m.residual(x, props)
    for iteration in range(1, max_iterations + 1):
        J = m.jacobian(x, props, r)
        with np.errstate(all="ignore"):
            dx = spsolve(J, -r)
        if not np.all(np.isfinite(dx)):
            raise SimulationError(
                "the column equations are singular at the current estimate",
                hint="check that every component can leave the column and that the specs "
                "are independent (e.g. not a distillate and a bottoms rate together)",
            )
        # Non-monotone line search: a full Newton step may raise the residual for a while
        # (the K-values are exponential in T) and still lead to fast convergence, so only
        # steps that make it much worse are cut back.
        alpha = _max_alpha(m, x, dx)
        norm = np.linalg.norm(r)
        for attempt in range(8):
            xn = _step(m, x, dx, alpha)
            try:
                props_n = m.all_props(xn)
                rn = m.residual(xn, props_n)
            except SimulationError:
                alpha *= 0.5
                continue
            if np.linalg.norm(rn) < _MAX_GROWTH * norm or attempt == 7:
                break
            alpha *= 0.5
        else:
            raise SimulationError(
                "the column solver left the range where properties can be evaluated",
                hint="check the column pressure and feed conditions, or start from specs "
                "closer to a converged case",
            )
        x, props, r = xn, props_n, rn
        if np.max(np.abs(r)) < m.setup.tolerance:
            return _solution(m, x, props, iteration)

    worst = int(np.argmax(np.abs(r)))
    if m.has_cond or m.has_reb:
        hint = ("the specifications may be infeasible for this column: try a larger reflux "
                "ratio, more stages, a distillate rate closer to the total feed of the "
                "components that should go overhead, or a feed stage nearer the middle")
    else:
        hint = ("check that the liquid and gas are fed at sensible temperatures and that the "
                "solvent can take up the solute; changing the solvent flow or the number of "
                "stages often helps")
    raise SimulationError(
        f"the column did not converge in {max_iterations} Newton iterations "
        f"(largest error: {m.describe_row(worst)}, {abs(r[worst]):.2g} scaled)", hint=hint,
    )


def _solution(m: _Model, x, props, iterations: int) -> ColumnSolution:
    l, v, T, s, Qc, Qr = m.unpack(x)
    l, v = np.maximum(l, 0.0), np.maximum(v, 0.0)
    K = np.array([p[0] for p in props])
    L, V = l.sum(axis=1), v.sum(axis=1)
    xm = l / L[:, None]
    ym = np.where(V[:, None] > 0, v / np.where(V > 0, V, 1.0)[:, None], K * xm)
    ym /= ym.sum(axis=1, keepdims=True)
    zero = np.zeros(m.C)
    h = np.array([p[1] for p in props])
    draws, pa_flows, pa_duties, pa_H = [None] * len(m.draws), [], [], []
    ratios = {key: r for j in m.withdrawing
              for key, r in m.withdrawal_ratios(j, L[j], V[j]).items()}
    for k, d in enumerate(m.draws):
        draws[k] = ratios["draw", k] * (l if d.phase == "liquid" else v)[d.stage]
    for k, pa in enumerate(m.pas):
        r, a = ratios["pa", k], pa.source
        H_in = r * L[a] * h[a]
        H_ret = r * L[a] * props[a][4][k] if pa.T is not None else H_in + (pa.duty or 0.0)
        pa_flows.append(r * l[a])
        pa_duties.append(float(H_ret - H_in))
        pa_H.append(float(H_ret))
    return ColumnSolution(
        T=T.copy(), P=m.P.copy(), L=L, V=V, x=xm, y=ym,
        h=h, Hv=np.array([p[2] for p in props]),
        distillate=s * l[0] if m.total else zero,
        vapor_distillate=zero if m.total else v[0].copy(),
        bottoms=l[-1].copy(),
        condenser_duty=float(Qc) if m.has_cond else 0.0,
        reboiler_duty=float(Qr) if m.has_reb else 0.0,
        iterations=iterations,
        side_draws=draws, pumparound_flows=pa_flows, pumparound_duties=pa_duties,
        pumparound_return_H=pa_H,
        generation=np.array([p[3] for p in props]) if m.reactions else None,
        restart=(_restart_key(m), x.copy(), m.F),
    )
