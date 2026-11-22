"""Dynamic distillation column: stage hold-ups that start from the steady MESH solution.

Each stage holds liquid and a little vapor. Its states are the component hold-ups N_ij and
the enthalpy H_j of its contents; a PH flash at the stage pressure gives temperature and
phase split, as for a vessel (units.Vessel). Pressures stay at their steady values (a
flow-driven model, as in the vessels).

- Trays: the liquid leaves over the weir (Francis: Q = 1.84 l_w h_ow^1.5, m3/s), the vapor
  as fast as it forms, with a short residence time. The tray holds the clear liquid on its
  active area; downcomers are not modelled.
- Condenser with reflux drum (stage 1): the condenser duty and the reflux and distillate
  flows are settings that controllers can move; vapor that the condenser does not condense
  leaves as vapor distillate.
- Reboiler sump (last stage): the reboiler duty and the bottoms flow are settings.
- Without a condenser or reboiler the top or bottom stage is a tray.

The start reproduces the steady column: tray liquid from the weir formula at the steady
liquid flow, drum and sump liquid from residence times, vapor from the steady vapor flow.
Side draws, pumparounds, stage duties and reactive stages are not modelled yet.
"""

from __future__ import annotations

import math

import numpy as np

from ..errors import SimulationError
from ..units.base import Param, ResultField
from .units import DYNAMIC_MODELS, DynamicModel, _pinned

ACTIVE = 0.76  # active area over column cross section (two downcomers of 12 %)
WEIR = 0.73  # weir length over diameter for that downcomer area
F_FACTOR = 2.0  # Pa^0.5, vapor load u sqrt(rho_v) used to size a column without a diameter


class DynamicColumn(DynamicModel):
    params = (
        Param("diameter", "Column diameter", "length", min=0.0,
              help="Inner diameter of the trays. Defaults to a size for a moderate vapor "
              "load (F-factor 2 Pa^0.5 on the active area) at the busiest stage."),
        Param("weir_height", "Weir height", "length", default=0.05, min=0.0,
              help="Height of the outlet weir; with the crest over it, it sets the liquid "
              "held on each tray."),
        Param("vapor_time", "Vapor residence time", "time", default=1.0, min=1e-3,
              help="How long vapor stays on a stage before it rises; short."),
        Param("drum_time", "Reflux drum residence time", "time", default=300.0, min=1.0,
              help="Liquid in the reflux drum at the start, as minutes of reflux plus "
              "distillate. The drum holds twice that when full."),
        Param("sump_time", "Sump residence time", "time", default=300.0, min=1.0,
              help="Liquid in the column sump at the start, as minutes of the liquid "
              "reaching it. The sump holds twice that when full."),
        Param("reflux_flow", "Reflux flow", "molar_flow", min=0.0,
              help="Liquid returned from the drum to the top tray. Defaults to the steady "
              "state; a controller often moves it."),
        Param("distillate_flow", "Distillate flow", "molar_flow", min=0.0,
              help="Liquid drawn from the drum. Defaults to the steady state; usually the "
              "drum level controller moves it."),
        Param("bottoms_flow", "Bottoms flow", "molar_flow", min=0.0,
              help="Liquid drawn from the sump. Defaults to the steady state; usually the "
              "sump level controller moves it."),
        Param("condenser_control", "Condenser", "choice", options=("total", "duty"),
              help="total: the cooling condenses all vapor reaching the condenser, as when "
              "the column pressure controller moves the coolant (the default for a total "
              "condenser); duty: the condenser duty is a setting, and vapor it leaves "
              "uncondensed vents as vapor distillate (the default for a partial condenser)."),
        Param("condenser_duty", "Condenser duty", "power",
              help="Heat added at the condenser (negative: removed) when the condenser runs "
              "on its duty. Defaults to the steady state."),
        Param("reboiler_duty", "Reboiler duty", "power",
              help="Heat added at the reboiler. Defaults to the steady state."),
        Param("control_stage", "Control stage", min=1.0,
              help="Stage whose temperature is reported as the control temperature. "
              "Defaults to the stage with the steepest steady temperature change."),
    )
    result_fields = (
        ResultField("control_T", "Control stage temperature", "temperature"),
        ResultField("top_T", "Top stage temperature", "temperature"),
        ResultField("bottom_T", "Bottom stage temperature", "temperature"),
        ResultField("drum_level", "Reflux drum level (fraction)", "fraction"),
        ResultField("sump_level", "Sump level (fraction)", "fraction"),
        ResultField("drum_holdup", "Reflux drum hold-up", "amount"),
        ResultField("sump_holdup", "Sump hold-up", "amount"),
        ResultField("reflux_flow", "Reflux flow", "molar_flow"),
        ResultField("distillate_flow", "Distillate flow", "molar_flow"),
        ResultField("bottoms_flow", "Bottoms flow", "molar_flow"),
        ResultField("boilup", "Boil-up", "molar_flow"),
        ResultField("vapor_distillate_flow", "Vapor distillate (uncondensed)", "molar_flow"),
        ResultField("condenser_duty", "Condenser duty", "power"),
        ResultField("reboiler_duty", "Reboiler duty", "power"),
    )

    # ---- start ---------------------------------------------------------------------------

    def initialize(self, inlets, outlets):
        u, pkg = self.unit, self.pkg
        sol = getattr(u, "_last_solution", None)
        if sol is None:
            raise SimulationError("the column has no steady solution to start from",
                                  where=u.id)
        self._check_supported()
        N, C = len(sol.T), pkg.n
        self.N, self.C, self.n = N, C, C + 1
        self.P = np.asarray(sol.P, dtype=float)
        self.cond = u.p["condenser"] != "none"
        self.reb = u.p["reboiler"] != "none"
        self.trays = [j for j in range(N) if not ((j == 0 and self.cond)
                                                   or (j == N - 1 and self.reb))]
        if not self.trays:
            raise SimulationError("the column has no trays between condenser and reboiler",
                                  where=u.id)
        L, V = np.asarray(sol.L, dtype=float), np.asarray(sol.V, dtype=float)
        x, y = np.asarray(sol.x, dtype=float), np.asarray(sol.y, dtype=float)
        D_liq = float(np.sum(sol.distillate))
        B = float(np.sum(sol.bottoms))
        vL = np.array([pkg.phase_props("liquid", sol.T[j], self.P[j], list(x[j])).V
                       for j in range(N)])

        # Tray geometry: given, or sized for the largest vapor load.
        D_col = self.s["diameter"]
        if D_col is None:
            area = 0.0
            for j in self.trays:
                vv = pkg.phase_props("vapor", sol.T[j], self.P[j], list(y[j]))
                rho = sum(yi * mw for yi, mw in zip(y[j], pkg.mws)) / 1000.0 / vv.V
                area = max(area, V[j] * vv.V * math.sqrt(rho) / F_FACTOR)
            D_col = math.sqrt(4.0 * area / ACTIVE / math.pi) if area > 0 else 1.0
            self.s["diameter"] = D_col
        self.active = ACTIVE * math.pi * D_col ** 2 / 4.0
        self.weir = WEIR * D_col

        # Settings that default to the steady state.
        defaults = {"reflux_flow": float(L[0]) if self.cond else 0.0,
                    "distillate_flow": D_liq if self.cond else 0.0,
                    "bottoms_flow": B if self.reb else 0.0,
                    "condenser_duty": float(sol.condenser_duty),
                    "reboiler_duty": float(sol.reboiler_duty)}
        for name, value in defaults.items():
            if self.s[name] is None:
                self.s[name] = value
        if self.s["condenser_control"] is None:
            self.s["condenser_control"] = "total" if u.p["condenser"] == "total" else "duty"
        if u.p["condenser"] == "partial" and self.s["condenser_control"] == "total":
            raise SimulationError(
                "a partial condenser cannot condense everything", where=u.id,
                hint="run it on its duty (condenser_control: duty), or use a total condenser")
        self.liquid_drum = self.cond and N > 1 and self.s["condenser_control"] == "total"
        if self.s["control_stage"] is None:
            dT = np.abs(np.diff(np.asarray(sol.T)))
            self.s["control_stage"] = float(int(np.argmax(dT)) + 2 if len(dT) else 1)
        if not 1 <= int(self.s["control_stage"]) <= N:
            raise SimulationError(f"the control stage must be between 1 and {N}", where=u.id)

        # Liquid on each stage: trays from the weir crest, drum and sump from residence times.
        nL = np.zeros(N)
        for j in self.trays:
            q = L[j] * vL[j]  # m3/s over the weir
            crest = (q / (1.84 * self.weir)) ** (2.0 / 3.0) if q > 0 else 0.0
            nL[j] = self.active * (self.s["weir_height"] + crest) / vL[j]
        self.capacity = {}
        if self.cond:
            nL[0] = max(L[0] + D_liq, 1e-9 * max(V.max(), 1e-9)) * self.s["drum_time"]
            self.capacity[0] = 2.0 * nL[0] * vL[0]
        if self.reb:
            inflow = L[N - 2] if N > 1 else B
            nL[-1] = max(inflow, 1e-9) * self.s["sump_time"]
            self.capacity[N - 1] = 2.0 * nL[-1] * vL[-1]
        nV = V * self.s["vapor_time"]
        H = nL * np.asarray(sol.h) + nV * np.asarray(sol.Hv)
        Nij = nL[:, None] * x + nV[:, None] * y

        self.T_last = np.asarray(sol.T, dtype=float).copy()
        self.slope = np.full(N, np.nan)  # dH/dT of each stage's contents, J/mol/K
        self._memo: list[tuple | None] = [None] * N
        self.T0 = np.asarray(sol.T, dtype=float)
        # One temperature series per stage, besides the class-wide results.
        self.result_fields = self.result_fields + tuple(
            ResultField(f"T_stage_{j + 1}", f"Temperature stage {j + 1}", "temperature")
            for j in range(N))
        return np.column_stack([Nij, H]).ravel()

    def _check_supported(self) -> None:
        p = self.unit.p
        missing = [label for name, label in (
            ("side_draws", "side draws"), ("pumparounds", "pumparounds"),
            ("stage_duties", "stage duties"), ("reactive_stages", "reactive stages"))
            if p.get(name)]
        if missing:
            raise SimulationError(
                f"the dynamic column does not model {', '.join(missing)} yet",
                where=self.unit.id, hint="simulate the column without them, or keep it "
                "steady (not under dynamics.units)")
        if (p.get("efficiency") or 1.0) < 1.0:
            raise SimulationError(
                "the dynamic column needs ideal stages (Murphree efficiency 1)",
                where=self.unit.id, hint="add stages instead of lowering the efficiency")

    def scale(self, x0):
        X = x0.reshape(self.N, self.n)
        total = X[:, :-1].sum(axis=1)
        floor = 1e-6 * max(total.max(), 1e-9)
        Ns = np.maximum(X[:, :-1], floor)
        Hs = np.maximum(np.abs(X[:, -1]), total * 8.314 * 300.0)
        return np.column_stack([Ns, Hs]).ravel()

    def pattern(self) -> np.ndarray:
        """Each stage's derivatives depend on its own state and its neighbours'."""
        n, N = self.n, self.N
        S = np.zeros((N * n, N * n), dtype=int)
        for j in range(N):
            lo, hi = max(j - 1, 0), min(j + 1, N - 1)
            S[j * n:(j + 1) * n, lo * n:(hi + 1) * n] = 1
        return S

    def result_states(self, name: str) -> np.ndarray | None:
        """Indices of the states a result depends on (None: any). A controller measuring
        one stage then couples to that stage only, which keeps the Jacobian cheap."""
        N, n = self.N, self.n
        stage = {"control_T": int(self.s["control_stage"]) - 1, "top_T": 0, "bottom_T": N - 1,
                 "drum_level": 0, "drum_holdup": 0, "reflux_flow": 0, "distillate_flow": 0,
                 "sump_level": N - 1, "sump_holdup": N - 1, "bottoms_flow": N - 1,
                 "boilup": N - 1}.get(name)
        if stage is None and name.startswith("T_stage_") and name[8:].isdigit():
            stage = int(name[8:]) - 1
        if stage is None or not 0 <= stage < N:
            return None
        return np.arange(stage * n, (stage + 1) * n)

    # ---- one instant ---------------------------------------------------------------------

    def _stage(self, j: int, s: np.ndarray):
        """Phase split of stage j's contents, cached on the exact state: the finite-difference
        Jacobian perturbs a few stages at a time and must see the others unchanged."""
        key = s.tobytes()
        memo = self._memo[j]
        if memo is not None and memo[0] == key:
            return memo[1]
        pkg, P = self.pkg, self.P[j]
        N = np.maximum(s[:-1], 0.0)
        total = N.sum()
        if total <= 0.0:
            raise SimulationError(f"stage {j + 1} ran dry", where=self.unit.id,
                                  hint="check the reflux, boil-up and level controllers")
        h = s[-1] / total
        if j == 0 and self.liquid_drum:
            # A condenser that condenses everything holds liquid only. Treating its drum as
            # liquid keeps it smooth: at its bubble point a flash would switch between
            # subcooled (no vapor) and two-phase on the slightest change.
            z = N / total
            T, vL = self._liquid_T(z, P, h)
            out = (T, z, z, total, 0.0, total * vL, h, h)
            self._memo[j] = (key, out)
            return out
        st = self._flash(j, list(N / total), P, h)
        VF = st.VF
        xs = np.asarray(st.x) if st.x else N / total
        ys = np.asarray(st.y) if st.y else N / total
        if 0.0 < VF < 1.0 and st.x and st.y:
            # Two phases: the split by the lever rule on the exact enthalpy. Near a pure
            # component H(T) is almost vertical, so T cannot carry the enthalpy to the last
            # digit; the phase enthalpies at that T can.
            hL = pkg.phase_props("liquid", st.T, P, list(xs)).H
            hV = pkg.phase_props("vapor", st.T, P, list(ys)).H
            VF = min(max((h - hL) / (hV - hL), 0.0), 1.0)
        else:
            hV = pkg.phase_props("vapor", st.T, P, list(ys)).H if VF > 0 else 0.0
            hL = (st.H - VF * hV) / (1.0 - VF) if VF < 1.0 else 0.0
        nL, nV = total * (1.0 - VF), total * VF
        vL = pkg.phase_props("liquid", st.T, P, list(xs)).V if nL > 0 else 0.0
        out = (st.T, xs, ys, nL, nV, nL * vL, hL, hV)
        self._memo[j] = (key, out)
        return out

    def _liquid_T(self, z: np.ndarray, P: float, h: float) -> tuple[float, float]:
        """Temperature and molar volume of a liquid of molar enthalpy h (Newton steps)."""
        pkg, zs = self.pkg, list(z)
        T = float(self.T_last[0])
        for _ in range(50):
            props = pkg.phase_props("liquid", T, P, zs)
            dh = pkg.phase_props("liquid", T + 0.01, P, zs).H - props.H
            step = (h - props.H) * 0.01 / dh if dh > 0 else 0.0
            step = min(max(step, -20.0), 20.0)
            T += step
            if abs(step) <= 1e-12 * T:
                self.T_last[0] = T
                return T, pkg.phase_props("liquid", T, P, zs).V
        raise SimulationError("could not find the reflux drum temperature",
                              where=self.unit.id)

    def _flash(self, j: int, zs: list[float], P: float, h: float):
        """The stage contents at molar enthalpy h, by secant steps on T with T-P flashes.

        Converged far tighter than the general PH flash: a stage holds little vapor, so its
        vapor fraction, and with it the vapor flow, magnifies any enthalpy error, and the
        integrator's finite-difference Jacobian divides by perturbations of 1e-7 and less.
        """
        pkg = self.pkg
        tol = 1e-10 * max(abs(h), 1e4)

        def f(T: float):
            st = pkg.flash(zs, T=T, P=P)
            return st.H - h, st

        T0 = float(self.T_last[j])
        f0, st = f(T0)
        if abs(f0) <= tol:
            return st
        slope = self.slope[j]
        if not np.isfinite(slope) or slope <= 0.0:
            fd, _ = f(T0 + 0.01)
            slope = max((fd - f0) / 0.01, 1.0)
        lo = hi = None
        width = math.inf
        T1 = T0 - f0 / slope
        for _ in range(100):
            T1 = min(max(T1, T0 - 20.0), T0 + 20.0, 5000.0)
            if lo is not None and hi is not None:
                # Secant steps stall on one end where H(T) is almost a step (a nearly pure
                # liquid boiling): bisect unless the last step at least halved the bracket.
                if not lo < T1 < hi or hi - lo > 0.5 * width:
                    T1 = 0.5 * (lo + hi)
                width = hi - lo
            f1, st = f(T1)
            if f1 < 0.0:
                lo = T1 if lo is None else max(lo, T1)
            else:
                hi = T1 if hi is None else min(hi, T1)
            narrow = lo is not None and hi is not None and hi - lo <= 1e-10 * T1
            if abs(f1) <= tol or narrow:  # narrow: H(T) almost vertical (see _stage)
                if T1 != T0 and f1 != f0:
                    self.slope[j] = (f1 - f0) / (T1 - T0)
                self.T_last[j] = T1
                return st
            s = (f1 - f0) / (T1 - T0) if T1 != T0 else slope
            if not np.isfinite(s) or s <= 0.0:
                s = slope
            T0, f0, T1 = T1, f1, T1 - f1 / s
        raise SimulationError(f"could not find the temperature on stage {j + 1}",
                              where=self.unit.id)

    def outlets(self, x, inlets=None):
        N, X = self.N, x.reshape(self.N, self.n)
        st = [self._stage(j, X[j]) for j in range(N)]
        T = np.array([s[0] for s in st])
        Lout = np.zeros(N)  # liquid leaving each stage downward (reflux from the drum)
        Vout = np.array([s[4] for s in st]) / self.s["vapor_time"]
        for j in self.trays:
            crest = max(st[j][5] / self.active - self.s["weir_height"], 0.0)
            q = 1.84 * self.weir * crest ** 1.5
            Lout[j] = q * st[j][3] / st[j][5] if st[j][5] > 0 else 0.0
        D = 0.0
        levels = {}
        for j, cap in self.capacity.items():
            levels[j] = st[j][5] / cap
        if self.cond:
            fade = min(1.0, levels[0] / 0.01)  # an empty drum delivers nothing
            Lout[0] = self.s["reflux_flow"] * fade
            D = self.s["distillate_flow"] * fade
        B = Lout[-1]
        if self.reb:
            B = self.s["bottoms_flow"] * min(1.0, levels[N - 1] / 0.01)
            Lout[-1] = 0.0
        Qc = self.s["condenser_duty"] if self.cond else 0.0
        if self.liquid_drum:
            # The cooling condenses the arriving vapor to drum liquid; nothing vents.
            Qc = Vout[1] * (st[0][6] - st[1][7])
            Vout[0] = 0.0
        self._Qc = Qc
        self._flows = (st, Lout, Vout, D, B)
        r = self.results
        r.clear()
        r.update(control_T=float(T[int(self.s["control_stage"]) - 1]), top_T=float(T[0]),
                 bottom_T=float(T[-1]), reflux_flow=float(Lout[0]) if self.cond else 0.0,
                 distillate_flow=D, bottoms_flow=B, boilup=float(Vout[-1]),
                 vapor_distillate_flow=float(Vout[0]),
                 condenser_duty=Qc,
                 reboiler_duty=self.s["reboiler_duty"] if self.reb else 0.0)
        if self.cond:
            r.update(drum_level=levels[0], drum_holdup=st[0][3])
        if self.reb:
            r.update(sump_level=levels[N - 1], sump_holdup=st[-1][3])
        for j in range(N):
            r[f"T_stage_{j + 1}"] = float(T[j])
        pkg, P = self.pkg, self.P
        return {
            "distillate": _pinned(pkg, D * st[0][1], T[0], P[0], 0.0),
            "vapor_distillate": _pinned(pkg, Vout[0] * st[0][2], T[0], P[0], 1.0),
            "bottoms": _pinned(pkg, B * st[-1][1], T[-1], P[-1], 0.0),
        }

    def derivatives(self, x, inlets):
        st, Lout, Vout, D, B = self._flows
        N, C = self.N, self.C
        dN = np.zeros((N, C))
        dH = np.zeros(N)
        u = self.unit
        ids = u.inlet_ids.get("in", [])
        for i, s in enumerate(inlets.get("in", [])):
            j = u._feed_stage(ids[i] if i < len(ids) else None)
            dN[j] += np.asarray(s.flows)
            dH[j] += s.enthalpy_flow
        for j in range(N):
            xj, yj, hL, hV = st[j][1], st[j][2], st[j][6], st[j][7]
            # Without a reboiler the bottom tray's weir flow is the bottoms (Lout[-1] = B).
            liq_out = Lout[j] + (D if j == 0 and self.cond else 0.0) + \
                (B if j == N - 1 and self.reb else 0.0)
            dN[j] -= liq_out * xj + Vout[j] * yj
            dH[j] -= liq_out * hL + Vout[j] * hV
            if j + 1 < N:  # liquid down to the next stage (the reflux from the drum)
                dN[j + 1] += Lout[j] * xj
                dH[j + 1] += Lout[j] * hL
            if j > 0:  # vapor up to the stage above
                dN[j - 1] += Vout[j] * yj
                dH[j - 1] += Vout[j] * hV
        if self.cond:
            dH[0] += self._Qc
        if self.reb:
            dH[-1] += self.s["reboiler_duty"]
        return np.column_stack([dN, dH]).ravel()


DYNAMIC_MODELS["DistillationColumn"] = DynamicColumn
