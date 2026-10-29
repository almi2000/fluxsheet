"""Distillation and absorption columns.

The models live in ``fluxsheet.columns``: a shortcut design method and one rigorous
equilibrium-stage solver that the distillation column and the absorber/stripper configure
differently.
"""

from __future__ import annotations

import math

import numpy as np

from ..columns import ColumnSetup, StageFeed, shortcut_design, solve_column
from ..columns.saturation import dew_point
from ..columns.stages import SPECS
from ..errors import SimulationError
from ..stream import Stream
from .base import IN, Param, Port, ResultField, UnitOperation, empty_stream, single_inlet
from .separators import _with_phase

CONDENSER_DUTY = ResultField("condenser_duty", "Condenser duty", "power")
REBOILER_DUTY = ResultField("reboiler_duty", "Reboiler duty", "power")
P_TOP = Param("P", "Top pressure", "pressure", min=1.0,
              help="Pressure of the top stage (the condenser, if any). Defaults to the feed "
              "pressure.")
DP_STAGE = Param("dP_stage", "Pressure drop per stage", "pressure_difference", min=0.0,
                 help="Pressure increase from each stage to the one below. Defaults to zero.")
EFFICIENCY = Param("efficiency", "Murphree efficiency", "fraction", min=0.01, max=1.0,
                   help="Vapor-phase Murphree efficiency of the stages between condenser and "
                   "reboiler. Defaults to 1 (equilibrium stages).")


def _product(pkg, flows, P: float, H: float, vf: float, T: float) -> Stream:
    """A product stream with the stage's enthalpy, so the column's energy balance closes."""
    flows = [float(f) for f in flows]
    if sum(flows) <= 0.0:
        return empty_stream(pkg, T, P)
    s = Stream.from_flash(pkg, flows, P=P, H=H)
    # At saturation a PH flash may land a hair on either side; pin the product's phase.
    return Stream(s.flows, _with_phase(s.state, vf, s.flows))


def _feed(stage: int, s: Stream) -> StageFeed:
    vap = s.total_flow * s.state.VF * np.asarray(s.state.y) if s.state.y else np.zeros(len(s.flows))
    return StageFeed(stage, np.asarray(s.flows, dtype=float), s.enthalpy_flow, vap, s.T)


class _StageColumn(UnitOperation):
    """Shared by the units built on the rigorous stage model."""

    def _pressures(self, P_top: float) -> np.ndarray:
        n = self.p["stages"]
        return P_top + np.arange(n) * (self.p.get("dP_stage") or 0.0)

    def _solve(self, pkg, setup: ColumnSetup):
        sol = solve_column(pkg, setup)
        stages = list(range(1, setup.stages + 1))
        self.add_profile("temperature", "Temperature", "Stage", stages,
                         {"T": ("temperature", sol.T.tolist())})
        self.add_profile("flows", "Liquid and vapor flow", "Stage", stages,
                         {"Liquid L": ("molar_flow", sol.L.tolist()),
                          "Vapor V": ("molar_flow", sol.V.tolist())})
        self.add_profile("liquid_composition", "Liquid composition", "Stage", stages,
                         {name: (None, sol.x[:, i].tolist()) for i, name in enumerate(pkg.names)})
        self.add_profile("vapor_composition", "Vapor composition", "Stage", stages,
                         {name: (None, sol.y[:, i].tolist()) for i, name in enumerate(pkg.names)})
        self.results["iterations"] = sol.iterations
        return sol


class ShortcutColumn(UnitOperation):
    type = "ShortcutColumn"
    label = "Shortcut column"
    category = "Columns"
    description = "Distillation design by Fenske, Underwood, Gilliland and Kirkbride: minimum "\
        "stages and reflux, the stages needed at your reflux, and the feed stage. A quick "\
        "start for a rigorous distillation column."
    inlets = (IN,)
    outlets = (Port("distillate", "Distillate"), Port("bottoms", "Bottoms"))
    params = (
        Param("light_key", "Light key", "component",
              help="The more volatile key component, mostly recovered in the distillate."),
        Param("heavy_key", "Heavy key", "component",
              help="The less volatile key component, mostly recovered in the bottoms."),
        Param("lk_recovery", "Light key to distillate", "fraction", default=0.99,
              min=0.0, max=1.0, help="Fraction of the light key feed leaving in the distillate."),
        Param("hk_recovery", "Heavy key to distillate", "fraction", default=0.01,
              min=0.0, max=1.0, help="Fraction of the heavy key feed leaving in the distillate."),
        Param("reflux_factor", "Reflux / minimum reflux", min=1.0,
              help="Reflux ratio as a multiple of the minimum; 1.1 to 1.5 is typical."),
        Param("reflux_ratio", "Reflux ratio", min=0.0, help="Absolute reflux ratio L/D."),
        Param("condenser", "Condenser", "choice", default="total", options=("total", "partial"),
              help="Total: liquid distillate. Partial: vapor distillate; the condenser is "
              "then an equilibrium stage."),
        Param("P", "Column pressure", "pressure", min=1.0,
              help="Defaults to the feed pressure."),
    )
    exclusive_groups = (("light_key",), ("heavy_key",), ("reflux_factor", "reflux_ratio"))
    required_groups = exclusive_groups
    result_fields = (
        ResultField("min_stages", "Minimum stages (Fenske)"),
        ResultField("min_reflux", "Minimum reflux ratio (Underwood)"),
        ResultField("reflux_ratio", "Reflux ratio"),
        ResultField("theoretical_stages", "Equilibrium stages (Gilliland)"),
        ResultField("stages", "Stages incl. condenser and reboiler"),
        ResultField("feed_stage", "Feed stage (Kirkbride)"),
        ResultField("alpha", "Relative volatility of the keys"),
        ResultField("q", "Feed quality q"),
        CONDENSER_DUTY,
        REBOILER_DUTY,
    )

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        P = self.p["P"] or s.P
        if s.total_flow <= 0.0:
            return {"distillate": empty_stream(pkg, s.T, P), "bottoms": empty_stream(pkg, s.T, P)}
        lk = self.component_index(pkg, self.p["light_key"], "light key")
        hk = self.component_index(pkg, self.p["heavy_key"], "heavy key")
        res = shortcut_design(
            pkg, s.flows, s.state.H, P, lk, hk, self.p["lk_recovery"], self.p["hk_recovery"],
            reflux_ratio=self.p["reflux_ratio"], reflux_factor=self.p["reflux_factor"],
        )
        D = res.distillate.sum()
        R = res.reflux_ratio
        x_d = list(res.distillate / D)
        total = self.p["condenser"] == "total"
        if total:
            # Saturated vapor of distillate composition is condensed to its bubble point.
            liquid = pkg.flash(x_d, P=P, VF=0.0)
            vapor = pkg.flash(x_d, P=P, VF=1.0)
            q_cond = (R + 1.0) * D * (liquid.H - vapor.H)
            distillate = Stream(list(res.distillate), liquid)
        else:
            # Vapor distillate at its dew point; the reflux is the liquid in equilibrium with
            # it, and the vapor from the top stage carries both.
            T_dew, K = dew_point(pkg, x_d, P)
            x_r = np.asarray(x_d) / K
            x_r /= x_r.sum()
            h_reflux = pkg.phase_props("liquid", T_dew, P, list(x_r)).H
            dew = pkg.flash(x_d, P=P, VF=1.0)
            y_top = (np.asarray(x_d) + R * x_r) / (R + 1.0)
            h_top = pkg.flash(list(y_top), P=P, VF=1.0).H
            q_cond = D * dew.H + R * D * h_reflux - (R + 1.0) * D * h_top
            distillate = Stream(list(res.distillate), dew)
        bottoms = Stream(list(res.bottoms), pkg.flash(list(res.bottoms / res.bottoms.sum()),
                                                      P=P, VF=0.0))
        q_reb = distillate.enthalpy_flow + bottoms.enthalpy_flow - s.enthalpy_flow - q_cond

        n_total = math.ceil(res.stages - 1e-9) + (1 if total else 0)
        feed_stage = round(res.rectifying_stages) + 1 + (1 if total else 0)
        self.results.update(
            min_stages=res.min_stages, min_reflux=res.min_reflux, reflux_ratio=R,
            theoretical_stages=res.stages, stages=n_total,
            feed_stage=min(max(feed_stage, 2), n_total - 1), alpha=float(res.alpha[lk]),
            q=res.q, condenser_duty=float(q_cond), reboiler_duty=float(q_reb),
        )
        return {"distillate": distillate, "bottoms": bottoms}


class DistillationColumn(_StageColumn):
    type = "DistillationColumn"
    label = "Distillation column"
    category = "Columns"
    description = "Rigorous distillation on equilibrium stages: material, equilibrium and "\
        "energy balances of every stage solved together. Stage 1 is the condenser, the last "\
        "stage the reboiler. Give one specification per condenser and reboiler."
    inlets = (IN,)
    outlets = (Port("distillate", "Distillate"), Port("bottoms", "Bottoms"),
               Port("vapor_distillate", "Vapor distillate"))
    params = (
        Param("stages", "Number of stages", "integer", default=10, min=2,
              help="Equilibrium stages, counting the condenser and the reboiler."),
        Param("feed_stage", "Feed stage", "integer", default=5, min=1,
              help="Stage the feed enters, counted from the top (condenser = 1)."),
        P_TOP,
        DP_STAGE,
        Param("condenser", "Condenser", "choice", default="total",
              options=("total", "partial", "none"),
              help="Total: liquid distillate. Partial: vapor distillate. None: the top "
              "vapor leaves as it is."),
        Param("reboiler", "Reboiler", "choice", default="kettle", options=("kettle", "none")),
        Param("reflux_ratio", "Reflux ratio", min=0.0,
              help="Reflux over distillate flow, L/D. Needs a condenser."),
        Param("boilup_ratio", "Boilup ratio", min=0.0,
              help="Vapor from the reboiler over bottoms flow, V/B. Needs a reboiler."),
        Param("distillate_rate", "Distillate rate", "molar_flow", min=0.0,
              help="Total top product flow (liquid or vapor distillate)."),
        Param("bottoms_rate", "Bottoms rate", "molar_flow", min=0.0),
        EFFICIENCY,
    )
    exclusive_groups = (("reflux_ratio", "boilup_ratio"), ("distillate_rate", "bottoms_rate"))
    result_fields = (CONDENSER_DUTY, REBOILER_DUTY, ResultField("reflux_ratio", "Reflux ratio"),
                     ResultField("boilup_ratio", "Boilup ratio"),
                     ResultField("iterations", "Newton iterations"))

    def __post_init__(self) -> None:
        super().__post_init__()
        has_cond = self.p["condenser"] != "none"
        has_reb = self.p["reboiler"] != "none"
        if self.p["reflux_ratio"] is not None and not has_cond:
            raise SimulationError("a reflux ratio needs a condenser", where=self.id,
                                  hint="choose a total or partial condenser, or specify a "
                                  "distillate or bottoms rate instead")
        if self.p["boilup_ratio"] is not None and not has_reb:
            raise SimulationError("a boilup ratio needs a reboiler", where=self.id,
                                  hint="choose a kettle reboiler, or specify a distillate or "
                                  "bottoms rate instead")
        given = [SPECS[name] for name in SPECS if self.p[name] is not None]
        needed = has_cond + has_reb
        if len(given) != needed:
            duties = {2: "one for the condenser and one for the reboiler",
                      1: "one, for the " + ("condenser" if has_cond else "reboiler"),
                      0: "none without condenser and reboiler"}[needed]
            raise SimulationError(
                f"{'under' if len(given) < needed else 'over'}-specified: this column takes "
                f"{duties}, but has {len(given)}{': ' + ', '.join(given) if given else ''}",
                where=self.id,
                hint={2: "set one of reflux ratio or boilup ratio and one of distillate rate "
                         "or bottoms rate; a reflux ratio and a distillate rate is typical",
                      1: "set exactly one ratio or rate that fits this configuration",
                      0: "remove the specifications, or add a condenser or reboiler"}[needed],
            )
        n, nf = self.p["stages"], self.p["feed_stage"]
        if n < 1 + has_cond + has_reb:
            raise SimulationError(f"{n} stages is too few", where=self.id,
                                  hint="count the condenser and reboiler as stages; use at "
                                  f"least {2 + has_cond + has_reb}")
        if not 1 + has_cond <= nf <= n:
            raise SimulationError(
                f"feed stage {nf} is outside the column", where=self.id,
                hint=f"use a stage from {1 + has_cond} to {n}"
                + (" (stage 1 is the condenser)" if has_cond else ""))

    def idle_outlets(self) -> tuple[str, ...]:
        """Outlets this configuration leaves empty; the solver does not warn about them."""
        return ("vapor_distillate",) if self.p["condenser"] == "total" else ("distillate",)

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        P = self._pressures(self.p["P"] or s.P)
        if s.total_flow <= 0.0:
            self.results.update(condenser_duty=0.0, reboiler_duty=0.0)
            return {port.name: empty_stream(pkg, s.T, P[0] if port.name != "bottoms" else P[-1])
                    for port in self.outlets}
        specs = {name: self.p[name] for name in SPECS if self.p[name] is not None}
        for name in ("distillate_rate", "bottoms_rate"):
            if name in specs and specs[name] >= s.total_flow:
                raise SimulationError(
                    f"{SPECS[name]} {specs[name] * 3.6:.4g} kmol/h is not below the feed "
                    f"{s.total_flow * 3.6:.4g} kmol/h",
                    hint="both products need some flow; a good first guess for the distillate "
                    "is the feed flow of the components lighter than the heavy key",
                )
        setup = ColumnSetup(
            stages=self.p["stages"], P=P, feeds=[_feed(self.p["feed_stage"] - 1, s)],
            condenser=self.p["condenser"], reboiler=self.p["reboiler"] != "none",
            specs=specs, efficiency=self.p["efficiency"] or 1.0,
        )
        sol = self._solve(pkg, setup)
        D = sol.distillate.sum() + sol.vapor_distillate.sum()
        B = sol.bottoms.sum()
        self.results.update(
            condenser_duty=sol.condenser_duty, reboiler_duty=sol.reboiler_duty,
            reflux_ratio=float(sol.L[0] / D) if self.p["condenser"] != "none" else None,
            boilup_ratio=float(sol.V[-1] / B) if setup.reboiler else None,
        )
        return {
            "distillate": _product(pkg, sol.distillate, P[0], sol.h[0], 0.0, sol.T[0]),
            "vapor_distillate": _product(pkg, sol.vapor_distillate, P[0], sol.Hv[0], 1.0,
                                         sol.T[0]),
            "bottoms": _product(pkg, sol.bottoms, P[-1], sol.h[-1], 0.0, sol.T[-1]),
        }


class Absorber(_StageColumn):
    type = "Absorber"
    label = "Absorber / Stripper"
    category = "Columns"
    description = "Countercurrent gas-liquid contactor on equilibrium stages, without "\
        "condenser or reboiler: liquid enters at the top, gas at the bottom. Absorbs solutes "\
        "from a gas into a solvent, or strips them out of a liquid."
    inlets = (Port("gas_in", "Gas in (bottom)"), Port("liquid_in", "Liquid in (top)"))
    outlets = (Port("gas_out", "Gas out (top)"), Port("liquid_out", "Liquid out (bottom)"))
    params = (
        Param("stages", "Number of stages", "integer", default=10, min=1,
              help="Equilibrium stages; the liquid enters stage 1, the gas the last stage."),
        Param("P", "Top pressure", "pressure", min=1.0,
              help="Defaults to the lower of the two inlet pressures."),
        DP_STAGE,
        EFFICIENCY,
    )
    result_fields = (ResultField("T_top", "Top temperature", "temperature"),
                     ResultField("T_bottom", "Bottom temperature", "temperature"),
                     ResultField("iterations", "Newton iterations"))

    def calculate(self, pkg, inlets):
        gas, liquid = single_inlet(inlets, "gas_in"), single_inlet(inlets, "liquid_in")
        P = self._pressures(self.p["P"] or min(gas.P, liquid.P))
        if gas.total_flow <= 0.0 or liquid.total_flow <= 0.0:
            # Nothing to contact (e.g. a recycle loop's first pass): streams pass through.
            self.results.update(T_top=liquid.T, T_bottom=gas.T)
            return {"gas_out": gas if gas.total_flow > 0 else empty_stream(pkg, gas.T, P[0]),
                    "liquid_out": liquid if liquid.total_flow > 0 else
                    empty_stream(pkg, liquid.T, P[-1])}
        n = self.p["stages"]
        setup = ColumnSetup(stages=n, P=P, feeds=[_feed(0, liquid), _feed(n - 1, gas)],
                            condenser="none", reboiler=False,
                            efficiency=self.p["efficiency"] or 1.0)
        sol = self._solve(pkg, setup)
        self.results.update(T_top=float(sol.T[0]), T_bottom=float(sol.T[-1]))
        return {
            "gas_out": _product(pkg, sol.vapor_distillate, P[0], sol.Hv[0], 1.0, sol.T[0]),
            "liquid_out": _product(pkg, sol.bottoms, P[-1], sol.h[-1], 0.0, sol.T[-1]),
        }
