"""Distillation and absorption columns.

The models live in ``fluxsheet.columns``: a shortcut design method and one rigorous
equilibrium-stage solver that the distillation column and the absorber/stripper configure
differently.
"""

from __future__ import annotations

import math

import numpy as np

from ..columns import (
    ColumnSetup,
    Pumparound,
    SideDraw,
    StageFeed,
    shortcut_design,
    solve_column,
)
from ..columns.saturation import dew_point
from ..columns.stages import SPECS
from ..errors import SimulationError
from ..reactions import ReactionSet
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


# Side draws, pumparounds, stage duties and reactive stages of the distillation column.
MAX_SIDE_DRAWS = 4
_STAGE = Param("stage", "Stage", "integer", min=1, help="Counted from the top (condenser = 1).")
SIDE_DRAWS = Param(
    "side_draws", "Side draws", "records",
    help=f"Liquid or vapor products from intermediate stages, at most {MAX_SIDE_DRAWS}; they "
    "leave through the outlets Side draw 1, 2, ... in this order. Give a flow, or a fraction "
    "of the stage's liquid (vapor) outflow.",
    fields=(_STAGE,
            Param("phase", "Phase", "choice", default="liquid", options=("liquid", "vapor")),
            Param("rate", "Flow", "molar_flow", min=0.0),
            Param("fraction", "Fraction of outflow", "fraction", min=0.0, max=0.99)),
)
PUMPAROUNDS = Param(
    "pumparounds", "Pumparounds", "records",
    help="Liquid drawn from one stage, cooled or heated, and returned to another (usually "
    "higher) stage. Give a flow or a fraction of the stage's liquid outflow, and a duty or a "
    "return temperature (neither: no heat exchanged).",
    fields=(Param("from_stage", "From stage", "integer", min=1),
            Param("to_stage", "To stage", "integer", min=1),
            Param("rate", "Flow", "molar_flow", min=0.0),
            Param("fraction", "Fraction of outflow", "fraction", min=0.0, max=0.99),
            Param("duty", "Duty", "power", help="Heat added; negative for a cooler."),
            Param("T", "Return temperature", "temperature", min=1.0)),
)
STAGE_DUTIES = Param(
    "stage_duties", "Side heaters and coolers", "records",
    help="Heat added to single stages (negative: removed), e.g. an intercooler.",
    fields=(_STAGE, Param("duty", "Duty", "power")),
)
COLUMN_REACTIONS = Param(
    "reactions", "Reactions", "string_list",
    help="Ids of flowsheet reactions with kinetics that run in the liquid on the reactive "
    "stages.")
REACTIVE_STAGES = Param(
    "reactive_stages", "Reactive stages", "records",
    help="Stages where the reactions run, each with its liquid hold-up (e.g. the catalyst "
    "volume wetted by liquid). Rates are those of the reactions at the stage's temperature "
    "and liquid composition.",
    fields=(_STAGE, Param("holdup", "Liquid hold-up", "volume", min=0.0)),
)

# Product rates given relative to the feed; converted to rates when the feed is known.
_FEED_RATIOS = {"distillate_to_feed": "distillate to feed ratio",
                "bottoms_to_feed": "bottoms to feed ratio"}


def _noncondensables(pkg, feeds, F: float) -> list[str]:
    """Feed components that stay gaseous in any ordinary condenser (critical point below
    250 K: hydrogen, nitrogen, methane, ...), if more than traces of them are fed."""
    total = sum(fd.flows for fd in feeds)
    return [name for name, f, Tc in zip(pkg.names, total, pkg.consts.Tcs)
            if Tc < 250.0 and f > 1e-4 * F]


def _product(pkg, flows, P: float, H: float, vf: float, T: float) -> Stream:
    """A product stream with the stage's enthalpy, so the column's energy balance closes."""
    flows = [float(f) for f in flows]
    total = sum(flows)
    if total <= 0.0:
        return empty_stream(pkg, T, P)
    if not pkg.supports_lle:
        # The stage already gives temperature, phase and composition; the same phase at the
        # same T has exactly the stage's enthalpy. A PH flash would only find T again, and on
        # wide-boiling mixtures it takes seconds.
        state = pkg.phase_state("vapor" if vf == 1.0 else "liquid", T, P,
                                [f / total for f in flows])
        return Stream(flows, state)
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
        # Inside a recycle loop the column is solved once per pass; the previous pass is the
        # best start for the next one.
        sol = solve_column(pkg, setup, start=getattr(self, "_last_solution", None))
        self._last_solution = sol
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
        "stage the reboiler. Give one specification per condenser and reboiler. Takes any "\
        "number of feeds, each on its own stage, side draws, pumparounds, side heaters and "\
        "kinetic reactive stages."
    inlets = (Port("in", "Feeds", multiple=True, separate=True),)
    outlets = (Port("distillate", "Distillate"), Port("bottoms", "Bottoms"),
               Port("vapor_distillate", "Vapor distillate"),
               *(Port(f"side{k}", f"Side draw {k}") for k in range(1, MAX_SIDE_DRAWS + 1)))
    reactions = ReactionSet()  # replaced by the solver with the flowsheet's reactions
    params = (
        Param("stages", "Number of stages", "integer", default=10, min=2,
              help="Equilibrium stages, counting the condenser and the reboiler."),
        Param("feed_stage", "Feed stage", "integer", default=5, min=1,
              help="Stage the feed enters, counted from the top (condenser = 1)."),
        Param("feed_stages", "Feed stage per stream", "stream_stages",
              help="With several feeds: the stage of each feed stream. Feeds not listed "
              "enter the feed stage above."),
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
        Param("distillate_to_feed", "Distillate to feed ratio", "fraction", min=0.0, max=1.0,
              help="D/F, molar. Follows the feed, so it suits a column inside a recycle loop, "
              "where the feed flow changes from one pass to the next."),
        Param("bottoms_to_feed", "Bottoms to feed ratio", "fraction", min=0.0, max=1.0,
              help="B/F, molar."),
        EFFICIENCY,
        SIDE_DRAWS,
        PUMPAROUNDS,
        STAGE_DUTIES,
        COLUMN_REACTIONS,
        REACTIVE_STAGES,
    )
    exclusive_groups = (("reflux_ratio", "boilup_ratio"),
                        ("distillate_rate", "bottoms_rate", "distillate_to_feed",
                         "bottoms_to_feed"))
    result_fields = (CONDENSER_DUTY, REBOILER_DUTY, ResultField("reflux_ratio", "Reflux ratio"),
                     ResultField("boilup_ratio", "Boilup ratio"),
                     ResultField("pumparound_duty", "Pumparound duty (total)", "power"),
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
        given = [label for name, label in {**SPECS, **_FEED_RATIOS}.items()
                 if self.p[name] is not None]
        for name, label in _FEED_RATIOS.items():
            if self.p[name] is not None and not 0.0 < self.p[name] < 1.0:
                raise SimulationError(f"{label} {self.p[name]:g} leaves a product empty",
                                      where=self.id, hint="use a value between 0 and 1")
        needed = has_cond + has_reb
        if len(given) != needed:
            duties = {2: "one for the condenser and one for the reboiler",
                      1: "one, for the " + ("condenser" if has_cond else "reboiler"),
                      0: "none without condenser and reboiler"}[needed]
            raise SimulationError(
                f"{'under' if len(given) < needed else 'over'}-specified: this column takes "
                f"{duties}, but has {len(given)}{': ' + ', '.join(given) if given else ''}",
                where=self.id,
                hint={2: "set one of reflux ratio or boilup ratio and one of distillate rate, "
                         "bottoms rate or their ratio to the feed; a reflux ratio and a "
                         "distillate rate is typical",
                      1: "set exactly one ratio or rate that fits this configuration",
                      0: "remove the specifications, or add a condenser or reboiler"}[needed],
            )
        n = self.p["stages"]
        if n < 1 + has_cond + has_reb:
            raise SimulationError(f"{n} stages is too few", where=self.id,
                                  hint="count the condenser and reboiler as stages; use at "
                                  f"least {2 + has_cond + has_reb}")
        stages = {"feed": self.p["feed_stage"], **(self.p["feed_stages"] or {})}
        for name, nf in stages.items():
            if not 1 + has_cond <= nf <= n:
                raise SimulationError(
                    f"feed stage {nf}{'' if name == 'feed' else f' of {name}'} is outside "
                    "the column", where=self.id,
                    hint=f"use a stage from {1 + has_cond} to {n}"
                    + (" (stage 1 is the condenser)" if has_cond else ""))
        self._check_side_streams(n, has_cond, has_reb)

    def _check_side_streams(self, n: int, has_cond: bool, has_reb: bool) -> None:
        def fail(message: str, hint: str) -> None:
            raise SimulationError(message, where=self.id, hint=hint)

        def check_stage(stage: int, what: str) -> None:
            if not 1 <= stage <= n:
                fail(f"{what}: stage {stage} is outside the column",
                     f"use a stage from 1 to {n}")

        def one_of(entry: dict, names: tuple[str, str], what: str, required: bool) -> None:
            given = [k for k in names if entry[k] is not None]
            if len(given) > 1 or (required and not given):
                fail(f"{what}: give {'either' if required else 'at most one of'} "
                     f"{' or '.join(names)}", f"set exactly {'one' if required else 'one or none'}"
                     f" of {' and '.join(names)}")

        draws = self.p["side_draws"] or []
        if len(draws) > MAX_SIDE_DRAWS:
            fail(f"{len(draws)} side draws; the column has outlets for {MAX_SIDE_DRAWS}",
                 "combine draws on the same stage, or use a Splitter after one draw")
        fractions: dict[tuple[int, str], float] = {}
        for k, d in enumerate(draws, 1):
            what = f"side draw {k}"
            if d["stage"] is None:
                fail(f"{what} has no stage", "give the stage it is drawn from")
            check_stage(d["stage"], what)
            if (has_cond and d["stage"] == 1) or (has_reb and d["stage"] == n):
                fail(f"{what} is on the {'condenser' if d['stage'] == 1 else 'reboiler'}",
                     "the distillate and bottoms are the products of condenser and reboiler; "
                     "draw from a stage in between")
            one_of(d, ("rate", "fraction"), what, True)
            key = (d["stage"], d["phase"])
            fractions[key] = fractions.get(key, 0.0) + (d["fraction"] or 0.0)
        for k, pa in enumerate(self.p["pumparounds"] or [], 1):
            what = f"pumparound {k}"
            for name in ("from_stage", "to_stage"):
                if pa[name] is None:
                    fail(f"{what} has no {name.replace('_', ' ')}",
                         "give the stage it is drawn from and the one it returns to")
                check_stage(pa[name], what)
            if pa["from_stage"] == pa["to_stage"]:
                fail(f"{what} returns to the stage it is drawn from",
                     "use a side heater or cooler on that stage instead")
            one_of(pa, ("rate", "fraction"), what, True)
            one_of(pa, ("duty", "T"), what, False)
            if pa["from_stage"] == 1 and self.p["condenser"] == "total" and pa["fraction"]:
                fail(f"{what} takes a fraction of the condenser's outflow, which the "
                     "distillate specification also sets", "give the pumparound a flow")
            if pa["rate"] == 0.0 and (pa["duty"] or pa["T"] is not None):
                fail(f"{what} has no flow to heat or cool",
                     "give it a flow, or use a side heater or cooler")
            key = (pa["from_stage"], "liquid")
            fractions[key] = fractions.get(key, 0.0) + (pa["fraction"] or 0.0)
        for (stage, phase), total in fractions.items():
            if total >= 0.99:
                fail(f"the {phase} withdrawn from stage {stage} adds up to a fraction of "
                     f"{total:g}", "leave some flow on the stage; fractions must sum below 1")
        for k, sd in enumerate(self.p["stage_duties"] or [], 1):
            if sd["stage"] is None or sd["duty"] is None:
                fail(f"side heater/cooler {k} needs a stage and a duty", "fill in both")
            check_stage(sd["stage"], f"side heater/cooler {k}")
            if (has_cond and sd["stage"] == 1) or (has_reb and sd["stage"] == n):
                fail(f"side heater/cooler {k} is on the "
                     f"{'condenser' if sd['stage'] == 1 else 'reboiler'}, whose duty the "
                     "specifications already set", "put it on a stage in between")
        reactive = self.p["reactive_stages"] or []
        for k, rs in enumerate(reactive, 1):
            if rs["stage"] is None or rs["holdup"] is None:
                fail(f"reactive stage {k} needs a stage and a liquid hold-up", "fill in both")
            check_stage(rs["stage"], f"reactive stage {k}")
        if reactive and not self.p["reactions"]:
            fail("reactive stages but no reactions",
                 "list the ids of the reactions that run on them")
        if self.p["reactions"] and not reactive:
            fail("reactions but no reactive stages",
                 "add the stages where they run, each with its liquid hold-up")

    def _feed_stage(self, stream_id: str | None) -> int:
        """Stage index (from 0) where a feed stream enters."""
        return (self.p["feed_stages"] or {}).get(stream_id, self.p["feed_stage"]) - 1

    def idle_outlets(self) -> tuple[str, ...]:
        """Outlets this configuration leaves empty; the solver does not warn about them."""
        top = ("vapor_distillate",) if self.p["condenser"] == "total" else ("distillate",)
        used = len(self.p["side_draws"] or [])
        return top + tuple(f"side{k}" for k in range(used + 1, MAX_SIDE_DRAWS + 1))

    def _side_setup(self, pkg) -> dict:
        """The side draws, pumparounds, stage duties and reactions as ColumnSetup fields."""
        draws = [SideDraw(d["stage"] - 1, d["phase"], d["rate"], d["fraction"])
                 for d in self.p["side_draws"] or []]
        pas = [Pumparound(pa["from_stage"] - 1, pa["to_stage"] - 1, pa["rate"], pa["fraction"],
                          pa["duty"], pa["T"]) for pa in self.p["pumparounds"] or []]
        duties: dict[int, float] = {}
        for sd in self.p["stage_duties"] or []:
            duties[sd["stage"] - 1] = duties.get(sd["stage"] - 1, 0.0) + sd["duty"]
        out = {"draws": draws, "pumparounds": pas, "stage_duties": duties}
        if not self.p["reactions"]:
            return out
        rxns = self.reactions.select(self.p["reactions"], self.id)
        for r in rxns:
            kin = r.kinetics
            if kin is None:
                raise SimulationError(
                    f"reaction {r.id} has no kinetics, which reactive stages need",
                    where=self.id, hint=f"add a rate law to {r.id} in the flowsheet's reactions")
            if kin.basis == "partial_pressure":
                raise SimulationError(
                    f"reaction {r.id} has a partial-pressure rate law, but reactive stages "
                    "react in the liquid", where=self.id, hint="use the concentration basis")
            if kin.equilibrium_limited and r.keq is None:
                raise SimulationError(
                    f"reaction {r.id} is equilibrium-limited in the liquid but has no user "
                    "equilibrium constant", where=self.id,
                    hint="add keq (ln K = A + B/T + C ln T + D T, K on mole fractions) to the "
                    "reaction; the built-in K refers to ideal gases")
        involved = sorted({i for r in rxns for i in r.involved})
        missing = [pkg.names[i] for i in involved if pkg.consts.Hfgs[i] is None]
        if missing:
            raise SimulationError(
                f"no formation enthalpy in the database for {', '.join(missing)}; the heat of "
                "reaction on the reactive stages needs it", where=self.id,
                hint="choose compounds that have formation data")
        holdup: dict[int, float] = {}
        for rs in self.p["reactive_stages"]:
            holdup[rs["stage"] - 1] = holdup.get(rs["stage"] - 1, 0.0) + rs["holdup"]
        hf = np.array([pkg.consts.Hfgs[i] if i in involved else 0.0 for i in range(pkg.n)])
        return {**out, "reactions": rxns, "holdup": holdup, "formation_enthalpies": hf}

    def calculate(self, pkg, inlets):
        streams = inlets["in"]
        ids = self.inlet_ids.get("in", [])
        F = sum(s.total_flow for s in streams)
        P = self._pressures(self.p["P"] or min(s.P for s in streams))
        if F <= 0.0:
            self.results.update(condenser_duty=0.0, reboiler_duty=0.0)
            T = streams[0].T
            return {port.name: empty_stream(pkg, T, P[0] if port.name != "bottoms" else P[-1])
                    for port in self.outlets}
        feeds = [_feed(self._feed_stage(ids[i] if i < len(ids) else None), s)
                 for i, s in enumerate(streams) if s.total_flow > 0.0]
        specs = {name: self.p[name] for name in SPECS if self.p[name] is not None}
        for name, rate in (("distillate_to_feed", "distillate_rate"),
                           ("bottoms_to_feed", "bottoms_rate")):
            if self.p[name] is not None:
                specs[rate] = self.p[name] * F
        for name in ("distillate_rate", "bottoms_rate"):
            if name in specs and specs[name] >= F:
                raise SimulationError(
                    f"{SPECS[name]} {specs[name] * 3.6:.4g} kmol/h is not below the feed "
                    f"{F * 3.6:.4g} kmol/h",
                    hint="both products need some flow; a good first guess for the distillate "
                    "is the feed flow of the components lighter than the heavy key",
                )
        setup = ColumnSetup(
            stages=self.p["stages"], P=P, feeds=feeds,
            condenser=self.p["condenser"], reboiler=self.p["reboiler"] != "none",
            specs=specs, efficiency=self.p["efficiency"] or 1.0, **self._side_setup(pkg),
        )
        try:
            sol = self._solve(pkg, setup)
        except SimulationError as exc:
            gases = _noncondensables(pkg, feeds, F)
            fed = sum(fd.flows for fd in feeds)
            gas_flow = sum(float(fed[pkg.names.index(g)]) for g in gases)
            # Only a fixed product rate tells how much room the distillate has.
            rate_given = "distillate_rate" in specs or "bottoms_rate" in specs
            D = specs.get("distillate_rate", F - specs.get("bottoms_rate", 0.0))
            cause = None
            if self.p["condenser"] == "total" and gases:
                cause = (f"the feed carries {' and '.join(gases)}, which a total condenser "
                         "cannot condense; use a partial condenser, or remove them upstream "
                         "(a flash or stabilizer column)")
            elif gases and rate_given and D < 1.05 * gas_flow:
                cause = (f"the distillate {D * 3.6:.4g} kmol/h leaves no room for the "
                         f"{' and '.join(gases)} fed ({gas_flow * 3.6:.4g} kmol/h), which "
                         "cannot leave with the bottoms; raise the distillate. Inside a "
                         "recycle loop the first passes may carry more of them than the "
                         "converged plant")
            if cause:
                # The light gases are the likely cause, not a certain one: keep the original.
                exc.hint = cause + (f". Otherwise: {exc.hint}" if exc.hint else "")
            raise
        D = sol.distillate.sum() + sol.vapor_distillate.sum()
        B = sol.bottoms.sum()
        self.results.update(
            condenser_duty=sol.condenser_duty, reboiler_duty=sol.reboiler_duty,
            reflux_ratio=float(sol.L[0] / D) if self.p["condenser"] != "none" else None,
            boilup_ratio=float(sol.V[-1] / B) if setup.reboiler else None,
        )
        out = {
            "distillate": _product(pkg, sol.distillate, P[0], sol.h[0], 0.0, sol.T[0]),
            "vapor_distillate": _product(pkg, sol.vapor_distillate, P[0], sol.Hv[0], 1.0,
                                         sol.T[0]),
            "bottoms": _product(pkg, sol.bottoms, P[-1], sol.h[-1], 0.0, sol.T[-1]),
        }
        for k, (d, flows) in enumerate(zip(setup.draws, sol.side_draws), 1):
            liquid = d.phase == "liquid"
            j = d.stage
            out[f"side{k}"] = _product(pkg, flows, P[j], sol.h[j] if liquid else sol.Hv[j],
                                       0.0 if liquid else 1.0, sol.T[j])
        for k in range(len(setup.draws) + 1, MAX_SIDE_DRAWS + 1):
            out[f"side{k}"] = empty_stream(pkg, sol.T[0], P[0])
        self._side_results(pkg, setup, sol)
        return out

    def _side_results(self, pkg, setup: ColumnSetup, sol) -> None:
        if setup.pumparounds:
            pas = []
            for pa, flows, duty, H in zip(setup.pumparounds, sol.pumparound_flows,
                                          sol.pumparound_duties, sol.pumparound_return_H):
                total = float(flows.sum())
                T_ret = pa.T
                if T_ret is None:
                    T_ret = float(sol.T[pa.source])
                    if total > 1e-9 * float(sol.L.max()):
                        try:
                            T_ret = pkg.flash(list(flows / total), P=float(setup.P[pa.target]),
                                              H=H / total).T
                        except SimulationError:
                            T_ret = None  # the column has converged; only the report lacks it
                pas.append({"from_stage": pa.source + 1, "to_stage": pa.target + 1,
                            "flow": total, "duty": duty, "T_return": float(T_ret)})
            self.results["pumparounds"] = pas
            self.results["pumparound_duty"] = float(sum(sol.pumparound_duties))
        if sol.generation is not None:
            stages = list(range(1, setup.stages + 1))
            self.add_profile("reaction", "Formed by reactions", "Stage", stages,
                             {name: ("molar_flow", sol.generation[:, i].tolist())
                              for i, name in enumerate(pkg.names)})


class Absorber(_StageColumn):
    type = "Absorber"
    label = "Absorber / Stripper"
    category = "Columns"
    description = "Countercurrent gas-liquid contactor on equilibrium stages, without "\
        "condenser or reboiler: liquid enters at the top, gas at the bottom. Absorbs solutes "\
        "from a gas into a solvent, or strips them out of a liquid."
    inlets = (Port("gas_in", "Gas in (bottom)", multiple=True),
              Port("liquid_in", "Liquid in (top)", multiple=True))
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
