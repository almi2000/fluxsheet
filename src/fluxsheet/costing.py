"""Equipment costing: purchased and bare-module costs, utility costs, total annualized cost.

A study-grade (order-of-magnitude, about ±30 %) estimate after the method of Turton et al.,
*Analysis, Synthesis and Design of Chemical Processes*, Appendix A. The purchased cost of a
piece of equipment of size A (volume, area or power) in carbon steel at ambient pressure is

    log10 Cp0 = K1 + K2 log10 A + K3 (log10 A)^2

in dollars of September 2001 (CEPCI 397), scaled to the requested year by CEPCI / 397. The
bare-module cost adds installation (piping, instruments, labour, ...) and corrects for the
material and the pressure: CBM = Cp0 (B1 + B2 FM FP), or Cp0 FBM for equipment without a
pressure correlation. Sizes above a correlation's range are split into identical units in
parallel; below it the cost follows the six-tenths rule from the smallest size, so it keeps
falling smoothly (with a note: the correlation does not cover it).

Sizes come from the solved flowsheet:
- distillation columns and absorbers: the diameter from the largest vapor load over the
  stages at an F-factor u sqrt(rho_V) = 2.0 Pa^0.5 (some 70–80 % of flood for sieve trays at
  0.6 m spacing and moderate pressure; the vapor an ideal gas, the liquid density ignored),
  the net area 88 % of the cross section; the
  height 0.6 m per tray plus 3 m for sump and disengagement. Equilibrium stages become trays
  through ``tray_efficiency`` unless the column has a Murphree efficiency already. Condensers
  are fixed-tube-sheet exchangers, reboilers kettles.
- flash drums: a vertical vessel with L/D = 3, holding five minutes of liquid when half
  full, and at least wide enough for the vapor at the Souders-Brown velocity
  (K = 0.07 m/s, no demister).
- heaters, coolers, the heat duty of drums, condensers and reboilers: an exchanger with
  A = |Q| / (U LMTD). The LMTD is the unit's own, with a utility, else that against the
  default utility chosen below. HeatExchanger: A = UA / U.
- pumps by shaft power (motor not included), compressors by shaft power (Turton's
  correlation is for fluid power, so this errs high by the efficiency) plus an electric
  drive, stirred reactors and tanks by volume.

Operating cost per year: each heating or cooling duty times a price in $/GJ, and the
electricity of pumps and compressors. Without a price, the cheapest utility hot or cold
enough for the duty is chosen (prices of Turton's Table 8.3, 4th ed.; replace them with
your own):

    heating: LP steam 160 °C 4.54, MP steam 184 °C 4.77, HP steam 254 °C 5.66 $/GJ,
             above that fired heat (natural gas 3.16 $/GJ at 80 % efficiency; the furnace
             itself is not costed),
    cooling: cooling water 30 → 40 °C 0.354, chilled water 5 → 15 °C 4.43,
             refrigerant −20 °C 7.89, −50 °C 13.11 $/GJ,

each with at least 10 K between utility and process at both ends. These prices are of
about 2016 while the equipment is escalated to the chosen CEPCI; set current prices for a
real estimate. A negative duty is cooling, a positive one heating, whatever the unit. The
total annualized cost is

    TAC = operating cost + f · capital,  capital = 1.18 Σ CBM (total module cost:
    contingency and fee), f = i (1 + i)^n / ((1 + i)^n − 1)

(capital recovery factor; ``annualization_factor`` overrides it).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .errors import SimulationError
from .flowsheet import FlowsheetModel
from .quantities import UnitError, to_si
from .solver import SolveResult
from .specs import UNIT_COSTS
from .study_models import CostingModel
from .utilities import lmtd

BASE_CEPCI = 397.0  # September 2001, the basis of Turton's correlations
R_GAS = 8.314462618
P_ATM = 101325.0
APPROACH = 10.0  # K, least temperature difference between a default utility and the process
MIN_LMTD = 5.0  # K, floor for an estimated exchanger LMTD
T_SLACK = 0.1  # K, so a process exactly at a utility's limit still qualifies
F_FACTOR = 2.0  # Pa^0.5, design vapor load of a tray column
TRAY_SPACING = 0.6  # m
COLUMN_EXTRA_HEIGHT = 3.0  # m, sump and vapor disengagement
NET_AREA = 0.88  # net (vapor flow) over total cross section of a tray column
SOUDERS_BROWN_K = 0.07  # m/s, vertical drum without demister
HOLDUP = 300.0  # s of liquid, drum half full
TOTAL_MODULE = 1.18  # contingency (15 %) and fee (3 %) on the bare-module cost

# Results a costed unit reports, readable by design specs, studies and optimizations.
UNIT_FIELDS = UNIT_COSTS
# Flowsheet totals, referenced as {"property": "tac"}.
SUMMARY_FIELDS = {
    "tac": "Total annualized cost ($/yr)",
    "capital_cost": "Capital cost ($)",
    "operating_cost": "Operating cost ($/yr)",
    "bare_module_cost": "Bare module cost ($)",
    "purchased_cost": "Purchased cost ($)",
}


# ---- correlations --------------------------------------------------------------------------

@dataclass(frozen=True)
class Correlation:
    """Turton's purchased-cost correlation of one kind of equipment (2001 $, CS)."""

    label: str
    K: tuple[float, float, float]
    size: str  # "m3", "m2" or "kW"
    lo: float
    hi: float

    def cost(self, A: float) -> tuple[float, int]:
        """Purchased cost of size ``A`` in base-year dollars, and the number in parallel."""
        n = max(1, math.ceil(A / self.hi - 1e-9))
        a = A / n
        lg = math.log10(max(a, self.lo))
        K1, K2, K3 = self.K
        cost = 10.0 ** (K1 + K2 * lg + K3 * lg * lg)
        if a < self.lo:  # six-tenths rule below the range
            cost *= (a / self.lo) ** 0.6
        return n * cost, n


VESSEL = Correlation("vertical vessel", (3.4974, 0.4485, 0.1074), "m3", 0.3, 520.0)
TRAYS = Correlation("sieve trays", (2.9949, 0.4465, 0.3961), "m2", 0.07, 12.3)
EXCHANGER = Correlation("fixed-tube-sheet exchanger", (4.3247, -0.3030, 0.1634), "m2", 10.0, 1000.0)
KETTLE = Correlation("kettle reboiler", (4.4646, -0.5277, 0.3955), "m2", 10.0, 100.0)
PUMP = Correlation("centrifugal pump", (3.3892, 0.0536, 0.1538), "kW", 1.0, 300.0)
COMPRESSOR = Correlation("centrifugal compressor", (2.2897, 1.3604, -0.1027), "kW", 450.0, 3000.0)
DRIVE = Correlation("electric drive", (2.4604, 1.4191, -0.1798), "kW", 75.0, 2600.0)
STIRRED = Correlation("jacketed agitated vessel", (4.1052, 0.5320, -0.0005), "m3", 0.1, 35.0)

# Material factors FM (Turton Fig. A.18, Table A.3), carbon steel / stainless steel.
FM_VESSEL = {"CS": 1.0, "SS": 3.1}
FM_EXCHANGER = {"CS": 1.0, "SS": 2.73}  # shell and tubes of the same material
FM_PUMP = {"CS": 1.6, "SS": 2.3}
FBM_TRAYS = {"CS": 1.0, "SS": 1.83}
FBM_COMPRESSOR = {"CS": 2.7, "SS": 5.8}
FBM_DRIVE = 1.5
FBM_STIRRED = 4.0  # Turton gives no material factor for it


def vessel_pressure_factor(P: float, D: float) -> float:
    """FP of a vessel at absolute pressure ``P`` (Pa) and diameter ``D`` (m), Turton eq. A.2:
    the wall thickness for the pressure over the 6.3 mm minimum, with 3.15 mm corrosion."""
    barg = (P - P_ATM) / 1e5
    if barg < -0.5:
        return 1.25  # vacuum vessel
    fp = ((barg + 1.0) * D / (2.0 * (850.0 - 0.6 * (barg + 1.0))) + 0.00315) / 0.0063
    return max(fp, 1.0)


def _log_factor(P: float, lo: float, hi: float, C: tuple[float, float, float]) -> float:
    barg = (P - P_ATM) / 1e5
    if barg <= lo:
        return 1.0
    lg = math.log10(min(barg, hi))
    return 10.0 ** (C[0] + C[1] * lg + C[2] * lg * lg)


def exchanger_pressure_factor(P: float) -> float:
    """FP of a shell-and-tube exchanger, both sides at ``P`` (Pa); 1 below 5 barg."""
    return _log_factor(P, 5.0, 140.0, (0.03881, -0.11272, 0.08183))


def pump_pressure_factor(P: float) -> float:
    """FP of a centrifugal pump discharging at ``P`` (Pa); 1 below 10 barg."""
    return _log_factor(P, 10.0, 100.0, (-0.3935, 0.3957, -0.00226))


# ---- utilities ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DefaultUtility:
    label: str
    T_in: float  # K
    T_out: float  # K
    price: float  # $/GJ

    def can(self, Q: float, T_proc: tuple[float, float]) -> bool:
        """Whether it keeps ``APPROACH`` from the process at both ends, countercurrent."""
        cold, hot = min(T_proc), max(T_proc)
        gap = APPROACH - T_SLACK
        if Q > 0:
            return self.T_in >= hot + gap and self.T_out >= cold + gap
        return self.T_in <= cold - gap and self.T_out <= hot - gap


HEATING = (
    DefaultUtility("LP steam (160 °C)", 433.15, 433.15, 4.54),
    DefaultUtility("MP steam (184 °C)", 457.15, 457.15, 4.77),
    DefaultUtility("HP steam (254 °C)", 527.15, 527.15, 5.66),
)
FIRED = DefaultUtility("fired heat (natural gas)", math.inf, math.inf, 3.16 / 0.8)
COOLING = (
    DefaultUtility("cooling water (30 → 40 °C)", 303.15, 313.15, 0.354),
    DefaultUtility("chilled water (5 → 15 °C)", 278.15, 288.15, 4.43),
    DefaultUtility("refrigerant (−20 °C)", 253.15, 253.15, 7.89),
    DefaultUtility("refrigerant (−50 °C)", 223.15, 223.15, 13.11),
)


def default_utility(Q: float, T_proc: tuple[float, float]) -> tuple[DefaultUtility, str | None]:
    """The cheapest default utility for a duty ``Q`` (W, > 0 heating) on a process going from
    ``T_proc[0]`` to ``T_proc[1]``, and a note when none really fits."""
    if Q > 0:
        for ut in HEATING:
            if ut.can(Q, T_proc):
                return ut, None
        return FIRED, None
    for ut in COOLING:
        if ut.can(Q, T_proc):
            return ut, None
    return COOLING[-1], (f"cooling down to {min(T_proc) - 273.15:.0f} °C needs a colder "
                         "refrigerant than −50 °C; priced as one at −50 °C")


def default_lmtd(ut: DefaultUtility, Q: float, T_proc: tuple[float, float]) -> float:
    """Countercurrent LMTD of a process going from T_proc[0] to T_proc[1] against ``ut``."""
    T1, T2 = T_proc
    d = lmtd(ut.T_out - T1, ut.T_in - T2) if Q > 0 else lmtd(T1 - ut.T_out, T2 - ut.T_in)
    return max(d or MIN_LMTD, MIN_LMTD)


def annualization_factor(c: CostingModel) -> float:
    """Capital recovery factor i(1+i)^n / ((1+i)^n - 1), or the given factor."""
    if c.annualization_factor is not None:
        return c.annualization_factor
    i, n = c.interest_rate, c.plant_life
    return 1.0 / n if i == 0.0 else i * (1 + i) ** n / ((1 + i) ** n - 1)


# ---- estimating ----------------------------------------------------------------------------

@dataclass
class _UnitCost:
    items: list[dict] = field(default_factory=list)
    energy: list[dict] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


class _Estimate:
    def __init__(self, model: FlowsheetModel, result: SolveResult, packages):
        self.model, self.result, self.packages = model, result, packages
        self.c: CostingModel = model.costing  # type: ignore[assignment]
        self.index = self.c.cepci / BASE_CEPCI
        self.seconds = self.c.hours_per_year * 3600.0
        try:
            self.U = to_si(self.c.U, "heat_transfer_coefficient")
        except (UnitError, ValueError) as exc:
            raise SimulationError(f"costing U: {exc}", where="costing") from exc
        if self.U <= 0:
            raise SimulationError("costing U must be positive", where="costing")
        unit_ids = {u.id for u in model.units}
        for uid in self.c.units:
            if uid not in unit_ids:
                raise SimulationError(f"costing names unit {uid!r}, which does not exist",
                                      where="costing", hint="remove it from costing.units")
        utility_ids = {u.id for u in model.utilities or ()}
        for uid in self.c.utility_prices:
            if uid not in utility_ids:
                raise SimulationError(f"costing prices utility {uid!r}, which does not exist",
                                      where="costing", hint="remove it from utility_prices")


    # -- unit context --

    def material(self, uid: str) -> str:
        own = self.c.units.get(uid)
        return own.material if own and own.material else self.c.material

    def unit_U(self, uid: str) -> float:
        own = self.c.units.get(uid)
        if own is None or own.U is None:
            return self.U
        try:
            U = to_si(own.U, "heat_transfer_coefficient")
        except (UnitError, ValueError) as exc:
            raise SimulationError(f"costing U: {exc}", where=uid) from exc
        if U <= 0:
            raise SimulationError("costing U must be positive", where=uid)
        return U

    def streams(self, uid: str) -> tuple[list[dict], dict[str, dict]]:
        """The unit's inlet streams and its outlets by port."""
        rs = self.result.streams
        ins = [rs[s.id] for s in self.model.streams if s.target and
               s.target.split(".")[0] == uid and s.id in rs]
        outs = {s.source.split(".")[1]: rs[s.id] for s in self.model.streams
                if s.source.split(".")[0] == uid and s.id in rs}
        for key, s in rs.items():  # unconnected outlets appear as "UNIT.port"
            if key.startswith(uid + ".") and key.count(".") == 1:
                outs.setdefault(key.split(".")[1], s)
        return ins, outs

    def item(self, cost: _UnitCost, what: str, corr: Correlation, size: float,
             bare: float | tuple[float, float, float], basis: str) -> None:
        """Add a piece of equipment; ``bare`` is FBM or (B1, B2, FM·FP)."""
        base, n = corr.cost(size)
        fbm = bare if isinstance(bare, float) else bare[0] + bare[1] * bare[2]
        cost.items.append({
            "item": what, "equipment": corr.label, "size": size, "size_unit": corr.size,
            "parallel": n, "purchased_cost": base * self.index,
            "bare_module_cost": base * fbm * self.index, "basis": basis,
        })
        if size < corr.lo:
            cost.notes.append(f"{what}: {size:.3g} {corr.size} is below the correlation's "
                              f"range ({corr.lo:g} {corr.size}); scaled down by the "
                              "six-tenths rule")

    def heat(self, cost: _UnitCost, uid: str, what: str, Q: float | None,
             T_proc: tuple[float, float], P: float, *, exchanger: Correlation | None = EXCHANGER,
             ures: dict | None = None) -> None:
        """Utility cost of duty ``Q`` (W) and, unless ``exchanger`` is None, its exchanger."""
        if Q is None or not math.isfinite(Q) or abs(Q) < 1.0:
            return
        ures = ures or {}
        utility_id = self.model_unit(uid).params.get("utility")
        ut, note = default_utility(Q, T_proc)
        if note:
            cost.notes.append(note)
        if utility_id is not None and utility_id in self.c.utility_prices:
            price, label = self.c.utility_prices[utility_id], f"utility {utility_id}"
        elif utility_id is not None:
            price, label = ut.price, f"utility {utility_id} (priced as {ut.label})"
        else:
            override = self.c.heating_price if Q > 0 else self.c.cooling_price
            price, label = (override, "heating" if Q > 0 else "cooling") if override \
                is not None else (ut.price, ut.label)
        cost.energy.append({
            "item": what, "kind": "heating" if Q > 0 else "cooling", "utility": label,
            "duty": Q, "price": price, "unit_price": "$/GJ",
            "cost": abs(Q) * self.seconds * price * 1e-9,
        })
        if exchanger is None:
            return
        if ut is FIRED and utility_id is None:
            cost.notes.append(f"{what}: heating above {HEATING[-1].T_in - 273.15 - APPROACH:.0f}"
                              " °C needs a fired heater, which is not costed")
            return
        area = ures.get("heat_transfer_area")
        if isinstance(area, (int, float)) and area > 0:
            basis = "area from the unit's U and LMTD"
        else:
            U = self.unit_U(uid)
            dT = ures.get("LMTD")
            if not (isinstance(dT, (int, float)) and dT > 0):
                dT = default_lmtd(ut, Q, T_proc)
                src = f"against {ut.label}"
            else:
                src = "the unit's"
            area = abs(Q) / (U * dT)
            basis = f"A = |Q| / (U LMTD), U {U:.0f} W/m²K, LMTD {dT:.1f} K ({src})"
        fp = exchanger_pressure_factor(P)
        self.item(cost, what, exchanger, area,
                  (1.63, 1.66, FM_EXCHANGER[self.material(uid)] * fp), basis)

    def electricity(self, cost: _UnitCost, what: str, power: float) -> None:
        if power > 0:
            cost.energy.append({
                "item": what, "kind": "electricity", "utility": "electricity",
                "duty": power, "price": self.c.electricity_price, "unit_price": "$/kWh",
                "cost": power / 1e3 * self.c.hours_per_year * self.c.electricity_price,
            })

    def model_unit(self, uid: str):
        return next(u for u in self.model.units if u.id == uid)

    def param(self, uid: str, name: str, quantity: str) -> float | None:
        raw = self.model_unit(uid).params.get(name)
        if raw is None:
            return None
        try:
            return to_si(raw, quantity)
        except (UnitError, ValueError):
            return None

    # -- unit kinds --

    def unit(self, uid: str, utype: str, ures: dict) -> _UnitCost | None:
        handler = _HANDLERS.get(utype)
        cost = _UnitCost()
        if handler is not None:
            handler(self, cost, uid, ures)
        elif isinstance(ures.get("duty"), (int, float)):  # utilities only
            ins, outs = self.streams(uid)
            self.heat(cost, uid, "duty", ures["duty"], _span(ins, outs), _pmax(ins, outs),
                      exchanger=None, ures=ures)
            if ures["duty"] and abs(ures["duty"]) >= 1.0:
                cost.notes.append("equipment not costed; utilities only")
        else:
            return None
        return cost if cost.items or cost.energy else None


def _span(ins: list[dict], outs: dict[str, dict]) -> tuple[float, float]:
    """Process temperatures at the inlet and the outlet of a unit."""
    out = [s for s in outs.values() if s["molar_flow"] > 0] or list(outs.values())
    T_in = ins[0]["T"] if ins else out[0]["T"]
    T_out = out[0]["T"] if out else T_in
    return T_in, T_out


def _pmax(ins: list[dict], outs: dict[str, dict]) -> float:
    return max([s["P"] for s in ins] + [s["P"] for s in outs.values()] + [P_ATM])


def _heater(e: _Estimate, cost: _UnitCost, uid: str, ures: dict) -> None:
    ins, outs = e.streams(uid)
    e.heat(cost, uid, "exchanger", ures.get("duty"), _span(ins, outs), _pmax(ins, outs),
           ures=ures)


def _exchanger(e: _Estimate, cost: _UnitCost, uid: str, ures: dict) -> None:
    UA = ures.get("UA")
    if not isinstance(UA, (int, float)) or UA <= 0:
        return
    ins, outs = e.streams(uid)
    U = e.unit_U(uid)
    fp = exchanger_pressure_factor(_pmax(ins, outs))
    e.item(cost, "exchanger", EXCHANGER, UA / U,
           (1.63, 1.66, FM_EXCHANGER[e.material(uid)] * fp), f"A = UA / U, U {U:.0f} W/m²K")


def _drum(e: _Estimate, cost: _UnitCost, uid: str, ures: dict) -> None:
    ins, outs = e.streams(uid)
    V_given = e.param(uid, "V", "volume")
    if V_given:  # a stirred tank of the given volume
        e.item(cost, "vessel", STIRRED, V_given, FBM_STIRRED, "given volume")
    else:
        vapor = outs.get("vapor")
        liquids = [s for port, s in outs.items() if port != "vapor" and s["molar_flow"] > 0]
        P = _pmax(ins, outs)
        T = (vapor or (liquids or ins)[0])["T"]
        pkg = e.packages(uid)
        mws = pkg.mws
        Q_liq, rho_L = 0.0, 800.0
        for s in liquids:
            zs = [s["mole_fractions"][c] for c in e.model.components]
            Vm = pkg.flash(zs, P=s["P"], VF=0.0).V  # saturated liquid molar volume
            Q_liq += s["molar_flow"] * Vm
            rho_L = s["mass_flow"] / max(s["molar_flow"] * Vm, 1e-30)
        V = 2.0 * Q_liq * HOLDUP
        basis = "5 min liquid at half full, L/D 3"
        if vapor is not None and vapor["molar_flow"] > 0:
            M = sum(vapor["mole_fractions"][c] * mw for c, mw in zip(e.model.components, mws))
            rho_V = P * M / 1e3 / (R_GAS * T)
            u = SOUDERS_BROWN_K * math.sqrt(max(rho_L - rho_V, 0.0) / rho_V)
            D_vap = math.sqrt(4.0 * vapor["mass_flow"] / (rho_V * u * math.pi)) if u > 0 else 0
            V_vap = math.pi / 4.0 * D_vap ** 2 * 3.0 * D_vap
            if V_vap > V:
                V, basis = V_vap, "vapor at the Souders-Brown velocity, L/D 3"
        if V <= 0:
            return
        D = (4.0 * V / (3.0 * math.pi)) ** (1.0 / 3.0)
        fp = vessel_pressure_factor(P, D)
        e.item(cost, "vessel", VESSEL, V, (2.25, 1.82, FM_VESSEL[e.material(uid)] * fp),
               f"{basis}: D {D:.2f} m")
    if isinstance(ures.get("duty"), (int, float)):
        e.heat(cost, uid, "heat exchange", ures["duty"], _span(ins, outs), _pmax(ins, outs),
               ures=ures)


def _reactor(e: _Estimate, cost: _UnitCost, uid: str, ures: dict) -> None:
    ins, outs = e.streams(uid)
    V = e.param(uid, "V", "volume")
    if V is None:
        L, d = e.param(uid, "length", "length"), e.param(uid, "diameter", "length")
        V = math.pi / 4.0 * d * d * L if L and d else None
    if V:
        e.item(cost, "vessel", STIRRED, V, FBM_STIRRED, "reactor volume")
    else:
        cost.notes.append("no volume: the vessel is not costed")
    # A jacketed vessel's heat transfer area is part of its cost.
    e.heat(cost, uid, "heat removal" if (ures.get("duty") or 0) < 0 else "heating",
           ures.get("duty"), _span(ins, outs), _pmax(ins, outs), exchanger=None, ures=ures)


def _pump(e: _Estimate, cost: _UnitCost, uid: str, ures: dict) -> None:
    power = ures.get("power") or 0.0
    if power <= 0:
        return
    ins, outs = e.streams(uid)
    fp = pump_pressure_factor(_pmax(ins, outs))
    e.item(cost, "pump", PUMP, power / 1e3, (1.89, 1.35, FM_PUMP[e.material(uid)] * fp),
           "shaft power; motor not included")
    e.electricity(cost, "pump", power)


def _compressor(e: _Estimate, cost: _UnitCost, uid: str, ures: dict) -> None:
    power = ures.get("power") or 0.0
    if power <= 0:
        return
    e.item(cost, "compressor", COMPRESSOR, power / 1e3, FBM_COMPRESSOR[e.material(uid)],
           "shaft power")
    e.item(cost, "drive", DRIVE, power / 1e3, FBM_DRIVE, "electric, explosion proof")
    e.electricity(cost, "compressor", power)
    if ures.get("cooler_duty"):  # multistage compressor: its intercoolers as one exchanger
        ins, outs = e.streams(uid)
        T_c = e.param(uid, "T_cooler", "temperature") or _span(ins, outs)[1]
        T_hot = ures.get("T_max") or T_c
        e.heat(cost, uid, "coolers", ures["cooler_duty"], (T_hot, T_c), _pmax(ins, outs))


def _column(e: _Estimate, cost: _UnitCost, uid: str, ures: dict) -> None:
    _, outs = e.streams(uid)
    profiles = ures.get("profiles") or {}
    comps = e.model.components
    mws = dict(zip(comps, e.packages(uid).mws))
    um = e.model_unit(uid)
    utype = um.type
    P_top = min(s["P"] for s in outs.values()) if outs else P_ATM
    P_bot = _pmax([], outs)
    if "flows" in profiles:  # rigorous: the largest vapor load over the stages
        series = {s["label"]: s["values"] for s in profiles["flows"]["series"]}
        Ts = profiles["temperature"]["series"][0]["values"]
        ys = {s["label"]: s["values"] for s in profiles["vapor_composition"]["series"]}
        V = series["Vapor V"]
        n = len(V)
        loads = []
        for j in range(n):
            P = P_top + (P_bot - P_top) * j / max(n - 1, 1)
            M = sum(ys[c][j] * mws[c] for c in comps) / 1e3
            rho = P * M / (R_GAS * Ts[j])
            loads.append(V[j] * M / math.sqrt(rho) if rho > 0 else 0.0)
        load = max(loads)
        condenser = um.params.get("condenser", "total") != "none" and utype != "Absorber"
        reboiler = um.params.get("reboiler", "kettle") != "none" and utype != "Absorber"
        stages = n
        T_top, T_top2, T_bot, T_bot2 = Ts[0], Ts[min(1, n - 1)], Ts[-1], Ts[max(n - 2, 0)]
        murphree = um.params.get("efficiency") is not None
    else:  # shortcut: the top vapor (R + 1) D at the distillate's conditions
        dist = outs.get("distillate")
        stages = ures.get("stages")
        R = ures.get("reflux_ratio")
        if dist is None or not stages or R is None:
            cost.notes.append("no stage count or reflux: the column is not costed")
            return
        M = sum(dist["mole_fractions"][c] * mws[c] for c in comps) / 1e3
        rho = dist["P"] * M / (R_GAS * dist["T"])
        load = (R + 1.0) * dist["molar_flow"] * M / math.sqrt(rho)
        stages = math.ceil(stages)
        condenser = reboiler = True
        bottoms = outs.get("bottoms", dist)
        T_top = T_top2 = dist["T"]
        T_bot = T_bot2 = bottoms["T"]
        murphree = False
    internal = stages - condenser - reboiler
    trays = internal if murphree else math.ceil(internal / e.c.tray_efficiency - 1e-9)
    if trays >= 1 and load > 0:
        area = load / F_FACTOR / NET_AREA  # m2, total cross section
        D = math.sqrt(4.0 * area / math.pi)
        H = TRAY_SPACING * trays + COLUMN_EXTRA_HEIGHT
        mat = e.material(uid)
        fp = vessel_pressure_factor(P_bot, D)
        e.item(cost, "shell", VESSEL, math.pi / 4.0 * D * D * H,
               (2.25, 1.82, FM_VESSEL[mat] * fp), f"D {D:.2f} m, H {H:.1f} m, {trays} trays")
        e.item(cost, "trays", TRAYS, area, FBM_TRAYS[mat] * trays * _quantity_factor(trays),
               f"{trays} sieve trays")
    if condenser:
        e.heat(cost, uid, "condenser", ures.get("condenser_duty"), (T_top2, T_top), P_top)
    if reboiler:
        e.heat(cost, uid, "reboiler", ures.get("reboiler_duty"), (T_bot2, T_bot), P_bot,
               exchanger=KETTLE)


def _quantity_factor(N: int) -> float:
    """Turton's tray quantity factor Fq: fewer than 20 trays cost more apiece."""
    if N >= 20:
        return 1.0
    lg = math.log10(N)
    return 10.0 ** (0.4771 + 0.08516 * lg - 0.3473 * lg * lg)


_HANDLERS = {
    "Heater": _heater,
    "HeatExchanger": _exchanger,
    "Flash": _drum,
    "MixingTank": _drum,
    "CSTR": _reactor,
    "PFR": _reactor,
    "StoichiometricReactor": _reactor,
    "YieldReactor": _reactor,
    "EquilibriumReactor": _reactor,
    "GibbsReactor": _reactor,
    "Pump": _pump,
    "Compressor": _compressor,
    "MultistageCompressor": _compressor,
    "DistillationColumn": _column,
    "Absorber": _column,
    "ShortcutColumn": _column,
}


def estimate(model: FlowsheetModel, result: SolveResult, packages) -> None:
    """Cost the solved flowsheet: fills ``result.costing`` and the cost results of each unit.
    ``packages(unit_id)`` gives the property package of a unit. Does nothing unless the
    model has a ``costing`` section."""
    if model.costing is None or not result.ok:
        return
    e = _Estimate(model, result, packages)
    c = model.costing
    units, not_costed = {}, []
    for um in model.units:
        ures = result.units.get(um.id)
        if ures is None or um.type in ("Feed", "Product"):
            continue
        own = c.units.get(um.id)
        if own is not None and own.exclude:
            continue
        try:
            cost = e.unit(um.id, um.type, ures)
        except SimulationError as exc:
            raise exc.located(um.id)
        if cost is None:
            if um.type not in ("Mixer", "Splitter", "Valve"):
                not_costed.append(um.id)
            continue
        purchased = sum(i["purchased_cost"] for i in cost.items)
        bare = sum(i["bare_module_cost"] for i in cost.items)
        operating = sum(x["cost"] for x in cost.energy)
        ures.update(purchased_cost=purchased, bare_module_cost=bare, operating_cost=operating)
        units[um.id] = {"type": um.type, "items": cost.items, "utilities": cost.energy,
                        "notes": cost.notes, "purchased_cost": purchased,
                        "bare_module_cost": bare, "operating_cost": operating}
    purchased = sum(u["purchased_cost"] for u in units.values())
    bare = sum(u["bare_module_cost"] for u in units.values())
    operating = sum(u["operating_cost"] for u in units.values())
    capital = TOTAL_MODULE * bare
    f = annualization_factor(c)
    result.costing = {
        "cepci": c.cepci, "hours_per_year": c.hours_per_year,
        "purchased_cost": purchased, "bare_module_cost": bare, "capital_cost": capital,
        "operating_cost": operating, "annualization_factor": f,
        "annualized_capital": f * capital, "tac": operating + f * capital,
        "units": units, "not_costed": not_costed,
    }
