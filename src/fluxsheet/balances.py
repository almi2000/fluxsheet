"""Independent material and energy checks on the converged flowsheet."""

from __future__ import annotations

import math

_REACTORS = {"StoichiometricReactor", "EquilibriumReactor", "GibbsReactor", "CSTR", "PFR",
             "BatchReactor"}
_COLUMNS = {"DistillationColumn", "ShortcutColumn", "Absorber"}
_NO_DUTY = {"Mixer", "Splitter", "Valve", "HeatExchanger"}
_DUTY = {"Heater", "Flash", "ComponentSplitter", "MixingTank", "Decanter",
         "ThinFilmEvaporator", "ExtractionColumn", "Pipe", "YieldReactor"} | _REACTORS
_MACHINES = {"Pump", "Compressor", "Turbine", "MultistageCompressor"}


def _reactive(unit):
    return unit.type in _REACTORS or bool(unit.p.get("reactions"))


def _record(out, scope, quantity, component, inputs, outputs, unit, absolute):
    terms = [*inputs, *outputs]
    if any(v is None or not math.isfinite(v) for v in terms):
        out.append({"scope": scope, "quantity": quantity, "component": component, "residual": None,
                        "unit": unit, "tolerance": None, "status": "not_checkable"})
        return
    residual = sum(inputs) - sum(outputs)
    tolerance = absolute + 1e-5 * max(sum(abs(v) for v in inputs),
                                       sum(abs(v) for v in outputs))
    out.append({"scope": scope, "quantity": quantity, "component": component,
                    "residual": float(residual), "unit": unit, "tolerance": float(tolerance),
                    "status": "ok" if abs(residual) <= tolerance else "warning"})


def _not_checkable(out, scope, quantity, unit):
    out.append({"scope": scope, "quantity": quantity, "component": None, "residual": None,
                    "unit": unit, "tolerance": None, "status": "not_checkable"})


def _mass(stream, mws):
    return sum(f * mw for f, mw in zip(stream.flows, mws)) / 1000.0


def _energy(stream, formation):
    if formation is None:
        return stream.enthalpy_flow
    if any(h is None and f > 0.0 for h, f in zip(formation, stream.flows)):
        return None
    return stream.enthalpy_flow + sum(f * (h or 0.0) for f, h in zip(stream.flows, formation))


def _utilities(unit, inlets, mws):
    r, kind = unit.results, unit.type
    if kind in _NO_DUTY:
        return [0.0]
    if kind in _COLUMNS:
        duties = [r.get("condenser_duty", 0.0), r.get("reboiler_duty", 0.0),
                  r.get("pumparound_duty", 0.0)]
        duties.extend(d["duty"] for d in (unit.p.get("stage_duties") or []))
        return duties
    if kind in _MACHINES:
        return [r.get("power", 0.0), r.get("cooler_duty", 0.0)]
    if kind in _DUTY:
        terms = [r.get("duty", 0.0)]
        if kind == "BatchReactor" and unit.p["vessel"] == "constant_volume":
            if not inlets or sum(s.total_flow for s in inlets) <= 0.0:
                return terms
            initial_pressure = unit.p.get("P") or inlets[0].P
            # The batch model initializes the charge at vessel pressure and feed T.
            # A different feed pressure also requires a filling-work balance.
            if any(not math.isclose(s.P, initial_pressure, rel_tol=1e-8, abs_tol=1e-3)
                   for s in inlets):
                return None
            if not r.get("cycle_time") or r.get("P_final") is None or r.get("max_volume") is None:
                return None
            terms.append((r["P_final"] - initial_pressure) * r["max_volume"] / r["cycle_time"])
        if kind == "Pipe":
            terms.append(-sum(_mass(s, mws) for s in inlets) * 9.80665
                         * (unit.p.get("elevation") or 0.0))
        return terms
    return None


def _materials(out, scope, inputs, outputs, pkg, reactive, yield_reactor=False):
    _record(out, scope, "mass", None, [_mass(s, pkg.mws) for s in inputs],
            [_mass(s, pkg.mws) for s in outputs], "kg/s", 1e-10)
    if yield_reactor:
        _not_checkable(out, scope, "elements", "mol/s")
        return
    if not reactive:
        for i, name in enumerate(pkg.names):
            _record(out, scope, "component", name, [s.flows[i] for s in inputs],
                    [s.flows[i] for s in outputs], "mol/s", 1e-10)
        return
    atoms = pkg.consts.atomss
    involved = [i for i in range(pkg.n)
                if any(s.flows[i] > 0 for s in inputs + outputs)]
    if any(not atoms[i] for i in involved):
        _not_checkable(out, scope, "elements", "mol/s")
        return
    elements = sorted({e for i in involved for e in atoms[i]})
    for element in elements:
        coeff = [a.get(element, 0) if a else 0 for a in atoms]
        _record(out, scope, "element", element,
                [sum(f * n for f, n in zip(s.flows, coeff)) for s in inputs],
                [sum(f * n for f, n in zip(s.flows, coeff)) for s in outputs],
                "mol/s", 1e-10)


def check_balances(solver, pkg, streams, crossings):
    checks = []
    active = [(uid, u) for uid, u in solver.units.items()
              if u.type not in ("Feed", "Product") and u.conserves_mass
              and uid in solver.balance_io]
    energy_terms = []
    reactive = any(_reactive(u) for _, u in active)
    has_yield = any(u.type == "YieldReactor" for _, u in active)
    for uid, unit in active:
        inputs, outputs = solver.balance_io[uid]
        reacts = _reactive(unit)
        _materials(checks, uid, inputs, outputs, pkg, reacts, unit.type == "YieldReactor")
        utilities = _utilities(unit, inputs, pkg.mws)
        energy_terms.append(utilities)
        if utilities is None:
            _not_checkable(checks, uid, "energy", "W")
            continue
        formation = pkg.consts.Hfgs if reacts or unit.type == "YieldReactor" else None
        _record(checks, uid, "energy", None,
                [_energy(s, formation) for s in inputs] + utilities,
                [_energy(s, formation) for s in outputs], "W", 0.1)

    # External feeds and products, plus boundaries at non-conserving logic blocks.
    by_id = {c.id: c for c in solver.connections}
    inputs, outputs = [], []
    for sid, stream in streams.items():
        conn = by_id.get(sid)
        src = conn.src_unit if conn else sid.rsplit(".", 1)[0]
        src_unit = solver.units.get(src)
        if src_unit and (src_unit.type == "Feed" or not src_unit.conserves_mass):
            inputs.append(stream)
        if conn is None or conn.dst_unit is None:
            outputs.append(stream)
        else:
            dst = solver.units[conn.dst_unit]
            if dst.type == "Product" or not dst.conserves_mass:
                outputs.append(stream)
    _materials(checks, "flowsheet", inputs, outputs, pkg, reactive, has_yield)
    if any(t is None for t in energy_terms):
        _not_checkable(checks, "flowsheet", "energy", "W")
    else:
        formation = pkg.consts.Hfgs if reactive or has_yield else None
        duties = [x for group in energy_terms for x in group]
        # A crossing into a virtual boundary or product occurs outside the
        # physical flowsheet being balanced. Its incoming stream is already
        # counted at the producer's state. Crossings into physical units do
        # contribute, including those from a virtual source.
        crossings_energy = [c["enthalpy_change"] for c in crossings
                            if c["to_unit"] in solver.units
                            and solver.units[c["to_unit"]].type != "Product"
                            and solver.units[c["to_unit"]].conserves_mass]
        _record(checks, "flowsheet", "energy", None,
                [_energy(s, formation) for s in inputs] + duties + crossings_energy,
                [_energy(s, formation) for s in outputs], "W", 0.1)
    return checks
