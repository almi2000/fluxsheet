"""Transient vessel inventories through boiling and subsequent cooling."""

import numpy as np
import pytest
from scipy.integrate import trapezoid

from fluxsheet.dynamics import simulate
from fluxsheet.thermo_package import get_package


def test_vessel_boil_and_cool_inventory_refines():
    pkg = get_package(("propane",), "peng-robinson")
    saturation = pkg.flash([1.0], P=1e5, VF=0).T
    inlet_h = pkg.flash([1.0], T=220, P=1e5).H
    errors = []
    for points, tolerance in ((121, 1e-5), (241, 1e-6), (481, 1e-7)):
        model = {
            "components": ["propane"], "property_method": "peng-robinson",
            "units": [
                {"id": "F", "type": "Feed", "params": {
                    "T": 220, "P": "1 bar", "flows": {"propane": 1}}},
                {"id": "T", "type": "MixingTank", "params": {"P": "1 bar"}},
            ],
            "streams": [{"id": "in", "from": "F.out", "to": "T.in"},
                        {"id": "liq", "from": "T.liquid"},
                        {"id": "vap", "from": "T.vapor"}],
            "dynamics": {
                "t_end": 60, "points": points, "rtol": tolerance,
                "units": {"T": {"diameter": 0.1, "height": 0.2,
                                "initial_level": 0.5, "vapor_time": 1}},
                "events": [{"time": 5, "unit": "T", "param": "duty", "value": 4000},
                           {"time": 15, "unit": "T", "param": "duty", "value": 0}],
            },
        }
        result = simulate(model)
        assert result["ok"], result.get("error")
        data = {s["key"]: np.asarray(s["values"]) for s in result["series"]}
        time = np.asarray(result["time"])
        temperature = data["unit.T.T"]
        liquid = data["unit.T.holdup"]
        vapor = data["unit.T.vapor_flow"]  # vapor residence time is 1 s
        flow_liquid = data["stream.liq.molar_flow"]
        flow_vapor = data["stream.vap.molar_flow"]
        h_liquid = np.array([pkg.phase_props("liquid", T, 1e5, [1]).H for T in temperature])
        h_vapor = np.array([pkg.phase_props("vapor", T, 1e5, [1]).H for T in temperature])
        total = liquid + vapor
        enthalpy = liquid * h_liquid + vapor * h_vapor
        mass_error = ((total[-1] - total[0])
                      - trapezoid(1 - flow_liquid - flow_vapor, time))
        # Integrate the rectangular duty pulse exactly to avoid sampling its jump endpoints.
        energy_error = ((enthalpy[-1] - enthalpy[0])
                        - trapezoid(inlet_h - flow_liquid * h_liquid - flow_vapor * h_vapor,
                                    time) - 4000 * 10)
        errors.append((abs(mass_error), abs(energy_error)))
        assert max(flow_vapor) > 0.1
        assert max(temperature) <= saturation + 1e-3
        assert temperature[-1] == pytest.approx(220, abs=0.2)
        assert flow_vapor[-1] == pytest.approx(0, abs=1e-6)
    for index in (0, 1):
        assert errors[2][index] < errors[1][index] < errors[0][index]
    assert errors[-1][0] < 0.001  # mol
    assert errors[-1][1] < 5  # J
