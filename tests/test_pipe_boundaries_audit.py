"""Grid refinement across two-phase heat transfer and elevation."""

import math

import pytest

from fluxsheet import solve


@pytest.mark.parametrize("vapor_fraction,ambient", [(0.1, 280), (0.9, 210)])
@pytest.mark.parametrize("elevation", [0, 0.5])
def test_two_phase_pipe_grid_and_energy(vapor_fraction, ambient, elevation):
    enthalpies = []
    for segments in (10, 40, 160):  # first order: each refinement quarters the change
        model = {
            "components": ["propane"], "property_method": "peng-robinson",
            "units": [
                {"id": "F", "type": "Feed", "params": {
                    "VF": vapor_fraction, "P": "2 bar", "flows": {"propane": 1}}},
                {"id": "P", "type": "Pipe", "params": {
                    "length": 3, "diameter": 0.1, "segments": segments,
                    "elevation": elevation, "U": 5, "T_ambient": ambient}},
            ],
            "streams": [{"id": "i", "from": "F.out", "to": "P.in"},
                        {"id": "o", "from": "P.out"}],
        }
        result = solve(model)
        assert result.ok, result.error
        inlet, outlet = result.streams["i"], result.streams["o"]
        assert math.isfinite(outlet["enthalpy_flow"])
        assert 0 <= outlet["vapor_fraction"] <= 1
        assert outlet["component_flows"]["propane"] >= 0
        assert outlet["enthalpy_flow"] == pytest.approx(
            inlet["enthalpy_flow"] + result.units["P"]["duty"]
            - inlet["mass_flow"] * 9.80665 * elevation, abs=1e-5)
        assert (outlet["vapor_fraction"] - vapor_fraction) * (ambient - inlet["T"]) > 0
        enthalpies.append(outlet["enthalpy_flow"])
    assert abs(enthalpies[2] - enthalpies[1]) < abs(enthalpies[1] - enthalpies[0])
