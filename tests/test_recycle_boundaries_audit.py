"""Recycle convergence at increasing recirculation ratios."""

import pytest

from fluxsheet import solve


@pytest.mark.parametrize("acceleration", ["none", "wegstein", "broyden"])
@pytest.mark.parametrize("fraction", [0, 0.5, 0.9, 0.99])
def test_two_phase_recycle_grid(acceleration, fraction):
    model = {
        "components": ["propane"], "property_method": "peng-robinson",
        "solver": {"acceleration": acceleration, "tolerance": 1e-8,
                   "max_iterations": 3000},
        "units": [
            {"id": "F", "type": "Feed", "params": {
                "VF": 0.5, "P": "1 bar", "flows": {"propane": 1}}},
            {"id": "M", "type": "Mixer"},
            {"id": "S", "type": "Splitter", "params": {"split": fraction}},
        ],
        "streams": [
            {"id": "f", "from": "F.out", "to": "M.in"},
            {"id": "m", "from": "M.out", "to": "S.in"},
            {"id": "r", "from": "S.out1", "to": "M.in"},
            {"id": "p", "from": "S.out2"},
        ],
    }
    result = solve(model)
    assert result.ok, result.error
    assert result.converged
    feed, product, loop = (result.streams[name] for name in ("f", "p", "m"))
    assert loop["molar_flow"] == pytest.approx(1 / (1 - fraction), rel=2e-6)
    assert product["molar_flow"] == pytest.approx(feed["molar_flow"], rel=2e-6)
    assert product["vapor_fraction"] == pytest.approx(0.5, abs=2e-6)
    assert product["enthalpy_flow"] == pytest.approx(feed["enthalpy_flow"], rel=2e-6, abs=1e-5)
