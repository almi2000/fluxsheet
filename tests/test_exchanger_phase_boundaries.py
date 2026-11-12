"""Compare phase-changing exchanger duties with independent PH quadrature."""

from itertools import pairwise

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

from fluxsheet import solve
from fluxsheet.thermo_package import get_package


@pytest.mark.parametrize("arrangement", ["countercurrent", "cocurrent"])
@pytest.mark.parametrize("changing_side", ["hot", "cold"])
@pytest.mark.parametrize("vapor_fraction", [0.001, 0.5, 0.999])
def test_phase_changing_exchanger_independent_ph_integral(arrangement, changing_side,
                                                       vapor_fraction):
    pkg = get_package(("propane",), "peng-robinson")
    hot = pkg.flash([1.0], T=320, P=3e5)
    cold = pkg.flash([1.0], T=200, P=1e5)
    hot_flow, cold_flow = (1.0, 10.0) if changing_side == "hot" else (10.0, 1.0)
    saturation = pkg.flash([1.0], P=3e5 if changing_side == "hot" else 1e5,
                           VF=vapor_fraction)
    duty = (hot_flow * (hot.H - saturation.H) if changing_side == "hot"
            else cold_flow * (saturation.H - cold.H))
    model = {
        "components": ["propane"], "property_method": "peng-robinson",
        "units": [
            {"id": "H", "type": "Feed", "params": {
                "T": 320, "P": 3e5, "flows": {"propane": hot_flow}}},
            {"id": "C", "type": "Feed", "params": {
                "T": 200, "P": 1e5, "flows": {"propane": cold_flow}}},
            {"id": "X", "type": "HeatExchanger", "params": {
                "duty": duty, "arrangement": arrangement}},
        ],
        "streams": [
            {"id": "hi", "from": "H.out", "to": "X.hot_in"},
            {"id": "ci", "from": "C.out", "to": "X.cold_in"},
            {"id": "ho", "from": "X.hot_out"},
            {"id": "co", "from": "X.cold_out"},
        ],
    }
    result = solve(model)
    assert result.ok, result.error
    hi, ci, ho, co = (result.streams[s] for s in ("hi", "ci", "ho", "co"))
    assert hi["enthalpy_flow"] - ho["enthalpy_flow"] == pytest.approx(duty, abs=1e-5)
    assert co["enthalpy_flow"] - ci["enthalpy_flow"] == pytest.approx(duty, abs=1e-5)
    assert (ho if changing_side == "hot" else co)["vapor_fraction"] == pytest.approx(
        vapor_fraction, abs=1e-5)

    def difference(q):
        hot_temperature = pkg.flash([1.0], P=3e5, H=hot.H - q / hot_flow).T
        cold_q = duty - q if arrangement == "countercurrent" else q
        cold_temperature = pkg.flash([1.0], P=1e5, H=cold.H + cold_q / cold_flow).T
        return hot_temperature - cold_temperature

    points = []
    for vf in (0.0, 1.0):
        hot_q = hot_flow * (hot.H - pkg.flash([1.0], P=3e5, VF=vf).H)
        cold_q = cold_flow * (pkg.flash([1.0], P=1e5, VF=vf).H - cold.H)
        points.extend([hot_q, duty - cold_q if arrangement == "countercurrent" else cold_q])
    breaks = sorted(p for p in points if 0 < p < duty)
    intervals = [0.0, *breaks, duty]
    def integrate(epsabs, epsrel):
        return sum(quad(lambda q: 1 / difference(q), a, b, epsabs=epsabs, epsrel=epsrel)[0]
                   for a, b in pairwise(intervals))

    ua = integrate(1e-6, 1e-8)
    assert ua == pytest.approx(integrate(1e-8, 1e-10), rel=1e-6)
    approach = min(difference(q) for q in intervals)
    for a, b in pairwise(intervals):
        minimum = minimize_scalar(difference, bounds=(a, b), method="bounded")
        approach = min(approach, minimum.fun)
    assert result.units["X"]["UA"] == pytest.approx(ua, rel=0.005)
    assert result.units["X"]["min_approach"] == pytest.approx(approach, abs=0.05)
    assert np.isfinite(ua)
