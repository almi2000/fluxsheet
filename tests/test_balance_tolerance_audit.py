"""Independent arithmetic checks for balance tolerances."""

import pytest

from fluxsheet.balances import _record


@pytest.mark.parametrize("flow", [0.0, 1e-12, 1e-8, 1e-4])
def test_component_tolerance_at_zero_and_trace_flow(flow):
    checks = []
    _record(checks, "unit", "component", "trace", [flow], [flow], "mol/s", 1e-10)
    assert checks[0]["status"] == "ok"
    assert checks[0]["tolerance"] == pytest.approx(1e-10 + 1e-5 * flow)


@pytest.mark.parametrize("residual,expected", [(0.999e-10, "ok"), (1.001e-10, "warning")])
def test_absolute_tolerance_boundary(residual, expected):
    checks = []
    _record(checks, "unit", "component", "trace", [residual], [0.0], "mol/s", 1e-10)
    assert checks[0]["status"] == expected


@pytest.mark.parametrize("difference,expected", [(19.9, "ok"), (20.2, "warning")])
def test_counterflow_energy_uses_absolute_terms(difference, expected):
    checks = []
    _record(checks, "unit", "energy", None, [1e6, -1e6], [difference], "W", 0.1)
    assert checks[0]["tolerance"] == pytest.approx(20.1)
    assert checks[0]["status"] == expected
