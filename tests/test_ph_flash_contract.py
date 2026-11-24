"""An isobaric flash must realize the requested enthalpy through phase changes."""

from types import SimpleNamespace

import pytest

from fluxsheet.errors import SimulationError
from fluxsheet.thermo_package import PropertyPackage, get_package


@pytest.mark.parametrize("trace_flow", [0.0, 6.782944151e-7])
@pytest.mark.parametrize("nearby", [False, True])
def test_condensing_water_with_optional_sucrose_trace_hits_specification(trace_flow, nearby):
    # The apple-juice five-effect study cooling probe. Pure-water saturation sets a
    # separate limit for T and the energy lever rule puts VF close to 0.306.
    pressure = 67000.0
    water_flow = 52.44554437001688  # mol/s
    flow = water_flow + trace_flow
    composition = [water_flow / flow, trace_flow / flow]
    duty = -1_498_766.0  # W
    pkg = get_package(("water", "sucrose"), "ideal")
    inlet = pkg.flash(composition, T=362.25, P=pressure)
    target = inlet.H + duty / flow

    if nearby:
        outlet = pkg.flash_ph_near(composition, pressure, target, inlet.T)
    else:
        outlet = pkg.flash(composition, P=pressure, H=target)

    water = get_package(("water",), "ideal")
    liquid = water.flash([1.0], P=pressure, VF=0.0)
    vapor = water.flash([1.0], P=pressure, VF=1.0)
    expected_vf = (target - liquid.H) / (vapor.H - liquid.H)
    assert (outlet.H - inlet.H) * flow == pytest.approx(duty, abs=1.0)
    assert outlet.phase == "two-phase"
    assert outlet.T == pytest.approx(liquid.T, abs=2e-3)
    assert outlet.VF == pytest.approx(expected_vf, abs=1e-5)
    assert outlet.x is not None and outlet.y is not None
    if trace_flow:
        assert 0.0 <= outlet.y[1] < composition[1] < outlet.x[1]
    combined = (1.0 - outlet.VF) * outlet.x[1] + outlet.VF * outlet.y[1]
    assert combined == pytest.approx(composition[1], abs=1e-12)


@pytest.mark.parametrize("nearby", [False, True])
def test_flash_reports_enthalpy_residual_when_all_candidates_miss_target(monkeypatch, nearby):
    pkg = PropertyPackage(["water"], "ideal")
    composition = [1.0]
    pressure = 101325.0
    raw = pkg._flasher.flash(zs=composition, P=pressure, T=300.0)
    actual = float(raw.H())
    target = actual + 1234.0
    # Model a backend or fallback that claims success but returns an unchanged
    # thermodynamic state. A plausible T/P is insufficient to satisfy P/H.
    monkeypatch.setattr(pkg, "_flasher", SimpleNamespace(flash=lambda **_: raw))
    monkeypatch.setattr(pkg, "_by_temperature", lambda *_: raw)
    monkeypatch.setattr(pkg, "_ph_by_vapor_fraction", lambda *_: raw)

    with pytest.raises(SimulationError) as failure:
        if nearby:
            pkg.flash_ph_near(composition, pressure, target, 300.0)
        else:
            pkg.flash(composition, P=pressure, H=target)

    assert f"requested H={target:.9g}" in failure.value.message
    assert f"calculated H={actual:.9g}" in failure.value.message
    assert "residual=-1234 J/mol" in failure.value.message
