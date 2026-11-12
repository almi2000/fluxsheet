"""Round trips and phase transitions at deliberately chosen boundaries."""

import math

import pytest

from fluxsheet.thermo_package import PropertyPackage


@pytest.mark.parametrize("method", ["peng-robinson", "srk", "nrtl"])
@pytest.mark.parametrize("boundary", ["bubble", "dew"])
@pytest.mark.parametrize("offset", [-0.05, 0.05])
def test_near_phase_boundary_round_trip(method, boundary, offset):
    pkg = PropertyPackage(["benzene", "toluene"], method)
    reverse = PropertyPackage(["toluene", "benzene"], method)
    zs = [0.5, 0.5]
    temperature = pkg.saturation_point(zs, boundary, P=1e5).T + offset
    state = pkg.flash(zs, T=temperature, P=1e5)
    restored = pkg.flash(zs, P=1e5, H=state.H)
    swapped = reverse.flash(zs[::-1], T=temperature, P=1e5)
    assert math.isfinite(state.H)
    assert 0 <= state.VF <= 1
    assert restored.T == pytest.approx(temperature, abs=1e-3)
    assert restored.VF == pytest.approx(state.VF, abs=1e-5)
    assert swapped.H == pytest.approx(state.H, abs=1e-3)
    assert swapped.VF == pytest.approx(state.VF, abs=1e-5)
    if boundary == "bubble" and offset < 0:
        assert state.VF == pytest.approx(0, abs=1e-8)
    if boundary == "dew" and offset > 0:
        assert state.VF == pytest.approx(1, abs=1e-8)


@pytest.mark.parametrize("method", ["peng-robinson", "srk", "nrtl"])
@pytest.mark.parametrize("trace", [0.0, 1e-12, 1e-8])
def test_trace_component_round_trip_and_order(method, trace):
    pkg = PropertyPackage(["benzene", "toluene"], method)
    reverse = PropertyPackage(["toluene", "benzene"], method)
    zs = [trace, 1 - trace]
    state = pkg.flash(zs, T=370, P=1e5)
    restored = pkg.flash(zs, P=1e5, H=state.H)
    swapped = reverse.flash(zs[::-1], T=370, P=1e5)
    assert restored.T == pytest.approx(370, abs=1e-3)
    assert restored.VF == pytest.approx(state.VF, abs=1e-5)
    assert swapped.H == pytest.approx(state.H, abs=1e-3)
    assert swapped.VF == pytest.approx(state.VF, abs=1e-5)


@pytest.mark.parametrize("method", ["peng-robinson", "srk", "nrtl"])
@pytest.mark.parametrize("vapor_fraction", [0.0, 1e-5, 0.5, 1 - 1e-5, 1.0])
def test_pure_saturation_ph_round_trip(method, vapor_fraction):
    pkg = PropertyPackage(["propane"], method)
    state = pkg.flash([1.0], P=1e5, VF=vapor_fraction)
    restored = pkg.flash([1.0], P=1e5, H=state.H)
    assert math.isfinite(state.H)
    assert restored.T == pytest.approx(state.T, abs=1e-3)
    assert restored.VF == pytest.approx(vapor_fraction, abs=1e-5)
