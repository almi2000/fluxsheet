"""A cached single-phase state must not cross into a metastable phase."""

import pytest

from fluxsheet.contents_flash import ContentsFlash
from fluxsheet.thermo_package import get_package


@pytest.mark.parametrize("method", ["peng-robinson", "srk", "nrtl"])
def test_single_phase_shortcut_does_not_cross_boiling_point(method):
    pkg = get_package(("propane",), method)
    saturation = pkg.flash([1.0], P=1e5, VF=0)
    start = pkg.flash([1.0], P=1e5, T=saturation.T - 0.5)
    metastable_h = pkg.phase_props("liquid", saturation.T + 0.5, 1e5, [1.0]).H
    expected = pkg.flash([1.0], P=1e5, H=metastable_h)
    flash = ContentsFlash(pkg, start.T)
    flash([1.0], 1e5, start.H)
    state = flash([1.0], 1e5, metastable_h)
    assert 0 < expected.VF < 1
    assert state.T == pytest.approx(expected.T, abs=1e-3)
    assert state.VF == pytest.approx(expected.VF, abs=1e-5)


def test_temperature_shortcut_does_not_cross_boiling_point():
    pkg = get_package(("propane",), "peng-robinson")
    saturation = pkg.flash([1.0], P=1e5, VF=0)
    flash = ContentsFlash(pkg, saturation.T - 0.5)
    flash.at_temperature([1.0], 1e5, saturation.T - 0.5)
    assert flash.at_temperature([1.0], 1e5, saturation.T + 0.5).phase == "vapor"


@pytest.mark.parametrize("method", ["peng-robinson", "srk", "nrtl"])
def test_single_phase_shortcut_does_not_cross_dew_point(method):
    pkg = get_package(("propane",), method)
    saturation = pkg.flash([1.0], P=1e5, VF=1)
    start = pkg.flash([1.0], P=1e5, T=saturation.T + 0.5)
    metastable_h = pkg.phase_props("vapor", saturation.T - 0.5, 1e5, [1.0]).H
    expected = pkg.flash([1.0], P=1e5, H=metastable_h)
    flash = ContentsFlash(pkg, start.T)
    flash([1.0], 1e5, start.H)
    state = flash([1.0], 1e5, metastable_h)
    assert 0 < expected.VF < 1
    assert state.T == pytest.approx(expected.T, abs=1e-3)
    assert state.VF == pytest.approx(expected.VF, abs=1e-5)


def test_single_phase_shortcut_retains_single_phase_round_trip():
    pkg = get_package(("propane",), "peng-robinson")
    start = pkg.flash([1.0], P=1e5, T=200)
    target = pkg.flash([1.0], P=1e5, T=201)
    flash = ContentsFlash(pkg, start.T)
    flash([1.0], 1e5, start.H)
    state = flash([1.0], 1e5, target.H)
    assert state.phase == "liquid"
    assert state.T == pytest.approx(target.T, abs=1e-6)
