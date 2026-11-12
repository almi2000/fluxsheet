"""Synthetic data exercise comparison logic, not experimental validation."""

import pytest
from scripts.equipment_metrics import column, dynamic_step, errors, interpolate, kinetics


def test_interpolation_and_bounds():
    assert interpolate([0, 2, 4], [0, 4, 8], [1, 4]) == [2, 8]
    with pytest.raises(ValueError, match="extrapolate"):
        interpolate([0, 1], [2, 3], [1.1])
    with pytest.raises(ValueError, match="strictly increasing"):
        interpolate([0, 0], [2, 3], [0])


def test_fixed_kinetics_and_column_limits():
    assert errors([0, .5], [.02, .53], .03, .08)["passed"]
    assert not errors([0, .5], [.02, .59], .03, .08)["passed"]
    result = errors([300, 320], [301, 323], 2, 5)
    assert result["maximum_error_index"] == 1
    assert result["passed"]
    with pytest.raises(ValueError, match="near zero"):
        errors([0], [1], .1, .2, relative=True, zero_floor=1e-6)
    with pytest.raises(ValueError, match="finite"):
        errors([float("nan")], [1], .1, .2)


def test_dynamic_detrending_and_noise_gate():
    result = dynamic_step([0, 1, 2], [1, 1, 1], [1, 1.5, 2], [1, 1.5, 2], 2)
    assert result["passed"]
    assert result["step_amplitude"] == 1
    assert not dynamic_step([0, 1, 2], [1, 1, 1], [1, 1.5, 2], [1, 1.8, 2.2], 2)["passed"]
    with pytest.raises(ValueError, match="five times"):
        dynamic_step([0, 1, 2], [0, 1, 0], [0, .2], [0, .2], 1)


def test_registered_limits_and_duty_zero():
    assert kinetics([.2, .5], [.21, .52])["passed"]
    assert not kinetics([.2, .5], [.3, .6])["passed"]
    assert column([.2, .5], [.21, .52], "composition")["passed"]
    assert column([100, 200], [105, 210], "duty")["passed"]
    with pytest.raises(ValueError, match="near zero"):
        column([0], [1], "duty", zero_floor=1e-6)
    with pytest.raises(ValueError, match="unsupported"):
        column([1], [1], "pressure")
