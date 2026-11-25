"""Verify the solid reference case against a sensible/latent heat calculation."""

import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from fluxsheet.cli import main
from fluxsheet.dynamics.simulator import simulate
from fluxsheet.errors import SimulationError
from fluxsheet.server import app
from fluxsheet.solver import solve
from fluxsheet.thermo_package import get_package

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "water_ice_freezing.json"

# Rounded data near 0 °C. See the IAPWS releases and constants source in architecture.md.
MELTING_T = 273.15  # K
FUSION = 6010.0  # J/mol
ICE_CP = 37.8  # J/(mol K)
WATER_CP = 76.0  # J/(mol K)
MW = 0.01801528  # kg/mol


def _model():
    return json.loads(EXAMPLE.read_text(encoding="utf-8"))


@pytest.mark.parametrize("pressure", [90_000.0, 101_325.0, 110_000.0])
def test_t_p_and_p_h_have_a_pure_ice_branch_and_a_latent_plateau(pressure):
    pkg = get_package(("water",), "water-ice")
    liquid = pkg.flash([1.0], T=MELTING_T, P=pressure)
    assert liquid.phase == "liquid" and liquid.x == (1.0,)
    assert liquid.liquid_fraction == 1.0

    at_melting_ice = pkg.flash([1.0], P=pressure, H=liquid.H - FUSION)
    assert at_melting_ice.phase == "solid"
    assert at_melting_ice.T == MELTING_T
    assert at_melting_ice.S == pytest.approx(liquid.S - FUSION / MELTING_T)
    assert at_melting_ice.V == pytest.approx(MW / 917.0)

    for fraction in (0.05, 0.5, 0.95):
        at_melting = pkg.flash([1.0], P=pressure, H=liquid.H - fraction * FUSION)
        assert at_melting.T == MELTING_T and at_melting.phase == "liquid-solid"
        assert at_melting.VF == 0.0
        assert at_melting.solid_fraction == pytest.approx(fraction)
        assert at_melting.liquid_fraction == pytest.approx(1.0 - fraction)

    for T, delta in (
        (263.15, -FUSION + ICE_CP * (263.15 - MELTING_T)),
        (293.15, WATER_CP * (293.15 - MELTING_T)),
    ):
        at_temperature = pkg.flash([1.0], T=T, P=pressure)
        assert at_temperature.H - liquid.H == pytest.approx(delta)
        at_enthalpy = pkg.flash([1.0], H=at_temperature.H, P=pressure)
        assert at_enthalpy.T == pytest.approx(T)
        assert at_enthalpy.solid_fraction == at_temperature.solid_fraction


@pytest.mark.parametrize(("duty", "fraction", "temperature"), [
    (-WATER_CP * 5 - FUSION / 2, 0.5, MELTING_T),
    (-WATER_CP * 5 - FUSION - ICE_CP * 5, 1.0, MELTING_T - 5),
])
def test_steady_flowsheet_freezes_water_and_closes_balances(duty, fraction, temperature):
    model = _model()
    model["units"][1]["params"]["duty"] = duty  # W for a water flow of 1 mol/s
    result = solve(model)
    assert result.ok, result.error
    feed = result.streams["WATER"]
    product = result.streams["ICE_WATER"]
    assert product["T"] == pytest.approx(temperature)
    assert product["mass_flow"] == pytest.approx(MW)
    assert product["component_flows"] == {"water": 1.0}
    assert product["solid_fraction"] == pytest.approx(fraction)
    assert product["solid_composition"] == {"water": 1.0}
    assert product["liquid_fraction"] == pytest.approx(1.0 - fraction)
    assert product["liquid_composition"] == ({"water": 1.0} if fraction < 1 else None)
    fractions = [product[f"{phase}_fraction"] for phase in ("liquid", "solid", "vapor")]
    assert sum(fractions) == pytest.approx(1.0)
    assert product["enthalpy_flow"] - feed["enthalpy_flow"] == pytest.approx(duty)
    assert result.units["COOLER"]["duty"] == pytest.approx(duty)
    assert result.specification_status == "ok"
    assert result.balance_status == "ok"
    checks = {(c["scope"], c["quantity"]): c for c in result.balance_checks}
    for scope in ("COOLER", "flowsheet"):
        for quantity in ("mass", "component", "energy"):
            assert checks[scope, quantity]["residual"] == pytest.approx(0.0)


def test_steady_flowsheet_can_melt_ice_at_a_specified_temperature():
    model = _model()
    model["units"][0]["params"]["T"] = "-5 degC"
    model["units"][1]["params"] = {"T": "5 degC"}
    result = solve(model)
    assert result.ok, result.error
    assert result.streams["WATER"]["solid_fraction"] == 1.0
    assert result.streams["ICE_WATER"]["phase"] == "liquid"
    assert result.units["COOLER"]["duty"] == pytest.approx(ICE_CP * 5 + FUSION + WATER_CP * 5)
    assert result.balance_status == "ok"


def test_water_ice_results_reach_the_api_and_cli(capsys):
    response = TestClient(app).post("/api/solve", json=_model())
    assert response.status_code == 200
    result = response.json()
    assert result["ok"], result["error"]
    assert result["streams"]["ICE_WATER"]["solid_fraction"] == pytest.approx(0.5)

    assert main(["run", str(EXAMPLE), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["streams"] == result["streams"]
    assert main(["run", str(EXAMPLE)]) == 0
    output = capsys.readouterr().out
    assert "Liq. frac." in output and "Solid frac." in output
    product_row = next(line for line in output.splitlines() if line.startswith("ICE_WATER"))
    assert product_row.split()[-2:] == ["0.5000", "0.5000"]


@pytest.mark.parametrize(("spec", "message"), [
    ({"P": 120_000.0, "T": MELTING_T}, "pressure"),
    ({"P": 100_000.0, "T": 320.0}, "temperature"),
    ({"P": 100_000.0, "H": -100.0}, "enthalpy"),
    ({"P": 100_000.0, "VF": 0.0}, "supports only T-P and P-H"),
])
def test_water_ice_rejects_unsupported_flash_specs(spec, message):
    pkg = get_package(("water",), "water-ice")
    with pytest.raises(SimulationError, match=message):
        pkg.flash([1.0], **spec)


def test_water_ice_does_not_silently_switch_to_a_fluid_only_unit_or_method():
    model = _model()
    model["units"][1]["property_method"] = "ideal"
    crossed = solve(model)
    assert not crossed.ok
    assert crossed.error["where"] == "WATER"
    assert "cannot re-flash across" in crossed.error["message"]

    model = _model()
    model["units"][1]["type"] = "Valve"
    model["units"][1]["params"] = {"P": "1 bar"}
    unsupported = solve(model)
    assert not unsupported.ok
    assert unsupported.error["where"] == "COOLER"
    assert "does not support Valve" in unsupported.error["message"]


def test_water_ice_requires_steady_state():
    model = _model()
    model["dynamics"] = {"t_end": "60 s", "units": {"COOLER": {}}}
    result = simulate(model)
    assert not result["ok"]
    assert result["error"]["where"] == "dynamics"
    assert "supports only steady-state" in result["error"]["message"]
