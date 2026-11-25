import json
from pathlib import Path

import numpy as np
import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError
from scipy.optimize import brentq

from fluxsheet.cli import main
from fluxsheet.flowsheet import FlowsheetModel
from fluxsheet.property_sections import check_property_data
from fluxsheet.server import app
from fluxsheet.solver import solve
from fluxsheet.thermo_package import default_package, get_package, package_options

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "glycerol_local_property_data.json"
PAIR_NAMES = {"ethanol/glycerol", "water/glycerol"}


def test_preflight_and_run_use_the_same_data_in_cli_and_editor(tmp_path, capsys):
    model = FlowsheetModel.load(EXAMPLE)
    path = tmp_path / "round_trip.json"
    model.dump(path)
    assert FlowsheetModel.load(path) == model
    with TestClient(app) as client:
        payload = json.loads(path.read_text(encoding="utf-8"))
        api_check = client.post("/api/property-data/check", json=payload).json()
        api_run = client.post("/api/solve", json=payload).json()
        assert client.post("/api/methods/suggest", json=payload).json()["method"] == "nrtl"
    assert main(["check-property-data", str(path), "--json"]) == 0
    assert json.loads(capsys.readouterr().out) == api_check
    assert main(["run", str(path), "--json"]) == 0
    assert json.loads(capsys.readouterr().out) == api_run
    assert api_run["balance_status"] == "ok"
    assert {s["pair"] for s in api_run["parameter_sources"] if s["source"] == "Flowsheet"} == PAIR_NAMES


def test_glycerol_fallbacks_and_strict_unit_override():
    model = FlowsheetModel.load(EXAMPLE)
    data = model.property_data
    model.property_data = None
    report = check_property_data(model)
    assert report["ok"]
    assert {s["pair"] for s in report["sections"][0]["parameter_sources"]
            if s["source"] == "none (ideal)"} == PAIR_NAMES
    result = solve(model)
    assert result.ok and result.balance_status == "ok"
    model.property_method = "ideal"
    model.units[0].property_method = "nrtl"
    model.property_data = data.model_copy(update={"binary_parameters": [], "unifac_groups": []})
    report = check_property_data(model)
    assert not report["ok"] and report["errors"][0]["where"] == "FEED"
    result = solve(model)
    assert not result.ok and not result.streams and not result.converged
    assert result.error == report["errors"][0]
    model.property_data.unifac_groups = data.unifac_groups
    report = check_property_data(model)
    assert not report["ok"]
    assert {s["pair"] for s in report["sections"][1]["parameter_sources"]
            if s["source"] == "UNIFAC estimate"} == PAIR_NAMES


def test_group_assignment_preserves_database_properties_and_is_local():
    model = FlowsheetModel.load(EXAMPLE)
    model.property_method = "unifac"
    pkg, _ = default_package(model)
    db = get_package(tuple(model.components), "ideal")
    assert pkg.props is db.props
    assert pkg.consts.Tcs == db.consts.Tcs
    assert pkg.consts.UNIFAC_Rs == db.consts.UNIFAC_Rs
    assert pkg.consts.UNIFAC_Dortmund_groups[-1] == {2: 2, 3: 1, 14: 2, 81: 1}
    assert not db.consts.UNIFAC_Dortmund_groups[-1]
    assert pkg.psats(355.0, [0.4, 0.3, 0.3]) == pytest.approx(db.psats(355.0, [0.4, 0.3, 0.3]))
    result = solve(model)
    assert result.ok and result.balance_status == "ok"
    model.property_data = None
    assert "group assignment for glycerol" in check_property_data(model)["errors"][0]["message"]


def test_cubic_default_interactions_are_reported_in_preflight():
    model = FlowsheetModel.load(EXAMPLE)
    model.property_method = "srk"
    report = check_property_data(model)
    assert not report["ok"]
    assert "kij = 0" in report["sections"][0]["warnings"][0]
    assert "ethanol/glycerol" in report["errors"][0]["message"]


def test_local_priority_and_order_against_independent_nrtl_equations():
    model = FlowsheetModel.load(EXAMPLE)
    # Synthetic local parameters override the packaged ethanol/water fit and the ChemSep set.
    param = model.property_data.binary_parameters[0]
    model.property_data.binary_parameters.append(param.model_copy(update={
        "components": ("water", "ethanol"), "a": (0.1, -0.2), "b": (100.0, 200.0),
        "alpha": 0.3, "source": "Synthetic regression case",
    }))
    names = tuple(model.components)
    pkg = get_package(names, "nrtl", prefer_fitted=False, **package_options(model))
    options = package_options(model)
    reordered = get_package(names[::-1], "nrtl", **options)
    assert all(src["source"] == "Flowsheet" for src in pkg.parameter_sources
               if src.get("kind") != "unifac")
    x = np.array([0.4, 0.3, 0.3])
    a = np.zeros((3, 3))
    b = np.zeros((3, 3))
    alpha = np.zeros((3, 3))
    for p in model.property_data.binary_parameters:
        i, j = (names.index(c) for c in p.components)
        a[i, j], a[j, i] = p.a
        b[i, j], b[j, i] = p.b
        alpha[i, j] = alpha[j, i] = p.alpha

    def partial_pressures(T):
        tau = a + b / T
        G = np.exp(-alpha * tau)
        col_weights = x @ G
        col_avg_tau = (x @ (G * tau)) / col_weights
        lngamma = col_avg_tau + np.sum(G * (x / col_weights) * (tau - col_avg_tau), axis=1)
        return x * np.exp(lngamma) * [vp(T) for vp in pkg.props.VaporPressures]

    P = 100000.0
    T = brentq(lambda T: sum(partial_pressures(T)) - P, 330.0, 400.0)
    expected_y = partial_pressures(T) / P
    for package, reverse in ((pkg, False), (reordered, True)):
        sat = package.saturation_point((x[::-1] if reverse else x).tolist(), "bubble", P=P)
        # The checked flash returns the liquid-side state within 1 mK of the root.
        assert sat.T == pytest.approx(T, abs=1.1e-3)
        assert (sat.y[::-1] if reverse else sat.y) == pytest.approx(expected_y, abs=1e-6)


@pytest.mark.parametrize(("change", "message"), [
    ({"T_range_K": [400, 300]}, "validity range"),
    ({"convention": "aspen"}, "thermo"),
    ({"components": ["ethanol", "methanol"]}, "unknown or repeated"),
    ({"components": ["glycerol", "water"]}, "unknown or repeated"),
])
def test_bad_binary_declarations_are_rejected(change, message):
    raw = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    raw["property_data"]["binary_parameters"][0].update(change)
    with pytest.raises(ValidationError, match=message):
        FlowsheetModel.model_validate(raw)
