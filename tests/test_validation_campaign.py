"""Campaign status and structured summaries on synthetic input, and the feature reference
campaign (docs/validation-campaign-features.json) run against its recorded results."""

import hashlib
import json
import math
import sys
from itertools import pairwise

import pytest
from scripts import validate_features, validate_models


def test_parse_legacy_summary_and_phase_errors():
    output = ("T=298.15 K aqueous: predicted=0.1 measured=0.2 absolute error=0.1\n"
              "T=298.15 K organic: predicted=0.9 measured=0.8 absolute error=0.1\n"
              "AAD=0.100000, max=0.100000 at 298.15 K (organic); limits 0.05/0.1: FAIL\n")
    metrics, phases = validate_models.parse_legacy("butanol_2010", output, 1, 2)
    assert metrics[0]["mean_absolute_error"] == 0.1
    assert metrics[0]["unit"] == "fraction"
    assert phases["aqueous"]["maximum_absolute_error"] == 0.1
    assert phases["organic"]["points"] == 1
    assert metrics[0]["quantity"] == "butanol mass fraction"
    assert metrics[0]["maximum_error_point"] == "298.15 K (organic)"
    with pytest.raises(ValueError, match="unexpected validator output"):
        validate_models.parse_legacy("butanol_2010", "bad output", 1, 2)


def test_parse_legacy_names_each_isotherm_without_the_verdict():
    output = ("T=333.15 K n=11 AAD=1.2000% max=2.5000% at x_benzene=0.5 "
              "(predicted=1.0 Pa, measured=2.0 Pa) PASS\n"
              "T=353.15 K n=11 AAD=6.0000% max=12.0000% at x_benzene=0.1 "
              "(predicted=1.0 Pa, measured=2.0 Pa) FAIL\n")
    metrics, phases = validate_models.parse_legacy("btx", output, 2, 11)
    assert [m["quantity"] for m in metrics] == ["bubble pressure, T=333.15 K",
                                                "bubble pressure, T=353.15 K"]
    assert [m["unit"] for m in metrics] == ["%", "%"]
    assert metrics[1]["maximum_absolute_error"] == 12.0
    assert metrics[0]["maximum_error_point"] == (
        "x_benzene=0.5 (predicted=1.0 Pa, measured=2.0 Pa)")
    assert not phases


def test_registered_thresholds_are_the_validators_limits():
    manifest = json.loads((validate_models.ROOT / "docs" / "validation-campaign.json")
                          .read_text(encoding="utf-8"))
    assert {s["id"]: s["thresholds"] for s in manifest["sources"]} == validate_models.THRESHOLDS
    assert set(validate_models.QUANTITIES) | {"butanol_2007"} == set(validate_models.SCRIPTS)
    for name, quantities in validate_models.QUANTITIES.items():
        assert len(quantities) == len(validate_models.THRESHOLDS[name]) == (
            validate_models.SCRIPTS[name][2])


@pytest.mark.parametrize("missing", ["doi", "sha256"])
def test_incomplete_source_is_an_invalid_campaign(tmp_path, monkeypatch, missing):
    source = {"id": "btx", "doi": "10.1234/example", "sha256": "0" * 64, "methods": ["ideal"],
              "thresholds": [[5, 10], [5, 10]]}
    del source[missing]
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"version": 1, "sources": [source]}))
    monkeypatch.setattr(sys, "argv", ["validate_models", "--manifest", str(manifest),
                                    "--cache-dir", str(tmp_path / "cache"),
                                    "--output-dir", str(tmp_path / "output")])
    with pytest.raises(SystemExit) as exc:
        validate_models.main()
    assert exc.value.code == 2


@pytest.mark.parametrize("mode,expected", [("pass", 0), ("download", 0), ("fail", 1),
                                            ("bad_hash", 2), ("missing", 2)])
def test_campaign_exit_status_and_reports(tmp_path, monkeypatch, mode, expected):
    raw = tmp_path / "raw.xml"
    raw.write_bytes(b"synthetic input")
    sha = hashlib.sha256(raw.read_bytes()).hexdigest()
    if mode == "bad_hash":
        sha = "0" * 64
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"version": 1, "sources": [
        {"id": "btx", "doi": "10.1234/example", "sha256": sha, "methods": ["ideal"],
         "thresholds": [[5, 10], [5, 10]]}]}))
    output = tmp_path / "output"
    flag = "--download" if mode == "download" else "--offline"
    monkeypatch.setattr(sys, "argv", ["validate_models", "--manifest", str(manifest), flag,
                                    "--cache-dir", str(tmp_path / "cache"),
                                    "--output-dir", str(output)])
    def fetch(*args, **kwargs):
        assert kwargs["offline"] == (mode != "download")
        if mode == "missing":
            raise validate_models.ArchiveError("not cached")
        return raw

    monkeypatch.setattr(validate_models, "fetch", fetch)
    monkeypatch.setattr(validate_models, "run_case", lambda source, method, path: {
        "source": "btx", "method": method, "passed": mode in ("pass", "download"), "metrics": [
            {"quantity": "pressure", "unit": "%", "points": 1,
             "mean_absolute_error": 1.0, "maximum_absolute_error": 1.0,
             "maximum_error_point": "x=0.5"}]})
    assert validate_models.main() == expected
    report = json.loads((output / "validation.json").read_text())
    if mode in ("bad_hash", "missing"):
        assert not report["results"]
        assert ("SHA-256 mismatch" if mode == "bad_hash" else "not cached") in (
            report["errors"][0]["error"])
    else:
        assert report["results"][0]["passed"] == (mode in ("pass", "download"))
        assert "x=0.5" in (output / "validation.md").read_text()


# --- Feature reference campaign (roadmap G5) --------------------------------------------------

@pytest.fixture(scope="module")
def feature_campaign():
    manifest = json.loads(validate_features.MANIFEST.read_text(encoding="utf-8"))
    validate_features.validate(manifest)
    results, errors = validate_features.run(manifest)
    assert not errors
    return {r["case"]: r for r in results}


def _metrics(result):
    return {m["quantity"]: m for m in result["metrics"]}


def test_ammonia_on_fugacities_against_larson(feature_campaign):
    # The numbers recorded in docs/model-validation.md; a model change that moves them must
    # update the page as well.
    result = feature_campaign["ammonia_larson"]
    m = _metrics(result)
    low = m["equilibrium ammonia mole percent, 10-300 atm"]
    high = m["equilibrium ammonia mole percent, 600-1000 atm"]
    assert low["passed"] and not high["passed"] and not result["passed"]
    assert (low["mean_absolute_error"], low["maximum_absolute_error"]) == pytest.approx(
        (0.258, 1.47), abs=0.01)
    assert (high["mean_absolute_error"], high["maximum_absolute_error"]) == pytest.approx(
        (4.40, 6.54), abs=0.01)
    # The ideal gas is not screened, and it is further from the measurement at every point
    # from 100 atm on.
    assert not m["equilibrium ammonia mole percent, 600-1000 atm, ideal gas"].get(
        "screened", True)
    for p in result["points"]:
        if p["P_atm"] >= 100:
            assert abs(p["property_method"] - p["measured"]) < abs(p["ideal_gas"] - p["measured"])


def test_equilibrium_methyl_acetate_column_against_mekala(feature_campaign):
    result = feature_campaign["methyl_acetate_mekala"]
    m = _metrics(result)
    T, x = m["stage temperature"], m["liquid mole fraction"]
    assert not T["passed"] and not x["passed"]
    assert (T["points"], x["points"]) == (8, 20)
    assert (T["mean_absolute_error"], T["maximum_absolute_error"]) == pytest.approx(
        (4.51, 13.09), abs=0.01)
    assert (x["mean_absolute_error"], x["maximum_absolute_error"]) == pytest.approx(
        (0.0872, 0.2815), abs=1e-3)
    assert x["maximum_error_point"].startswith("stage 30, water")


def test_feature_reference_values_are_consistent():
    manifest = json.loads(validate_features.MANIFEST.read_text(encoding="utf-8"))
    cases = {c["id"]: c for c in manifest["cases"]}
    assert set(cases) == set(validate_features.RUNNERS)
    ammonia = cases["ammonia_larson"]
    pressures = [p for s in ammonia["series"] for p in s["pressures_atm"]]
    for row in ammonia["measured"].values():
        assert sorted(map(int, row)) == sorted(pressures)
        # More ammonia at higher pressure, less at higher temperature.
        assert [row[str(p)] for p in sorted(pressures)] == sorted(row.values())
    for T_low, T_high in pairwise(sorted(ammonia["measured"], key=float)):
        for p in pressures:
            assert ammonia["measured"][T_low][str(p)] > ammonia["measured"][T_high][str(p)]
    column = cases["methyl_acetate_mekala"]
    for row in column["measured_liquid_mole_fraction"].values():
        assert sum(row.values()) == pytest.approx(1.0, abs=2e-5)
    # The column limits are the preregistered equipment-campaign limits.
    assert column["thresholds"] == {"temperature": [2, 5], "composition": [0.02, 0.05]}


@pytest.mark.parametrize("manifest", [
    {"version": 2, "cases": [{"id": "ammonia_larson", "references": ["x"]}]},
    {"version": 1, "cases": []},
    {"version": 1, "cases": [{"id": "unknown", "references": ["x"]}]},
    {"version": 1, "cases": [{"id": "ammonia_larson", "references": []}]},
    {"version": 1, "cases": [{"id": "ammonia_larson", "references": ["x"]}] * 2},
])
def test_invalid_feature_campaign_exits_2(tmp_path, manifest):
    path = tmp_path / "features.json"
    path.write_text(json.dumps(manifest))
    with pytest.raises(SystemExit) as exc:
        validate_features.main(["--manifest", str(path), "--output-dir", str(tmp_path / "o")])
    assert exc.value.code == 2


@pytest.mark.parametrize("outcome, expected", [("pass", 0), ("fail", 1), ("error", 2)])
def test_feature_campaign_exit_status_and_reports(tmp_path, monkeypatch, outcome, expected):
    def runner(case):
        if outcome == "error":
            raise ValueError("did not converge")
        limit = [1.0, 2.0]
        errors = [("a", 0.5), ("b", 0.5 if outcome == "pass" else 3.0)]
        return [validate_features.metric("q", "K", errors, limit)], []

    monkeypatch.setattr(validate_features, "RUNNERS", {"ammonia_larson": runner})
    path = tmp_path / "features.json"
    path.write_text(json.dumps({"version": 1, "cases": [
        {"id": "ammonia_larson", "feature": "G2", "references": ["x"]}]}))
    out = tmp_path / "o"
    assert validate_features.main(["--manifest", str(path), "--output-dir", str(out)]) == expected
    report = json.loads((out / "features.json").read_text())
    text = (out / "features.md").read_text()
    if outcome == "error":
        assert "did not converge" in report["errors"][0]["error"] and "Not run" in text
    else:
        assert report["results"][0]["passed"] == (outcome == "pass")
        assert ("| pass |" if outcome == "pass" else "| fail |") in text


def test_batch_vessel_cost_by_hand_from_turton_table_a1():
    # Turton Table A.1, jacketed agitated reactor: log10 Cp = K1 + K2 log10 V + K3 (log10 V)^2
    # with K = 4.1052, 0.5320, -0.0005 (2001 $, CEPCI 397), typed here from the table rather
    # than imported, at 2 m3 of vessel.
    from fluxsheet.costing import STIRRED

    log_v = math.log10(2.0)
    by_hand = 10 ** (4.1052 + 0.5320 * log_v - 0.0005 * log_v ** 2)
    assert STIRRED.cost(2.0)[0] == pytest.approx(by_hand, rel=1e-12)
