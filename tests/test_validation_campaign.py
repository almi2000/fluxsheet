"""Offline campaign status and structured summary checks use synthetic input."""

import hashlib
import json
import sys

import pytest
from scripts import validate_models


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
