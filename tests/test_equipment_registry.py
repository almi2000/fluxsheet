"""Qualification-only campaigns must never masquerade as successful validation."""

import json
import subprocess
import sys
from pathlib import Path

import pytest
from scripts import validate_models
from scripts.equipment_registry import qualify

REGISTER = Path(__file__).parents[1] / "docs" / "equipment-campaign.json"


def test_registered_sources_are_open():
    records = json.loads(REGISTER.read_text(encoding="utf-8"))["equipment_sources"]
    results = qualify(records)
    assert len(results) == 7
    assert {r["family"] for r in results} == {"dynamics", "kinetics", "column"}
    assert all(r["passed"] is None and r["screening"] == "not_checkable" for r in results)


@pytest.mark.parametrize("field,value", [("admission", "qualified"), ("sha256", "bad"),
    ("missing_inputs", []), ("independence", "independent"), ("family", "other")])
def test_invalid_qualification_rejected(field, value):
    record = json.loads(REGISTER.read_text(encoding="utf-8"))["equipment_sources"][0]
    record[field] = value
    with pytest.raises(ValueError):
        qualify([record])


def test_duplicate_ids_rejected():
    record = json.loads(REGISTER.read_text(encoding="utf-8"))["equipment_sources"][0]
    with pytest.raises(ValueError, match="duplicate"):
        qualify([record], [record["id"]])


def test_offline_qualification_report(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["validate_models", "--manifest", str(REGISTER),
        "--offline", "--cache-dir", str(tmp_path / "cache"),
        "--output-dir", str(tmp_path / "output")])
    def no_fetch(*args, **kwargs):
        pytest.fail("qualification-only campaign must not download or read raw data")
    monkeypatch.setattr(validate_models, "fetch", no_fetch)
    assert validate_models.main() == 2
    report = json.loads((tmp_path / "output" / "validation.json").read_text())
    assert len(report["errors"]) == 7
    assert not report["results"]
    assert len(report["equipment_sources"]) == 7
    markdown = (tmp_path / "output" / "validation.md").read_text()
    assert "not_checkable" in markdown
    assert "Missing inputs:" in markdown
    assert "shared_fit_data" in markdown


def test_v1_cannot_silently_ignore_equipment(tmp_path, monkeypatch):
    data = json.loads(REGISTER.read_text(encoding="utf-8"))
    data["version"] = 1
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(data))
    monkeypatch.setattr(sys, "argv", ["validate_models", "--manifest", str(manifest),
        "--cache-dir", str(tmp_path / "cache"), "--output-dir", str(tmp_path / "output")])
    with pytest.raises(SystemExit) as exc:
        validate_models.main()
    assert exc.value.code == 2


def test_direct_cli_returns_two(tmp_path):
    process = subprocess.run([sys.executable, str(REGISTER.parents[1] / "scripts" /
        "validate_models.py"), "--manifest", str(REGISTER), "--offline", "--cache-dir",
        str(tmp_path / "cache"), "--output-dir", str(tmp_path / "output")],
        capture_output=True, text=True, check=False, timeout=30)
    assert process.returncode == 2, process.stderr
    assert "7 unavailable or invalid" in process.stdout


def test_mixed_campaign_preserves_comparisons(tmp_path, monkeypatch):
    data = json.loads(REGISTER.read_text(encoding="utf-8"))
    raw = tmp_path / "raw.xml"
    raw.write_bytes(b"synthetic")
    data["sources"] = [{"id": "btx", "doi": "10.1234/test", "sha256": validate_models.digest(raw),
                        "methods": ["ideal"], "thresholds": [[5, 10], [5, 10]]}]
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(data))
    monkeypatch.setattr(sys, "argv", ["validate_models", "--manifest", str(manifest),
        "--cache-dir", str(tmp_path / "cache"), "--output-dir", str(tmp_path / "output")])
    monkeypatch.setattr(validate_models, "fetch", lambda *a, **kw: raw)
    monkeypatch.setattr(validate_models, "run_case", lambda *a: {
        "source": "btx", "method": "ideal", "passed": False, "metrics": []})
    assert validate_models.main() == 2
    report = json.loads((tmp_path / "output" / "validation.json").read_text())
    assert report["results"][0]["passed"] is False
    assert len(report["equipment_sources"]) == 7
