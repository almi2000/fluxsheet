import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from fluxsheet.server import app

client = TestClient(app)
EXAMPLE = Path(__file__).parent.parent / "examples" / "ethanol_flash.json"


def test_compound_search():
    r = client.get("/api/compounds", params={"q": "toluene"})
    assert r.status_code == 200
    assert r.json()[0]["cas"] == "108-88-3"


def test_unit_types_describe_ports_and_params():
    types = {t["type"]: t for t in client.get("/api/unit-types").json()}
    assert {"Feed", "Heater", "Flash", "Mixer"} <= set(types)
    assert [p["name"] for p in types["Flash"]["outlets"]] == ["vapor", "liquid"]


def test_suggest_method():
    r = client.post("/api/methods/suggest", json={"components": ["water", "acetone"]})
    assert r.json()["method"] == "unifac"


def test_solve_example():
    r = client.post("/api/solve", json=json.loads(EXAMPLE.read_text()))
    body = r.json()
    assert body["ok"], body["error"]
    assert set(body["streams"]) == {"S1", "S2", "S3", "S4"}


def test_solve_reports_unknown_compound():
    model = json.loads(EXAMPLE.read_text())
    model["components"] = ["water", "unobtainium"]
    body = client.post("/api/solve", json=model).json()
    assert not body["ok"] and "unobtainium" in body["error"]["message"]


def test_editor_index_is_revalidated_and_assets_are_cached():
    from fluxsheet.server import _dist

    if not (_dist / "index.html").is_file():
        pytest.skip("web editor not built")
    assert client.get("/").headers["cache-control"] == "no-cache"
    asset = next((_dist / "assets").iterdir()).name
    assert "immutable" in client.get(f"/assets/{asset}").headers["cache-control"]
