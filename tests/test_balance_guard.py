"""The test-suite guard must reject a solver error even through an imported alias."""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.parametrize("failure", ["mass", "energy"])
def test_unexpected_balance_error_fails_pytest(tmp_path, failure):
    shutil.copyfile(Path(__file__).with_name("conftest.py"), tmp_path / "conftest.py")
    (tmp_path / "test_injected.py").write_text('''
from fluxsheet import solve as imported_solve
from fluxsheet.stream import Stream
from fluxsheet.units.mixers import Splitter

def test_apparently_successful_solve(monkeypatch):
    original = Splitter.calculate
    def broken(self, pkg, inlets):
        outputs = original(self, pkg, inlets)
        if FAILURE == "mass":
            outputs["out1"].flows[0] *= 0.5
        else:
            stream = outputs["out1"]
            outputs["out1"] = Stream.from_flash(pkg, stream.flows, T=stream.T + 10, P=stream.P)
        return outputs
    monkeypatch.setattr(Splitter, "calculate", broken)
    model = {"components": ["propane"], "property_method": "peng-robinson",
             "units": [{"id": "F", "type": "Feed", "params": {
                 "T": 300, "P": "2 bar", "flows": {"propane": 1}}},
                       {"id": "S", "type": "Splitter", "params": {"split": 0.5}}],
             "streams": [{"id": "feed", "from": "F.out", "to": "S.in"},
                         {"id": "A", "from": "S.out1"},
                         {"id": "B", "from": "S.out2"}]}
    assert imported_solve(model).ok
'''.replace("FAILURE", repr(failure)), encoding="utf-8")
    result = subprocess.run([sys.executable, "-m", "pytest", "-q", "--tb=short", str(tmp_path)],
                            capture_output=True, text=True, timeout=90, check=False)
    assert result.returncode == 1, result.stdout + result.stderr
    assert "unexpected balance diagnostic" in result.stdout
    assert "teardown" in result.stdout
