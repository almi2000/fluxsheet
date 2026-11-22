"""Plug-ins (fluxsheet.plugins): unit operations, property methods and compounds from other
packages through entry points, with their failures reported instead of raised."""

import json
import math
import os
import subprocess
import sys
import textwrap
import tomllib
from importlib.metadata import EntryPoint
from pathlib import Path

import numpy as np
import pytest
from fastapi.testclient import TestClient

from fluxsheet import plugins, solve
from fluxsheet import thermo_package as tp
from fluxsheet import units as unit_registry
from fluxsheet.cli import main as cli
from fluxsheet.compounds import lookup, search
from fluxsheet.diagrams import phase_envelope
from fluxsheet.thermo_package import get_package, is_eos
from fluxsheet.units import UNIT_TYPES

ROOT = Path(__file__).parent.parent
EXAMPLE = ROOT / "examples" / "plugin"
sys.path.insert(0, str(EXAMPLE / "src"))
from fluxsheet_example_plugin import compounds as example_compounds

# ---- fixtures ---------------------------------------------------------------------------------


@pytest.fixture
def registries():
    """Undo what a test's plug-ins registered."""
    units, methods = dict(UNIT_TYPES), dict(tp.METHODS)
    plugin_methods, compounds = dict(tp.PLUGIN_METHODS), dict(tp.PLUGIN_COMPOUNDS)
    loaded, modules = list(plugins.PLUGINS), set(sys.modules)
    get_package.cache_clear()
    yield
    for d, saved in ((UNIT_TYPES, units), (tp.METHODS, methods),
                     (tp.PLUGIN_METHODS, plugin_methods), (tp.PLUGIN_COMPOUNDS, compounds)):
        d.clear()
        d.update(saved)
    plugins.PLUGINS[:] = loaded
    lookup.cache_clear()  # it remembers plug-in compounds by name
    for name in set(sys.modules) - modules:
        if name.startswith("plugtest_"):
            del sys.modules[name]
    get_package.cache_clear()


@pytest.fixture
def package_dir(tmp_path, monkeypatch):
    """Write modules of a would-be plug-in package and put them on sys.path."""
    monkeypatch.syspath_prepend(str(tmp_path))

    def write(name, source):
        path = tmp_path / f"{name}.py"
        path.write_text(textwrap.dedent(source), encoding="utf-8")
        return name

    return write


def entry(group, value, name="test"):
    return EntryPoint(name, value, plugins.GROUPS[group])


def load_units(value, name="test"):
    return plugins.load("units", unit_registry._register_plugin, [entry("units", value, name)])[0]


def load_methods(value):
    return plugins.load("property_methods", tp._register_methods,
                        [entry("property_methods", value)])[0]


def load_compounds(value):
    return plugins.load("components", tp._register_compounds, [entry("components", value)])[0]


TRIM = """
    from fluxsheet.stream import Stream
    from fluxsheet.units.base import DUTY, IN, OUT, Param, UnitOperation, single_inlet

    class TrimHeater(UnitOperation):
        type = "TrimHeater"
        label = "Trim heater"
        category = "Heat Exchangers"
        description = "Raises the temperature by a fixed amount."
        inlets = (IN,)
        outlets = (OUT,)
        params = (Param("rise", "Temperature rise", "temperature_difference", default=10.0),)
        result_fields = (DUTY,)
        energy_results = ENERGY

        def calculate(self, pkg, inlets):
            s = single_inlet(inlets)
            out = Stream.from_flash(pkg, s.flows, T=s.T + self.p["rise"], P=s.P)
            self.results["duty"] = out.enthalpy_flow - s.enthalpy_flow
            return {"out": out}
"""


def trim_flowsheet(**extra):
    return {
        "components": ["water"], "property_method": "ideal",
        "units": [{"id": "F", "type": "Feed",
                   "params": {"T": "300 K", "P": "1 bar", "flows": {"water": "1 mol/s"}}},
                  {"id": "T", "type": "TrimHeater", "params": {"rise": "15 K"}},
                  {"id": "P", "type": "Product"}],
        "streams": [{"id": "S1", "from": "F.out", "to": "T.in"},
                    {"id": "S2", "from": "T.out", "to": "P.in"}],
        **extra,
    }


# ---- units ------------------------------------------------------------------------------------


def test_a_module_of_units_is_registered_and_solves(registries, package_dir):
    package_dir("plugtest_trim", TRIM.replace("ENERGY", '("duty",)'))
    plugin = load_units("plugtest_trim")
    assert plugin.error is None and plugin.provides == ["TrimHeater"]
    cls = UNIT_TYPES["TrimHeater"]
    assert cls.plugin == "test" and cls.schema()["plugin"] == "test"
    assert UNIT_TYPES["Heater"].schema()["plugin"] is None
    r = solve(trim_flowsheet())
    assert r.ok, r.error
    assert r.streams["S2"]["T"] == pytest.approx(315.0)
    assert r.units["T"]["duty"] > 0
    # The declared energy terms close the unit's energy balance.
    energy = {c["scope"]: c["status"] for c in r.balance_checks if c["quantity"] == "energy"}
    assert energy == {"T": "ok", "flowsheet": "ok"}


def test_a_unit_without_energy_terms_is_not_energy_checkable(registries, package_dir,
                                                             expected_balance_exceptions):
    package_dir("plugtest_trim2", TRIM.replace("ENERGY", "None"))
    load_units("plugtest_trim2")
    for scope in ("T", "flowsheet"):
        expected_balance_exceptions[(scope, "energy", None, "not_checkable")] = (
            "the plug-in does not declare its energy terms")
    assert solve(trim_flowsheet()).ok


@pytest.mark.parametrize("value", ["plugtest_cls:TrimHeater", "plugtest_cls:LIST",
                                   "plugtest_cls:provide"])
def test_a_class_a_list_or_a_function_gives_units(registries, package_dir, value):
    package_dir("plugtest_cls", TRIM.replace("ENERGY", '("duty",)') + """
    LIST = [TrimHeater]

    def provide():
        return [TrimHeater]
    """)
    assert load_units(value).provides == ["TrimHeater"]


@pytest.mark.parametrize("source, message", [
    ("""
     from fluxsheet.units.heat import Heater
     class MyHeater(Heater):
         type = "Heater"
     """, "unit type already defined: Heater"),
    ("x = 1\n", "defines no unit operation"),
    ("raise ImportError('needs numpy 3')\n", "ImportError: needs numpy 3"),
])
def test_a_broken_unit_plugin_is_reported_and_skipped(registries, package_dir, source, message):
    builtin = UNIT_TYPES["Heater"]
    package_dir("plugtest_bad", source)
    plugin = load_units("plugtest_bad")
    assert message in plugin.error and not plugin.provides
    assert UNIT_TYPES["Heater"] is builtin
    assert plugins.report()[-1]["error"] == plugin.error


def test_two_plugins_cannot_share_a_unit_type(registries, package_dir):
    package_dir("plugtest_one", TRIM.replace("ENERGY", "()"))
    package_dir("plugtest_two", TRIM.replace("ENERGY", "()"))
    assert load_units("plugtest_one", "one").error is None
    assert "already defined: TrimHeater" in load_units("plugtest_two", "two").error
    assert UNIT_TYPES["TrimHeater"].__module__ == "plugtest_one"


def test_a_class_that_is_no_unit_operation_is_rejected(registries, package_dir):
    package_dir("plugtest_notunit", "class Thing:\n    type = 'Thing'\n")
    assert "expected a module, UnitOperation" in load_units("plugtest_notunit:Thing").error


# ---- property methods ---------------------------------------------------------------------------

AROMATICS = ["n-hexane", "benzene", "toluene"]


def test_an_activity_method_runs_through_the_activity_liquid(registries):
    plugin = load_methods("fluxsheet_example_plugin.thermo:REGULAR_SOLUTION")
    assert plugin.provides == ["regular-solution"]
    method = tp.METHODS["regular-solution"]
    assert method.plugin == "test" and method.__dict__["plugin"] == "test"
    assert not is_eos("regular-solution")
    pkg = get_package(tuple(AROMATICS), "regular-solution")
    s = pkg.flash([0.2, 0.3, 0.5], VF=0.5, P=101325.0)
    assert s.VF == pytest.approx(0.5) and 350.0 < s.T < 385.0
    # Modified Raoult's law with gammas computed independently from the Hildebrand model:
    # ln gamma_i = V_i (delta_i - delta_mix)^2 / RT.
    c = pkg.consts
    x, y = np.array(s.x), np.array(s.y)
    V, d = np.array(c.Vml_STPs), np.array(c.solubility_parameters)
    phi = x * V / (x @ V)
    gamma = np.exp(V * (d - phi @ d) ** 2 / (8.314462618 * s.T))
    psat = np.array([vp(s.T) for vp in pkg.props.VaporPressures])
    np.testing.assert_allclose(y * 101325.0, x * gamma * psat, rtol=1e-6)


def test_the_example_regular_solution_needs_its_data(registries):
    load_methods("fluxsheet_example_plugin.thermo:REGULAR_SOLUTION")
    load_compounds("fluxsheet_example_plugin.compounds:COMPOUNDS")
    model = {"components": ["n-hexane", "light-naphtha-cut"],
             "property_method": "regular-solution",
             "units": [{"id": "F", "type": "Feed", "params": {
                 "T": "300 K", "P": "1 bar",
                 "flows": {"n-hexane": "1 mol/s", "light-naphtha-cut": "1 mol/s"}}}],
             "streams": [{"id": "S", "from": "F.out"}]}
    r = solve(model)
    assert not r.ok
    assert "no solubility parameter or liquid volume for light-naphtha-cut" in r.error["message"]
    assert "fluxsheet-example-plugin" not in r.error["hint"] and "'test'" in r.error["hint"]


def test_an_equation_of_state_plugin_matches_the_builtin_phase_equilibrium(registries,
                                                                           package_dir):
    package_dir("plugtest_eos", """
        from thermo import PRMIX
        from fluxsheet.plugins import EquationOfState
        from fluxsheet.thermo_package import ip_matrix

        PR = EquationOfState("plain-pr", "Plain Peng-Robinson", "PR without volume shift.",
                             PRMIX, kijs=lambda c: ip_matrix("ChemSep PR", c.CASs, "kij"))
    """)
    plugin = load_methods("plugtest_eos:PR")
    assert plugin.provides == ["plain-pr"], plugin.error
    assert is_eos("plain-pr")
    comps = ("methane", "propane", "n-pentane")
    zs, spec = [0.3, 0.4, 0.3], {"T": 300.0, "P": 20e5}
    own, builtin = get_package(comps, "plain-pr"), get_package(comps, "peng-robinson")
    a, b = own.flash(zs, **spec), builtin.flash(zs, **spec)
    # The built-in method's volume translation changes no phase equilibrium.
    assert a.VF == pytest.approx(b.VF, rel=1e-8)
    np.testing.assert_allclose(a.y, b.y, rtol=1e-8)
    assert own._shifts is None and builtin._shifts is not None
    # An equation of state may draw a phase envelope; the activity models may not.
    curves = phase_envelope(list(comps), zs, "plain-pr", points=6)
    assert curves["bubble"]


def test_kijs_of_the_wrong_shape_are_an_error(registries, package_dir):
    package_dir("plugtest_kij", """
        from thermo import PRMIX
        from fluxsheet.plugins import EquationOfState
        BAD = EquationOfState("bad-pr", "Bad PR", "", PRMIX, kijs=lambda c: [[0.0]])
    """)
    load_methods("plugtest_kij:BAD")
    with pytest.raises(Exception, match="kijs must be a 2 x 2 matrix"):
        get_package(("methane", "ethane"), "bad-pr")


@pytest.mark.parametrize("source, message", [
    (("from fluxsheet_example_plugin.thermo import REGULAR_SOLUTION as R\n"
      "from dataclasses import replace\nX = replace(R, key='nrtl')\n"),
     "property method already defined: nrtl"),
    ("X = {'key': 'dict-method'}\n", "expected an ActivityMethod or EquationOfState"),
])
def test_a_bad_method_plugin_is_reported(registries, package_dir, source, message):
    package_dir("plugtest_meth", source)
    assert message in load_methods("plugtest_meth:X").error
    assert tp.METHODS["nrtl"].plugin is None


# ---- compounds --------------------------------------------------------------------------------


def test_plugin_compounds_are_found_and_simulated(registries):
    plugin = load_compounds("fluxsheet_example_plugin.compounds:COMPOUNDS")
    assert plugin.provides == ["light-naphtha-cut"]
    found = lookup("Light-Naphtha-Cut")
    assert found.name == "light-naphtha-cut" and found.formula == "C8H18"
    assert found.cas.startswith("custom:")
    assert math.isclose(found.mw, 114.23, rel_tol=1e-3)
    assert search("light-nap")[0].name == "light-naphtha-cut"
    model = {"components": ["n-octane", "light-naphtha-cut"], "property_method": "peng-robinson",
             "units": [{"id": "F", "type": "Feed", "params": {
                 "VF": 0.0, "P": "1 atm",
                 "flows": {"n-octane": "1 mol/s", "light-naphtha-cut": "1 mol/s"}}},
                       {"id": "P", "type": "Product"}],
             "streams": [{"id": "S", "from": "F.out", "to": "P.in"}]}
    r = solve(model)
    assert r.ok, r.error
    # Its given boiling point (120 °C) lies below n-octane's (125.6 °C).
    T = r.streams["S"]["T"]
    assert 120 + 273.15 < T < 125.7 + 273.15


def test_a_flowsheet_definition_wins_over_the_plugin(registries):
    load_compounds("fluxsheet_example_plugin.compounds:COMPOUNDS")
    own = dict(example_compounds.COMPOUNDS[0], formula="C9H20",
               joback_groups={"-CH3": 2, "-CH2-": 7})
    custom = (tp.CustomCompound.from_dict(own),)
    (compound,) = tp.resolve(["light-naphtha-cut"], custom)
    assert compound.formula == "C9H20"
    (compound,) = tp.resolve(["light-naphtha-cut"])
    assert compound.formula == "C8H18"


@pytest.mark.parametrize("compounds, message", [
    ([{"name": "ethanol", "formula": "C2H6O"}], "'ethanol' is already in the compound database"),
    ([{"name": "x-cut", "formula": "C8H18"}] * 2, "compound already defined by a plug-in"),
    ([{"name": "no-formula"}], "has no formula"),
])
def test_a_bad_compound_plugin_is_reported(registries, package_dir, compounds, message):
    package_dir("plugtest_comp", f"X = {compounds!r}\n")
    plugin = load_compounds("plugtest_comp:X")
    assert message in plugin.error
    assert not set(tp.PLUGIN_COMPOUNDS) - {"light-naphtha-cut"}


# ---- discovery, CLI and API ----------------------------------------------------------------------


def test_installed_entry_points_are_read_and_can_be_switched_off(monkeypatch):
    calls = []

    def fake(group):
        calls.append(group)
        return [EntryPoint("b", "m:b", group), EntryPoint("a", "m:a", group)]

    monkeypatch.setattr(plugins.metadata, "entry_points", lambda group: fake(group))
    assert [ep.name for ep in plugins.entry_points("units")] == ["a", "b"]
    assert calls == ["fluxsheet.units"]
    for value in ("1", "true", "yes"):
        monkeypatch.setenv(plugins.DISABLE_ENV, value)
        assert plugins.disabled() and plugins.entry_points("units") == []
    monkeypatch.setenv(plugins.DISABLE_ENV, "0")
    assert not plugins.disabled()


def test_cli_lists_plugins_and_fails_on_a_broken_one(registries, package_dir, capsys):
    package_dir("plugtest_cli", TRIM.replace("ENERGY", "()"))
    load_units("plugtest_cli", "trim")
    assert cli(["plugins"]) == 0
    assert "TrimHeater" in capsys.readouterr().out
    package_dir("plugtest_cli_bad", "raise RuntimeError('boom')\n")
    load_units("plugtest_cli_bad", "broken")
    assert cli(["plugins", "--json"]) == 1
    listed = {p["name"]: p for p in json.loads(capsys.readouterr().out)}
    assert listed["broken"]["error"] == "RuntimeError: boom"
    assert listed["trim"]["provides"] == ["TrimHeater"]


def test_server_reports_plugins_units_and_methods(registries, package_dir):
    from fluxsheet.server import app

    package_dir("plugtest_api", TRIM.replace("ENERGY", "()"))
    load_units("plugtest_api", "trim")
    load_methods("fluxsheet_example_plugin.thermo:REGULAR_SOLUTION")
    load_compounds("fluxsheet_example_plugin.compounds:COMPOUNDS")
    client = TestClient(app)
    listed = client.get("/api/plugins").json()[-3:]
    assert [(p["group"], p["provides"]) for p in listed] == [
        ("units", ["TrimHeater"]), ("property_methods", ["regular-solution"]),
        ("components", ["light-naphtha-cut"])]
    types = {t["type"]: t for t in client.get("/api/unit-types").json()}
    assert types["TrimHeater"]["plugin"] == "trim" and types["Flash"]["plugin"] is None
    methods = {m["key"]: m for m in client.get("/api/methods").json()}
    assert methods["regular-solution"]["plugin"] == "test" and methods["nrtl"]["plugin"] is None
    assert client.get("/api/compounds/light-naphtha-cut").json()["formula"] == "C8H18"
    assert client.get("/api/compounds", params={"q": "light-n"}).json()[0]["name"] == (
        "light-naphtha-cut")


# ---- the example plug-in, installed --------------------------------------------------------------


def test_example_plugin_as_installed(tmp_path):
    """The example's pyproject entry points, discovered from an installed distribution's
    metadata in a fresh interpreter, run both example flowsheets."""
    project = tomllib.loads((EXAMPLE / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    info = tmp_path / "fluxsheet_example_plugin-0.1.0.dist-info"
    info.mkdir()
    (info / "METADATA").write_text(
        f"Metadata-Version: 2.1\nName: {project['name']}\nVersion: {project['version']}\n")
    (info / "entry_points.txt").write_text("".join(
        f"[{group}]\n" + "".join(f"{k} = {v}\n" for k, v in eps.items()) + "\n"
        for group, eps in project["entry-points"].items()))
    script = textwrap.dedent(f"""
        import json, sys
        from fluxsheet import plugins, solve
        from fluxsheet.cli import main
        assert main(["plugins"]) == 0
        out = {{}}
        for name in ("naphtha_cooler", "regular_solution_flash"):
            with open(r"{EXAMPLE}/" + name + ".json", encoding="utf-8") as f:
                r = solve(json.load(f))
            assert r.ok, r.error
            out[name] = {{"streams": r.streams, "units": r.units,
                          "checks": [c["status"] for c in r.balance_checks]}}
        print("RESULT" + json.dumps(out))
    """)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(
        [str(ROOT / "src"), str(tmp_path), str(EXAMPLE / "src")]))
    env.pop(plugins.DISABLE_ENV, None)
    done = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True,
                          env=env, timeout=300, encoding="utf-8", check=False)
    assert done.returncode == 0, done.stderr
    listing, result = done.stdout.split("RESULT")
    for provided in ("DewPointCooler", "regular-solution", "light-naphtha-cut"):
        assert provided in listing
    assert "fluxsheet-example-plugin 0.1.0" in listing
    out = json.loads(result)
    for case in out.values():
        assert set(case["checks"]) == {"ok"}
    cooler = out["naphtha_cooler"]
    assert cooler["streams"]["COOLED"]["vapor_fraction"] == 1.0
    assert cooler["units"]["COOLER"]["duty"] < 0
    flash = out["regular_solution_flash"]
    assert 0.0 < flash["units"]["FLASH"]["vapor_fraction"] < 1.0


def test_example_dew_point_cooler_stops_five_kelvin_above_the_dew_point(registries):
    load_units("fluxsheet_example_plugin.units", "example")
    load_compounds("fluxsheet_example_plugin.compounds:COMPOUNDS")
    model = json.loads((EXAMPLE / "naphtha_cooler.json").read_text(encoding="utf-8"))
    r = solve(model)
    assert r.ok, r.error
    out = r.streams["COOLED"]
    pkg = get_package(tuple(model["components"]), "peng-robinson")
    zs = [out["mole_fractions"][c] for c in model["components"]]
    dew = pkg.saturation_point(zs, "dew", P=2.5e5).T
    assert out["T"] == pytest.approx(dew + 5.0, abs=1e-6)
    assert out["P"] == pytest.approx(2.5e5)
