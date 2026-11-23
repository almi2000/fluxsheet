"""Plug-ins (fluxsheet.plugins): unit operations, property methods, compounds and dynamic
models from other packages through entry points, with their failures reported instead of
raised; plug-in units' own cost estimates; Henry components with plug-in activity methods."""

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
from fluxsheet.dynamics import units as dynamic_registry
from fluxsheet.dynamics.simulator import simulate
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
    dynamic = dict(dynamic_registry.DYNAMIC_MODELS)
    plugin_dynamic = dict(dynamic_registry.PLUGIN_MODELS)
    loaded, modules = list(plugins.PLUGINS), set(sys.modules)
    get_package.cache_clear()
    yield
    for d, saved in ((UNIT_TYPES, units), (tp.METHODS, methods),
                     (tp.PLUGIN_METHODS, plugin_methods), (tp.PLUGIN_COMPOUNDS, compounds),
                     (dynamic_registry.DYNAMIC_MODELS, dynamic),
                     (dynamic_registry.PLUGIN_MODELS, plugin_dynamic)):
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


def load_dynamic(value, name="test"):
    return plugins.load("dynamic_models", dynamic_registry._register_plugin,
                        [entry("dynamic_models", value, name)])[0]


def load_example():
    """Everything the example plug-in registers, as if installed."""
    load_units("fluxsheet_example_plugin.units", "example")
    load_compounds("fluxsheet_example_plugin.compounds:COMPOUNDS")
    load_methods("fluxsheet_example_plugin.thermo:REGULAR_SOLUTION")
    load_dynamic("fluxsheet_example_plugin.dynamics:MODELS", "example")


def example(name):
    return json.loads((EXAMPLE / f"{name}.json").read_text(encoding="utf-8"))


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


# ---- dynamic models ------------------------------------------------------------------------


def series(out):
    return {x["key"]: np.array(x["values"]) for x in out["series"]}


@pytest.mark.parametrize("pressure_driven", [False, True])
def test_a_plugin_dynamic_model_follows_its_analytic_step_response(registries, pressure_driven):
    load_example()
    model = example("naphtha_cooler")
    model["dynamics"]["pressure_driven"] = pressure_driven
    out = simulate(model)
    assert out["ok"], out.get("error")
    s, t = series(out), np.array(out["time"])
    # The example's lag: dh/dt = (h_ss - h) / tau with a constant flow, so after the step of
    # the superheat at t0 the duty is Q1 + (Q0 - Q1) exp(-(t - t0) / tau), Q0 and Q1 the
    # steady duties at the old and the new superheat.
    Q0 = solve(model).units["COOLER"]["duty"]
    model["units"][1]["params"]["superheat"] = "15 K"
    after = solve(model)
    Q1 = after.units["COOLER"]["duty"]
    assert Q1 - Q0 > 1e4  # a step of about 52 kW less cooling
    t0, tau = 30.0, 40.0
    expected = np.where(t < t0, Q0, Q1 + (Q0 - Q1) * np.exp(-(t - t0) / tau))
    np.testing.assert_allclose(s["unit.COOLER.duty"], expected, atol=1e-3 * (Q1 - Q0))
    # At the end (6.75 time constants) the outlet is 0.1 % of the step from the new steady
    # state, which lies 15 K above the dew point.
    assert s["stream.COOLED.T"][-1] == pytest.approx(after.streams["COOLED"]["T"], abs=0.02)
    start = solve(example("naphtha_cooler")).streams["COOLED"]["T"]
    assert s["unit.COOLER.T_out"][0] == pytest.approx(start, abs=1e-6)


LAG = "from fluxsheet.dynamics.units import Lag\n"


@pytest.mark.parametrize("source, message", [
    (LAG + "X = {'Heater': Lag}\n", "dynamic model already defined: Heater"),
    (LAG + "X = {'DistillationColumn': Lag}\n", "already defined: DistillationColumn"),
    (LAG + "X = {'NoSuchUnit': Lag}\n", "unknown unit type 'NoSuchUnit'"),
    ("X = {'Mixer': int}\n", "is not a DynamicModel subclass"),
    (LAG + "X = Lag\n", "expected a dict {unit type: DynamicModel subclass}"),
    ("X = {}\n", "provides no dynamic model"),
    ("raise ImportError('no scipy')\n", "ImportError: no scipy"),
])
def test_a_bad_dynamic_model_plugin_is_reported(registries, package_dir, source, message):
    package_dir("plugtest_dyn", source)
    plugin = load_dynamic("plugtest_dyn:X")
    assert message in plugin.error and not plugin.provides
    assert dynamic_registry.DYNAMIC_MODELS["Heater"] is dynamic_registry.Lag
    assert "Mixer" not in dynamic_registry.DYNAMIC_MODELS
    assert plugins.report()[-1]["error"] == plugin.error


def test_a_plugin_may_give_a_builtin_unit_a_model_but_only_once(registries, package_dir):
    package_dir("plugtest_dyn_mixer", LAG + "def models():\n    return {'Mixer': Lag}\n")
    one = load_dynamic("plugtest_dyn_mixer:models", "one")
    assert one.error is None and one.provides == ["Mixer"]
    assert dynamic_registry.DYNAMIC_MODELS["Mixer"] is dynamic_registry.Lag
    assert dynamic_registry.PLUGIN_MODELS == {"Mixer": "one"}
    two = load_dynamic("plugtest_dyn_mixer:models", "two")
    assert "dynamic model already defined: Mixer (plug-in one)" in two.error


def test_the_builtin_dynamic_models_are_reserved():
    # Plug-ins load before dynamics.column registers itself when that module is imported
    # first, so its unit type must be reserved by name.
    builtin = {t for t in dynamic_registry.DYNAMIC_MODELS
               if t not in dynamic_registry.PLUGIN_MODELS}
    assert builtin == set(dynamic_registry.BUILTIN_MODELS)


# ---- costing --------------------------------------------------------------------------------

COSTED = TRIM.replace("ENERGY", '("duty",)') + """
    from fluxsheet import costing

    def _cost(c):
        HOOK

    TrimHeater.estimate_cost = staticmethod(_cost)
"""
DRUM = """c.equipment("drum", costing.VESSEL, 4.0,
                    (2.25, 1.82, costing.FM_VESSEL[c.material]), "a 4 m3 drum")
        c.electricity("stirrer", 2000.0)
        c.note("inlet at %.1f K" % c.temperatures[0])"""
COSTING = {"cepci": 600.0, "material": "SS", "electricity_price": 0.1, "hours_per_year": 8000}


def test_a_plugin_unit_prices_its_equipment_as_the_correlation_by_hand(registries,
                                                                      package_dir):
    package_dir("plugtest_cost", COSTED.replace("HOOK", DRUM))
    load_units("plugtest_cost")
    r = solve(trim_flowsheet(costing=COSTING))
    assert r.ok, r.error
    unit = r.costing["units"]["T"]
    # Turton Table A.1, vertical vessel: log10 Cp0 = 3.4974 + 0.4485 lg V + 0.1074 lg^2 V
    # (2001 $, CEPCI 397); bare module (B1 + B2 FM FP) with B 2.25, 1.82, FM 3.1 (stainless),
    # FP 1 at 1 bar.
    lg = math.log10(4.0)
    purchased = 10 ** (3.4974 + 0.4485 * lg + 0.1074 * lg * lg) * 600.0 / 397.0
    assert unit["purchased_cost"] == pytest.approx(purchased, rel=1e-12)
    assert unit["bare_module_cost"] == pytest.approx(purchased * (2.25 + 1.82 * 3.1), rel=1e-12)
    assert unit["operating_cost"] == pytest.approx(2.0 * 8000 * 0.1)  # 2 kW, 8000 h, 0.1 $/kWh
    assert unit["notes"] == ["inlet at 300.0 K"]
    assert r.units["T"]["bare_module_cost"] == unit["bare_module_cost"]
    assert r.costing["bare_module_cost"] == unit["bare_module_cost"]
    assert r.costing["not_costed"] == []


@pytest.mark.parametrize("hook, outcome", [
    ("pass", "not costed"),
    ("raise KeyError('diameter')",
     "the cost estimate of Trim heater failed: KeyError: 'diameter'"),
    ("c.equipment('drum', costing.VESSEL, -1.0, 4.0, 'x')", "the size to cost must be positive"),
])
def test_a_plugin_cost_hook_that_adds_nothing_or_fails(registries, package_dir, hook, outcome):
    package_dir("plugtest_cost_bad", COSTED.replace("HOOK", hook))
    load_units("plugtest_cost_bad")
    r = solve(trim_flowsheet(costing={}))
    if outcome == "not costed":
        assert r.ok and r.costing["not_costed"] == ["T"] and not r.costing["units"]
        return
    assert not r.ok and outcome in r.error["message"]
    assert r.error["where"] == "T"
    if "KeyError" in outcome:
        assert "from the plug-in 'test'" in r.error["hint"]


def test_a_plugin_unit_without_a_cost_hook_gets_its_utility_only(registries, package_dir):
    package_dir("plugtest_nocost", TRIM.replace("ENERGY", '("duty",)'))
    load_units("plugtest_nocost")
    r = solve(trim_flowsheet(costing={}))
    assert r.ok, r.error
    unit = r.costing["units"]["T"]
    assert not unit["items"] and unit["utilities"][0]["kind"] == "heating"
    assert "equipment not costed; utilities only" in unit["notes"]


def test_example_cooler_costs_its_exchanger_by_hand(registries):
    load_example()
    r = solve(example("naphtha_cooler"))
    assert r.ok, r.error
    c, unit = r.costing, r.costing["units"]["COOLER"]
    Q = r.units["COOLER"]["duty"]
    T_in, T_out = r.streams["HOT"]["T"], r.streams["COOLED"]["T"]
    # Cooling water 30 -> 40 °C countercurrent against the gas cooled from T_in to T_out.
    d1, d2 = T_in - 313.15, T_out - 303.15
    area = abs(Q) / (500.0 * (d1 - d2) / math.log(d1 / d2))  # default U 500 W/m2K
    assert area < 10.0  # below Turton's range: six-tenths rule from 10 m2
    # Fixed-tube-sheet exchanger, Turton Table A.1: K 4.3247, -0.3030, 0.1634 at 10 m2.
    purchased = (10 ** (4.3247 - 0.3030 + 0.1634) * (area / 10.0) ** 0.6
                 * c["cepci"] / 397.0)
    (item,) = unit["items"]
    assert item["size"] == pytest.approx(area, rel=1e-9)
    assert item["purchased_cost"] == pytest.approx(purchased, rel=1e-9)
    # Bare module: B1 1.63, B2 1.66, FM 1 (carbon steel), FP 1 below 5 barg.
    assert item["bare_module_cost"] == pytest.approx(purchased * (1.63 + 1.66), rel=1e-9)
    (cw,) = unit["utilities"]
    assert cw["cost"] == pytest.approx(abs(Q) * c["hours_per_year"] * 3600 * 0.354e-9)


# ---- Henry components with plug-in activity methods --------------------------------------------

GASES = ("hydrogen", "methane", "n-heptane", "toluene")


def _recording(package_dir, henry=True, key="rs-recording"):
    package_dir("plugtest_henry", f"""
        from dataclasses import replace
        from fluxsheet_example_plugin.thermo import REGULAR_SOLUTION, regular_solution
        CALLS = []

        def recording(consts, components):
            CALLS.append(list(components))
            return regular_solution(consts, components)

        X = replace(REGULAR_SOLUTION, key={key!r}, gibbs_excess=recording, henry={henry})
    """)
    assert load_methods("plugtest_henry:X").error is None
    return sys.modules["plugtest_henry"].CALLS


def test_henry_with_a_plugin_method_follows_henrys_law_by_hand(registries, package_dir):
    from thermo.interaction_parameters import IPDB

    calls = _recording(package_dir)
    pkg = get_package(GASES, "rs-recording")
    assert calls == [[2, 3]]  # the solvents only
    assert pkg.henry_components == ["hydrogen", "methane"]
    T, P = 313.15, 10e5
    s = pkg.flash([0.3, 0.1, 0.3, 0.3], T=T, P=P)
    assert 0.2 < s.VF < 0.6
    x, y, c = np.array(s.x), np.array(s.y), pkg.consts
    R = 8.314462618
    # Solvents: modified Raoult's law with the Hildebrand gammas of the gas-free liquid.
    xs = x[2:] / x[2:].sum()
    V, d = np.array(c.Vml_STPs[2:]), np.array(c.solubility_parameters[2:])
    phi = xs * V / (xs @ V)
    gamma = np.exp(V * (d - phi @ d) ** 2 / (R * T))
    psat = np.array([pkg.props.VaporPressures[i](T) for i in (2, 3)])
    np.testing.assert_allclose(y[2:] * P, x[2:] * gamma * psat, rtol=1e-6)
    # Gases: y P = x H, ln H = sum_j x_j ln H_ij / sum_j x_j over the solvents, H_ij from the
    # table thermo ships (Peng-Robinson predictions for these pairs).
    table = IPDB.tables["PRTranslated Henry Best"]
    for i in (0, 1):
        lnH = []
        for j in (2, 3):
            k = table[f"{c.CASs[i]} {c.CASs[j]}"]
            lnH.append(k["A"] + k["B"] / T + k["C"] * math.log(T) + k["D"] * T
                       + k["E"] / T ** 2 + k["F"] * T ** 2)
        H = math.exp((x[2] * lnH[0] + x[3] * lnH[1]) / (x[2] + x[3]))
        assert y[i] * P == pytest.approx(x[i] * H, rel=1e-6)


def test_the_example_hydrogen_separator(registries):
    load_example()
    r = solve(example("hydrogen_separator"))
    assert r.ok, r.error
    assert any("hydrogen, methane are Henry components" in w for w in r.warnings)
    oil = r.streams["L1"]["mole_fractions"]
    assert oil["hydrogen"] < 0.01 and oil["methane"] < 0.01  # dilute: Henry's law holds
    assert not any("makes up more than" in w for w in r.warnings)


NRTL_PLUGIN = """
    from thermo import NRTL
    from fluxsheet.plugins import ActivityMethod
    from fluxsheet.thermo_package import ip_matrix

    def nrtl(consts, idx, every=False):
        idx = list(range(consts.N)) if every else idx
        b = ip_matrix("ChemSep NRTL", consts.CASs, "bij")
        alpha = ip_matrix("ChemSep NRTL", consts.CASs, "alphaij")
        n = len(idx)
        return NRTL(T=298.15, xs=[1.0 / n] * n, tau_as=[[0.0] * n for _ in idx],
                    tau_bs=[[b[i][j] for j in idx] for i in idx],
                    alpha_cs=[[alpha[i][j] for j in idx] for i in idx])

    X = ActivityMethod("plain-nrtl", "ChemSep NRTL", "NRTL with ChemSep parameters.", nrtl,
                       henry=True)
    WRONG = ActivityMethod("wrong-nrtl", "Wrong NRTL", "All compounds despite Henry.",
                           lambda c, idx: nrtl(c, idx, every=True), henry=True)
"""


def test_henry_with_a_plugin_method_matches_the_builtin_nrtl(registries, package_dir):
    package_dir("plugtest_nrtl", NRTL_PLUGIN)
    assert load_methods("plugtest_nrtl:X").error is None
    comps = ("methane", "water", "ethanol")
    own = get_package(comps, "plain-nrtl")
    builtin = get_package(comps, "nrtl", prefer_fitted=False)  # ChemSep parameters only
    assert [p["source"] for p in builtin.parameter_sources if p.get("kind") != "henry"] == [
        "ChemSep"]
    assert own.henry_components == builtin.henry_components == ["methane"]
    zs, spec = [0.02, 0.6, 0.38], {"T": 345.0, "P": 1.0e5}
    a, b = own.flash(zs, **spec), builtin.flash(zs, **spec)
    assert 0.0 < b.VF < 1.0
    assert a.VF == pytest.approx(b.VF, rel=1e-9)
    np.testing.assert_allclose(a.x, b.x, rtol=1e-9)
    np.testing.assert_allclose(a.y, b.y, rtol=1e-9)
    assert a.H == pytest.approx(b.H, rel=1e-9)


def test_a_plugin_model_of_the_wrong_size_is_an_error(registries, package_dir):
    package_dir("plugtest_nrtl_wrong", NRTL_PLUGIN)
    assert load_methods("plugtest_nrtl_wrong:WRONG").error is None
    with pytest.raises(Exception, match="returned a model of 3 compounds for 2"):
        get_package(("methane", "water", "ethanol"), "wrong-nrtl")


def test_a_plugin_method_without_henry_gets_every_compound(registries, package_dir):
    calls = _recording(package_dir, henry=False, key="rs-plain")
    pkg = get_package(("n-hexane", "toluene"), "rs-plain", henry_components=("n-hexane",))
    assert calls == [[0, 1]] and pkg.henry_components == []
    assert any("this plug-in method has no Henry components" in w for w in pkg.warnings)


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


def test_cli_and_api_list_dynamic_models_of_plugins(registries, capsys):
    from fluxsheet.server import app

    load_units("fluxsheet_example_plugin.units", "example")
    load_dynamic("fluxsheet_example_plugin.dynamics:MODELS", "example")
    assert cli(["plugins"]) == 0
    (line,) = [x for x in capsys.readouterr().out.splitlines() if x.startswith("dynamic_models")]
    assert line.split()[1] == "example" and line.endswith("DewPointCooler")
    client = TestClient(app)
    assert client.get("/api/plugins").json()[-1]["provides"] == ["DewPointCooler"]
    types = {t["type"]: t for t in client.get("/api/unit-types").json()}
    dynamic = types["DewPointCooler"]["dynamic"]
    assert [p["name"] for p in dynamic["params"]] == ["time_constant"]
    assert [f["name"] for f in dynamic["result_fields"]] == ["T_out", "duty"]


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
    metadata in a fresh interpreter, run the example flowsheets, the cooler's cost estimate
    and its dynamic model."""
    project = tomllib.loads((EXAMPLE / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    assert set(project["entry-points"]) == set(plugins.GROUPS.values())  # one of each kind
    info = tmp_path / f"fluxsheet_example_plugin-{project['version']}.dist-info"
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
        from fluxsheet.dynamics import simulate
        assert main(["plugins"]) == 0
        out = {{}}
        for name in ("naphtha_cooler", "regular_solution_flash", "hydrogen_separator"):
            with open(r"{EXAMPLE}/" + name + ".json", encoding="utf-8") as f:
                model = json.load(f)
            r = solve(model)
            assert r.ok, r.error
            out[name] = {{"streams": r.streams, "units": r.units, "costing": r.costing,
                          "checks": [c["status"] for c in r.balance_checks]}}
        with open(r"{EXAMPLE}/naphtha_cooler.json", encoding="utf-8") as f:
            dyn = simulate(json.load(f))
        assert dyn["ok"], dyn.get("error")
        out["dynamic"] = {{s["key"]: s["values"] for s in dyn["series"]}}
        print("RESULT" + json.dumps(out))
    """)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(
        [str(ROOT / "src"), str(tmp_path), str(EXAMPLE / "src")]))
    env.pop(plugins.DISABLE_ENV, None)
    done = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True,
                          env=env, timeout=300, encoding="utf-8", check=False)
    assert done.returncode == 0, done.stderr
    listing, result = done.stdout.split("RESULT")
    lines = {line.split()[1]: line for line in listing.splitlines() if line.strip()}
    for name, group, provided in (("example_units", "units", "DewPointCooler"),
                                  ("regular_solution", "property_methods", "regular-solution"),
                                  ("naphtha", "components", "light-naphtha-cut"),
                                  ("cooler_lag", "dynamic_models", "DewPointCooler")):
        assert lines[name].startswith(group) and lines[name].endswith(provided)
        assert f"fluxsheet-example-plugin {project['version']}" in lines[name]
    out = json.loads(result)
    dynamic = out.pop("dynamic")
    for case in out.values():
        assert set(case["checks"]) == {"ok"}
    cooler = out["naphtha_cooler"]
    assert cooler["streams"]["COOLED"]["vapor_fraction"] == 1.0
    assert cooler["units"]["COOLER"]["duty"] < 0
    assert cooler["costing"]["units"]["COOLER"]["items"][0]["item"] == "exchanger"
    assert cooler["costing"]["not_costed"] == []
    flash = out["regular_solution_flash"]
    assert 0.0 < flash["units"]["FLASH"]["vapor_fraction"] < 1.0
    assert out["hydrogen_separator"]["streams"]["L1"]["mole_fractions"]["hydrogen"] < 0.01
    # The superheat steps up at 30 s: the duty moves towards less cooling, and slowly.
    duty = dynamic["unit.COOLER.duty"]
    assert duty[0] == pytest.approx(cooler["units"]["COOLER"]["duty"], rel=1e-6)
    assert duty[0] < duty[10] < duty[-1]


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
