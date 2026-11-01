"""Phase behaviour checked against physics, not against hand-picked numbers.

Example-based tests only cover the conditions someone thought of; a compressed hydrogen feed
came out as a liquid for months because none of them had a light gas at high pressure. The
checks here derive the right answer from first principles over a grid of conditions:

- a pure component is vapor above its critical temperature, below it vapor under its vapor
  pressure and liquid above it (the dense supercritical region, where the name is a matter
  of convention, and the immediate neighbourhood of the saturation curve are left out);
- a mixture is vapor above its dew point, liquid below its bubble point, two-phase between;
- a bubble point is never a dew point in disguise.
"""

import pytest

from fluxsheet.errors import SimulationError
from fluxsheet.thermo_package import PropertyPackage

COMPOUNDS = ["hydrogen", "helium", "nitrogen", "methane", "ethane", "propane", "butane",
             "hexane", "decane", "benzene", "toluene", "biphenyl", "water", "methanol",
             "ethanol", "ammonia", "carbon dioxide", "acetone", "acetic acid"]
REDUCED_T = [0.45, 0.6, 0.75, 0.9, 0.97, 1.1, 1.6, 3.0, 8.0]


def pure_cases(pkg):
    """(T, P, expected phase) for one compound, from Tc, Pc and the vapor pressure."""
    Tc, Pc, Tm = pkg.consts.Tcs[0], pkg.consts.Pcs[0], pkg.consts.Tms[0] or 0.0
    for Tr in REDUCED_T:
        T = Tr * Tc
        if T < Tm + 5.0 or T > 1500.0:
            continue
        if Tr < 1.0:
            ps = pkg.props.VaporPressures[0](T)
            for factor, phase in ((0.3, "vapor"), (3.0, "liquid"), (30.0, "liquid")):
                P = factor * ps
                if 1e3 <= P <= 1e8 and not (P > Pc and Tr > 0.9):
                    yield T, P, phase
        else:
            for Pr in (0.05, 0.5, 2.0, 10.0):
                if Pr * Pc <= 1e8 and not (Pr > 1.0 and Tr < 1.5):
                    yield T, Pr * Pc, "vapor"


@pytest.mark.parametrize("method", ["peng-robinson", "srk", "ideal"])
@pytest.mark.parametrize("compound", COMPOUNDS)
def test_pure_component_phase_follows_critical_point_and_vapor_pressure(compound, method):
    pkg = PropertyPackage([compound], method)
    wrong = [f"T={T:.0f} K, P={P / 1e5:.3g} bar: {phase} instead of {expected}"
             for T, P, expected in pure_cases(pkg)
             if (phase := pkg.flash([1.0], T=T, P=P).phase) != expected]
    assert not wrong, f"{compound} ({method}): " + "; ".join(wrong)


MIXTURES = [
    (("benzene", "toluene"), [0.5, 0.5], [1e5, 5e5, 20e5], ["peng-robinson", "srk", "nrtl"]),
    (("methane", "propane", "butane"), [0.6, 0.25, 0.15], [10e5, 40e5], ["peng-robinson", "srk"]),
    (("ethane", "heptane"), [0.4, 0.6], [5e5, 30e5], ["peng-robinson", "srk"]),
    (("hydrogen", "toluene"), [0.8, 0.2], [1e5, 30e5, 70e5], ["peng-robinson", "srk"]),
    (("hydrogen", "methane", "benzene", "toluene", "biphenyl"),
     [0.65, 0.16, 0.14, 0.045, 0.005], [1e5, 34e5], ["peng-robinson", "srk"]),
    (("nitrogen", "oxygen", "argon"), [0.78, 0.21, 0.01], [1e5, 10e5], ["peng-robinson"]),
    (("carbon dioxide", "methanol"), [0.3, 0.7], [10e5, 50e5], ["peng-robinson"]),
    (("water", "ethanol"), [0.6, 0.4], [1e5, 5e5], ["nrtl", "unifac", "peng-robinson"]),
    (("ammonia", "water"), [0.3, 0.7], [2e5, 15e5], ["peng-robinson", "nrtl"]),
]


def saturation_T(pkg, zs, kind, P):
    try:
        return pkg.saturation_point(zs, kind, P=P).T
    except SimulationError:  # no such point at this pressure (e.g. hydrogen-rich): skip
        return None


@pytest.mark.parametrize("components, zs, pressures, methods", MIXTURES,
                         ids=["/".join(m[0]) for m in MIXTURES])
def test_mixture_phase_agrees_with_its_bubble_and_dew_points(components, zs, pressures, methods):
    wrong = []
    for method in methods:
        pkg = PropertyPackage(list(components), method)
        for P in pressures:
            Tb, Td = saturation_T(pkg, zs, "bubble", P), saturation_T(pkg, zs, "dew", P)
            checks = []
            if Td:
                checks.append((Td * 1.03 + 2.0, "vapor"))
            if Tb:
                checks.append((Tb * 0.97 - 2.0, "liquid"))
            if Tb and Td and Td - Tb > 4.0:
                checks.append((0.5 * (Tb + Td), "two-phase"))
            for T, expected in checks:
                if (phase := pkg.flash(zs, T=T, P=P).phase) != expected:
                    wrong.append(f"{method}, P={P / 1e5:g} bar, T={T:.0f} K: {phase} "
                                 f"instead of {expected}")
    assert not wrong, "; ".join(wrong)


HDA = ["hydrogen", "methane", "benzene", "toluene", "biphenyl"]


@pytest.mark.parametrize("method", ["peng-robinson", "srk"])
def test_bubble_point_is_never_a_dew_point_in_disguise(method):
    # Hydrogen-rich: no bubble point at 1 bar. The solver used to return the dew point.
    pkg = PropertyPackage(HDA, method)
    zs = [0.65, 0.16, 0.14, 0.045, 0.005]
    pkg.saturation_point(zs, "dew", P=1e5)
    with pytest.raises(SimulationError, match="no bubble point|did not converge"):
        pkg.saturation_point(zs, "bubble", P=1e5)


@pytest.mark.parametrize("method", ["peng-robinson", "srk"])
def test_bubble_point_with_dissolved_traces_lies_below_the_dew_point(method):
    # Traces of hydrogen and methane in aromatics at 20 bar: a real bubble point, whose
    # incipient vapor is enriched in the light gases.
    pkg = PropertyPackage(HDA, method)
    zs = [0.0005, 0.002, 0.6, 0.3775, 0.02]
    dew = pkg.saturation_point(zs, "dew", P=20e5)
    bubble = pkg.saturation_point(zs, "bubble", P=20e5)
    assert bubble.T < dew.T - 1.0
    assert max(abs(a - b) for a, b in zip(bubble.x, zs)) < 1e-3
    assert bubble.y[0] > 10 * zs[0] and bubble.y[1] > 10 * zs[1]
