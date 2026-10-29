"""Chemical reactions: stoichiometry, rate laws and equilibrium constants.

Reactions are defined once per flowsheet and referenced by id from reactor units::

    "reactions": [
      {"id": "SMR",
       "stoichiometry": {"methane": -1, "water": -1, "carbon monoxide": 1, "hydrogen": 3},
       "kinetics": {"k0": 1.2e9, "Ea": "240 kJ/mol", "orders": {"methane": 1},
                    "equilibrium_limited": true}}
    ]

Stoichiometric coefficients are negative for reactants and positive for products. Every
reaction must conserve each element; this is checked against the compound database.

Rates are per unit reactor volume and per unit stoichiometric coefficient: component i is
formed at ``nu_i * r`` in mol/(m3 s). The power law is ``r = k(T) prod(c_i^a_i)`` with
``k(T) = k0 exp(-Ea/(R T))``, or ``k(T) = k0 exp(-Ea/R (1/T - 1/T_ref))`` when ``T_ref`` is
given. On the ``"concentration"`` basis c_i is in mol/m3, on the ``"partial_pressure"``
basis it is p_i in Pa; k0 carries whatever units make r come out in mol/(m3 s). Orders
default to the reactant coefficients (elementary reaction).

A reversible reaction either has its own reverse rate law (``r = r_fwd - r_rev``, reverse
orders default to the product coefficients) or is ``equilibrium_limited``:
``r = r_fwd (1 - Q/K)``, with Q the activity quotient (gas: a_i = y_i P/P°, P° = 1 bar;
liquid: ideal solution, a_i = x_i).

Equilibrium constants come from the user correlation ``ln K = A + B/T + C ln T + D T`` when
given. Otherwise they follow from ideal-gas formation data: ΔG° and ΔH° at 298.15 K, carried
to T by integrating the van 't Hoff equation with ideal-gas heat capacities (ΔH°(T) by
Kirchhoff's law), i.e. ``ln K = -ΔG°(T)/(R T)``. This K refers to the ideal-gas standard
state; liquid-phase equilibria need a user correlation.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog

from .errors import SimulationError
from .quantities import UnitError, to_si
from .thermo_package import PropertyPackage

R = 8.314462618  # J/(mol K)
T0 = 298.15  # K, reference temperature of the formation data
P_STD = 1e5  # Pa, standard pressure of the ideal-gas standard state
BASES = ("concentration", "partial_pressure")


@dataclass(frozen=True)
class RateLaw:
    k0: float
    Ea: float  # J/mol
    orders: tuple[float, ...]  # per component
    T_ref: float | None = None

    def k(self, T: float) -> float:
        if self.T_ref is None:
            return self.k0 * math.exp(-self.Ea / (R * T))
        return self.k0 * math.exp(-self.Ea / R * (1.0 / T - 1.0 / self.T_ref))

    def rate(self, T: float, c: np.ndarray) -> float:
        r = self.k(T)
        for ci, a in zip(c, self.orders):
            if a > 0.0:
                r *= max(ci, 0.0) ** a
            elif a < 0.0:
                r *= max(ci, 1e-30) ** a
        return r


@dataclass(frozen=True)
class Kinetics:
    forward: RateLaw
    basis: str  # "concentration" or "partial_pressure"
    reverse: RateLaw | None = None
    equilibrium_limited: bool = False


@dataclass(frozen=True)
class Reaction:
    id: str
    nu: np.ndarray  # stoichiometric coefficient per flowsheet component
    names: tuple[str, ...]
    kinetics: Kinetics | None = None
    keq: tuple[float, float, float, float] | None = None  # A, B, C, D
    order: tuple[int, ...] = ()  # component indices in the order they were written

    @property
    def involved(self) -> list[int]:
        return [i for i, v in enumerate(self.nu) if v != 0.0]

    @property
    def first_reactant(self) -> int:
        return next(i for i in (self.order or range(len(self.nu))) if self.nu[i] < 0)

    @property
    def equation(self) -> str:
        def side(sign: float) -> str:
            terms = []
            for i in self.order or range(len(self.nu)):
                name, v = self.names[i], self.nu[i]
                if v * sign > 0:
                    a = abs(v)
                    terms.append(name if a == 1 else f"{a:g} {name}")
            return " + ".join(terms)

        return f"{side(-1)} -> {side(1)}"

    def ln_keq(self, pkg: PropertyPackage, T: float) -> float:
        """ln K at T, from the user correlation or from ideal-gas formation data."""
        if self.keq is not None:
            A, B, C, D = self.keq
            return A + B / T + C * math.log(T) + D * T
        g = gibbs_over_rt(pkg, T, self.involved, why=f"the equilibrium constant of {self.id}")
        return -float(sum(self.nu[i] * g[i] for i in self.involved))

    def heat_of_reaction(self, pkg: PropertyPackage, T: float = T0) -> float:
        """Ideal-gas standard heat of reaction at T, J/mol of extent."""
        h = 0.0
        for i in self.involved:
            hf = pkg.consts.Hfgs[i]
            if hf is None:
                raise SimulationError(f"no formation enthalpy for {self.names[i]!r}",
                                      where=f"reaction {self.id}")
            h += self.nu[i] * (hf + _cp_integral(pkg, i, T0, T))
        return h

    def rate(self, pkg: PropertyPackage, T: float, P: float, zs: np.ndarray, Vm: float,
             phase: str) -> float:
        """Reaction rate in mol/(m3 s) at T, P and mole fractions zs of a phase with molar
        volume Vm (m3/mol)."""
        kin = self.kinetics
        assert kin is not None
        c = zs / Vm if kin.basis == "concentration" else zs * P
        r = kin.forward.rate(T, c)
        if kin.reverse is not None:
            r -= kin.reverse.rate(T, c)
        elif kin.equilibrium_limited and r != 0.0:
            ln_p = math.log(P / P_STD) if phase == "gas" else 0.0
            ln_q = 0.0
            for i in self.involved:
                if zs[i] <= 0.0 and self.nu[i] > 0:
                    return r  # no product yet: Q = 0
                ln_q += self.nu[i] * (math.log(max(zs[i], 1e-30)) + ln_p)
            r *= 1.0 - math.exp(min(ln_q - self.ln_keq(pkg, T), 50.0))
        return r


class ReactionSet:
    """The reactions of a flowsheet, by id."""

    def __init__(self, reactions: Sequence[Reaction] = ()):
        self.by_id = {r.id: r for r in reactions}

    def __len__(self) -> int:
        return len(self.by_id)

    def select(self, ids: list[str] | None, where: str) -> list[Reaction]:
        available = ", ".join(self.by_id) or "none defined yet"
        if not ids:
            raise SimulationError(
                "no reactions selected", where=where,
                hint=f"list the reaction ids this reactor runs (available: {available})",
            )
        out = []
        for rid in ids:
            if rid not in self.by_id:
                raise SimulationError(
                    f"unknown reaction {rid!r}", where=where,
                    hint=f"define it in the flowsheet's reactions, or use one of: {available}",
                )
            if any(r.id == rid for r in out):
                raise SimulationError(f"reaction {rid!r} is listed twice", where=where,
                                      hint="list each reaction once")
            out.append(self.by_id[rid])
        return out


# ---- parsing and validation ----------------------------------------------------------------

def parse_reactions(models, pkg: PropertyPackage) -> ReactionSet:
    """Validate the flowsheet's reaction definitions (``ReactionModel``s) against ``pkg``."""
    reactions: list[Reaction] = []
    for m in models or ():
        where = f"reaction {m.id}"
        if not m.id or m.id != m.id.strip():
            raise SimulationError(f"reaction id {m.id!r} must be non-empty without edge spaces",
                                  hint="name it e.g. R1")
        if any(r.id == m.id for r in reactions):
            raise SimulationError("two reactions share this id", where=where,
                                  hint="give every reaction a unique id")
        nu = np.zeros(pkg.n)
        order = []
        for name, coeff in m.stoichiometry.items():
            if name not in pkg.names:
                raise SimulationError(
                    f"{name!r} is not in the flowsheet's component list", where=where,
                    hint=f"add it to the components, or use one of: {', '.join(pkg.names)}",
                )
            nu[pkg.names.index(name)] = float(coeff)
            order.append(pkg.names.index(name))
        if not (nu < 0).any() or not (nu > 0).any():
            raise SimulationError(
                "a reaction needs at least one reactant and one product", where=where,
                hint="give reactants negative and products positive coefficients, "
                "e.g. {\"A\": -1, \"B\": 1}",
            )
        _check_element_balance(pkg, m.id, nu)
        reactions.append(Reaction(
            id=m.id, nu=nu, names=tuple(pkg.names),
            kinetics=_parse_kinetics(m.kinetics, nu, pkg, where) if m.kinetics else None,
            keq=(m.keq.A, m.keq.B, m.keq.C, m.keq.D) if m.keq else None,
            order=tuple(i for i in order if nu[i] != 0.0),
        ))
    return ReactionSet(reactions)


def _check_element_balance(pkg: PropertyPackage, rid: str, nu: np.ndarray) -> None:
    where = f"reaction {rid}"
    totals: dict[str, list[float]] = {}  # element -> [in reactants, in products]
    for i, v in enumerate(nu):
        if v == 0.0:
            continue
        atoms = pkg.consts.atomss[i]
        if not atoms:
            raise SimulationError(
                f"the database has no formula for {pkg.names[i]!r}, so the element balance "
                "cannot be checked", where=where, hint="use a compound with a known formula",
            )
        for el, count in atoms.items():
            side = totals.setdefault(el, [0.0, 0.0])
            side[0 if v < 0 else 1] += abs(v) * count
    off = {el: (a, b) for el, (a, b) in totals.items() if abs(a - b) > 1e-9 * max(a, b, 1.0)}
    if off:
        detail = "; ".join(f"{el}: {a:g} in reactants, {b:g} in products"
                           for el, (a, b) in sorted(off.items()))
        raise SimulationError(
            f"the reaction does not conserve elements ({detail})", where=where,
            hint="adjust the stoichiometric coefficients so every element appears equally on "
            "both sides (reactants negative, products positive)",
        )


def _parse_kinetics(k, nu: np.ndarray, pkg: PropertyPackage, where: str) -> Kinetics:
    if k.basis not in BASES:
        raise SimulationError(f"unknown rate basis {k.basis!r}", where=where,
                              hint=f"use one of: {', '.join(BASES)}")
    if k.reverse is not None and k.equilibrium_limited:
        raise SimulationError(
            "the reverse reaction is given twice: a reverse rate law and equilibrium_limited",
            where=where, hint="keep one: either the reverse rate law or the equilibrium limit",
        )
    forward = _parse_rate_law(k, np.maximum(-nu, 0.0), pkg, where, "forward")
    reverse = (_parse_rate_law(k.reverse, np.maximum(nu, 0.0), pkg, where, "reverse")
               if k.reverse is not None else None)
    return Kinetics(forward, k.basis, reverse, k.equilibrium_limited)


def _parse_rate_law(law, default_orders: np.ndarray, pkg: PropertyPackage, where: str,
                    label: str) -> RateLaw:
    try:
        Ea = to_si(law.Ea, "molar_energy")
        T_ref = to_si(law.T_ref, "temperature") if law.T_ref is not None else None
    except UnitError as exc:
        raise SimulationError(f"{label} rate law: {exc}", where=where) from exc
    if law.k0 < 0:
        raise SimulationError(f"{label} rate constant k0 must not be negative", where=where)
    orders = default_orders.copy()
    if law.orders is not None:
        orders[:] = 0.0
        for name, a in law.orders.items():
            if name not in pkg.names:
                raise SimulationError(
                    f"{label} reaction order given for {name!r}, which is not a component",
                    where=where, hint=f"use one of: {', '.join(pkg.names)}",
                )
            orders[pkg.names.index(name)] = float(a)
    return RateLaw(float(law.k0), Ea, tuple(float(a) for a in orders), T_ref)


# ---- thermochemistry -----------------------------------------------------------------------

def _cp_integral(pkg: PropertyPackage, i: int, T1: float, T2: float, over_T: bool = False):
    cp = pkg.props.HeatCapacityGases[i]
    try:
        value = (cp.T_dependent_property_integral_over_T(T1, T2) if over_T
                 else cp.T_dependent_property_integral(T1, T2))
    except (ValueError, TypeError, AttributeError):  # no correlation for this compound
        value = None
    if value is None:
        raise SimulationError(f"no ideal-gas heat capacity for {pkg.names[i]!r}",
                              hint="pick a compound with heat capacity data")
    return float(value)


def missing_gibbs_data(pkg: PropertyPackage, indices: Sequence[int]) -> list[str]:
    c = pkg.consts
    return [pkg.names[i] for i in indices if c.Hfgs[i] is None or c.Gfgs[i] is None]


def gibbs_over_rt(pkg: PropertyPackage, T: float, indices: Sequence[int],
                  why: str = "chemical equilibrium") -> np.ndarray:
    """Standard ideal-gas Gibbs energy of formation G°_i(T)/(R T) at 1 bar for ``indices``.

    G°(T) = H°(T) - T S°(T) with H°(T) = ΔHf°(298) + ∫Cp dT and the formation entropy
    S°(T) = (ΔHf° - ΔGf°)/298.15 + ∫Cp/T dT. Element reference terms cancel in any
    element-balanced reaction, so differences of these values give exact ΔG°(T).
    Other entries are NaN.
    """
    missing = missing_gibbs_data(pkg, indices)
    if missing:
        raise SimulationError(
            f"no ideal-gas formation data (ΔHf°, ΔGf°) for {', '.join(missing)}, needed for "
            f"{why}",
            hint="give the reaction a user equilibrium constant (keq), or leave these "
            "compounds out",
        )
    c = pkg.consts
    g = np.full(pkg.n, np.nan)
    for i in indices:
        H = c.Hfgs[i] + _cp_integral(pkg, i, T0, T)
        S = (c.Hfgs[i] - c.Gfgs[i]) / T0 + _cp_integral(pkg, i, T0, T, over_T=True)
        g[i] = (H - T * S) / (R * T)
    return g


def element_matrix(pkg: PropertyPackage, indices: Sequence[int]) -> tuple[list[str], np.ndarray]:
    """Elements and the (elements x species) atom-count matrix for ``indices``."""
    elements = sorted({el for i in indices for el in pkg.consts.atomss[i]})
    A = np.array([[pkg.consts.atomss[i].get(el, 0) for i in indices] for el in elements],
                 dtype=float)
    return elements, A


# ---- chemical equilibrium ------------------------------------------------------------------

def equilibrium_extents(n0: np.ndarray, nu: np.ndarray, ln_k: np.ndarray,
                        ln_p: float = 0.0, tol: float = 1e-10) -> np.ndarray:
    """Extents of reactions ``nu`` (reactions x species) at equilibrium from amounts ``n0``.

    Each reaction satisfies ``sum_i nu_i ln(x_i exp(ln_p)) = ln K``, i.e. activities
    a_i = x_i P/P° for an ideal gas (``ln_p = ln(P/P°)``) or a_i = x_i for an ideal solution
    (``ln_p = 0``). This is the minimum of the convex Gibbs energy over the extents with all
    amounts non-negative, found by damped Newton steps from a strictly interior point.
    Species that no extent can make positive (e.g. a product whose co-reactant is absent)
    are held at zero.
    """
    n0 = np.asarray(n0, dtype=float)
    scale = n0.sum()
    m = nu.shape[0]
    if scale <= 0.0 or m == 0:
        return np.zeros(m)
    b = n0 / scale
    active = np.flatnonzero(np.any(nu != 0.0, axis=0))
    bound = 1e3  # in units of the total feed; any physical extent is far smaller

    def lp(c, A_ub, b_ub, nvar):
        return linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=[(-bound, bound)] * nvar,
                       method="highs")

    # Species that must stay at zero, and the subspace of extents that keeps them there.
    pinned = []
    for i in active:
        if b[i] > 0.0:
            continue
        res = lp(-nu[:, i], -nu[:, active].T, b[active], m)
        if res.status != 0 or -res.fun <= 1e-12:
            pinned.append(i)
    Z = null_space(nu[:, pinned].T) if pinned else np.eye(m)
    if Z.shape[1] == 0:
        return np.zeros(m)
    nz = Z.T @ nu  # independent reactions that leave the pinned species alone
    nz[:, pinned] = 0.0
    lnk = Z.T @ ln_k
    free = np.array([i for i in active if i not in pinned], dtype=int)
    if free.size == 0:
        return np.zeros(m)
    k = nz.shape[0]

    # Strictly interior starting point: maximize the smallest free amount.
    A_ub = np.hstack([-nz[:, free].T, np.ones((len(free), 1))])
    res = linprog(np.r_[np.zeros(k), -1.0], A_ub=A_ub, b_ub=b[free],
                  bounds=[(-bound, bound)] * k + [(None, 1.0)], method="highs")
    if res.status != 0 or res.x[-1] <= 0.0:
        raise SimulationError("could not find a feasible starting point for the equilibrium",
                              hint="check that the feed contains the reactants")
    eta = res.x[:-1]
    present = b > 0.0
    present[free] = True
    dnu = nz.sum(axis=1)

    def amounts(e):
        return b + nz.T @ e

    def gibbs(e):
        n = amounts(e)
        nt = n[present].sum()
        p = n[present]
        return -lnk @ e + float(np.sum(p * (np.log(p / nt) + ln_p)))

    for _ in range(200):
        n = amounts(eta)
        nt = n[present].sum()
        grad = -lnk + nz[:, free] @ (np.log(n[free] / nt) + ln_p)
        if np.max(np.abs(grad)) < tol:
            break
        H = (nz[:, free] / n[free]) @ nz[:, free].T - np.outer(dnu, dnu) / nt
        step = np.linalg.lstsq(H, -grad, rcond=None)[0]
        dn = nz[:, free].T @ step
        shrinking = dn < 0.0
        alpha = 1.0
        if shrinking.any():
            alpha = min(1.0, 0.99 * float(np.min(-n[free][shrinking] / dn[shrinking])))
        g0, slope = gibbs(eta), float(grad @ step)
        while alpha > 1e-14 and gibbs(eta + alpha * step) > g0 + 1e-4 * alpha * slope:
            alpha *= 0.5
        if np.max(np.abs(alpha * step)) <= 1e-15 * (1.0 + np.max(np.abs(eta))):
            # No further progress in double precision (e.g. a species at ~1e-15 of the
            # feed); accept if the remaining error in ln Q is small.
            if np.max(np.abs(grad)) < 1e-6:
                break
            raise SimulationError("chemical equilibrium did not converge",
                                  hint="check the temperature and the equilibrium constants")
        eta = eta + alpha * step
    else:
        raise SimulationError("chemical equilibrium did not converge",
                              hint="check the temperature and the equilibrium constants")
    return (Z @ eta) * scale
