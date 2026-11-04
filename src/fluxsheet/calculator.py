"""Calculator blocks: set unit parameters from formulas over flowsheet values.

A calculator names some variables, each a stream property or unit result (as in design specs)
or the current value of a unit parameter, and assigns expressions over them to parameters,
e.g. "make-up flow = methanol flow in the purge". Everything is in SI units (K, Pa, mol/s, W).

Expressions are plain arithmetic: numbers, the variable names, ``+ - * / ** %``, comparisons,
``and``/``or``/``not``, ``a if cond else b``, the constants ``pi`` and ``e`` and the functions
in ``FUNCTIONS``. They are parsed with ``ast`` and evaluated by walking the tree; nothing is
passed to ``eval``, so a flowsheet file cannot run code.

Assignments that only read parameters are applied once, before the flowsheet is solved. When
an expression reads a result, the flowsheet is solved, the expressions evaluated and the
parameters updated, until they stop changing (successive substitution, at most
``MAX_PASSES`` solves).
"""

from __future__ import annotations

import ast
import math
import operator
from dataclasses import dataclass
from typing import Any

from .errors import SimulationError
from .flowsheet import FlowsheetModel
from .quantities import UnitError, to_si
from .solver import Solver, SolveResult
from .specs import Variable, _check_component, _unit_class, set_param
from .study_models import Assignment, CalculatorModel, ParamRef

MAX_PASSES = 30
TOLERANCE = 1e-7  # relative change of every assigned value between passes

FUNCTIONS = {
    "abs": abs, "min": min, "max": max, "sqrt": math.sqrt, "exp": math.exp, "log": math.log,
    "log10": math.log10, "sin": math.sin, "cos": math.cos, "tan": math.tan,
    "floor": math.floor, "ceil": math.ceil,
}
CONSTANTS = {"pi": math.pi, "e": math.e}
_BINARY = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
           ast.Div: operator.truediv, ast.Pow: operator.pow, ast.Mod: operator.mod}
_COMPARE = {ast.Lt: operator.lt, ast.LtE: operator.le, ast.Gt: operator.gt,
            ast.GtE: operator.ge, ast.Eq: operator.eq, ast.NotEq: operator.ne}
_UNARY = {ast.USub: operator.neg, ast.UAdd: operator.pos, ast.Not: operator.not_}


# ---- expressions -------------------------------------------------------------------------

class Expression:
    """A parsed, checked arithmetic expression; see the module docstring."""

    def __init__(self, text: str, names: set[str], where: str):
        self.text, self.where = text, where
        try:
            self.tree = ast.parse(text.strip(), mode="eval").body
        except SyntaxError as exc:
            raise SimulationError(f"cannot read the expression {text!r}: {exc.msg}",
                                  where=where) from None
        self.names = set()
        self._check(self.tree, names)

    def _check(self, node: ast.AST, names: set[str]) -> None:
        if isinstance(node, ast.Constant):
            if isinstance(node.value, bool) or not isinstance(node.value, (int, float)):
                raise self._error(f"{node.value!r} is not a number")
        elif isinstance(node, ast.Name):
            if node.id not in names and node.id not in CONSTANTS:
                raise self._error(f"unknown name {node.id!r}",
                                  f"define it under 'variables'; known: "
                                  f"{', '.join(sorted(names | set(CONSTANTS)))}")
            self.names.add(node.id)
        elif isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name) or node.func.id not in FUNCTIONS:
                raise self._error("only these functions can be called: "
                                  + ", ".join(FUNCTIONS))
            if node.keywords:
                raise self._error("functions take no keyword arguments")
            for arg in node.args:
                self._check(arg, names)
        elif isinstance(node, ast.BinOp) and type(node.op) in _BINARY:
            self._check(node.left, names)
            self._check(node.right, names)
        elif isinstance(node, ast.UnaryOp) and type(node.op) in _UNARY:
            self._check(node.operand, names)
        elif isinstance(node, ast.Compare) and all(type(o) in _COMPARE for o in node.ops):
            for sub in (node.left, *node.comparators):
                self._check(sub, names)
        elif isinstance(node, ast.BoolOp):
            for sub in node.values:
                self._check(sub, names)
        elif isinstance(node, ast.IfExp):
            for sub in (node.test, node.body, node.orelse):
                self._check(sub, names)
        else:
            raise self._error(f"{type(node).__name__} is not allowed in an expression",
                              "use numbers, variable names, arithmetic and the listed functions")

    def evaluate(self, values: dict[str, float]) -> float:
        try:
            value = self._eval(self.tree, values)
        except (ArithmeticError, ValueError) as exc:
            raise self._error(f"cannot evaluate {self.text!r}: {exc}") from None
        if isinstance(value, bool) or not math.isfinite(value):
            raise self._error(f"{self.text!r} does not give a finite number")
        return float(value)

    def _eval(self, node: ast.AST, v: dict[str, float]) -> Any:
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            return v[node.id] if node.id in v else CONSTANTS[node.id]
        if isinstance(node, ast.Call):
            return FUNCTIONS[node.func.id](*(self._eval(a, v) for a in node.args))
        if isinstance(node, ast.BinOp):
            return _BINARY[type(node.op)](self._eval(node.left, v), self._eval(node.right, v))
        if isinstance(node, ast.UnaryOp):
            return _UNARY[type(node.op)](self._eval(node.operand, v))
        if isinstance(node, ast.Compare):
            left = self._eval(node.left, v)
            for op, comp in zip(node.ops, node.comparators):
                right = self._eval(comp, v)
                if not _COMPARE[type(op)](left, right):
                    return False
                left = right
            return True
        if isinstance(node, ast.BoolOp):
            values = (self._eval(sub, v) for sub in node.values)
            return all(values) if isinstance(node.op, ast.And) else any(values)
        if isinstance(node, ast.IfExp):
            branch = node.body if self._eval(node.test, v) else node.orelse
            return self._eval(branch, v)
        raise AssertionError("unchecked node")

    def _error(self, message: str, hint: str | None = None) -> SimulationError:
        return SimulationError(message, where=self.where, hint=hint)


# ---- variables and targets ---------------------------------------------------------------

def _param_quantity(model: FlowsheetModel, unit: str, param: str, component: str | None,
                    where: str) -> str | None:
    """The quantity of a continuous unit parameter (per component for flows and fractions)."""
    unit_cls = _unit_class(unit, model, where)
    params = {p.name: p for p in unit_cls.params}
    if param not in params:
        raise SimulationError(f"{unit} ({unit_cls.label}) has no parameter {param!r}",
                              where=where, hint=f"use one of: {', '.join(params)}")
    q = params[param].quantity
    if q in ("composition", "component_fractions"):
        _check_component(component, model, where)
        return "molar_flow" if q == "composition" else "fraction"
    if q in ("integer", "choice", "component", "string_list", "stream_stages", "utility",
             "records"):
        raise SimulationError(f"{param} of {unit} is not a number a calculator can set or read",
                              where=where)
    return q


def _read_param(model: FlowsheetModel, ref: ParamRef | Assignment, where: str) -> float | None:
    q = _param_quantity(model, ref.unit, ref.param, ref.component, where)
    raw = next(u for u in model.units if u.id == ref.unit).params.get(ref.param)
    if ref.component is not None:
        raw = (raw or {}).get(ref.component)
    if raw is None:
        return None
    try:
        return to_si(raw, q) if q else float(raw)
    except (UnitError, ValueError) as exc:
        raise SimulationError(f"{ref.param} of {ref.unit}: {exc}", where=where) from exc


@dataclass
class _Calculator:
    model: CalculatorModel
    results: dict[str, Variable]  # variables read from the solution
    params: dict[str, ParamRef]  # variables read from the parameters
    expressions: list[tuple[Assignment, Expression]]

    @classmethod
    def build(cls, calc: CalculatorModel, model: FlowsheetModel) -> _Calculator:
        results, params = {}, {}
        for name, ref in calc.variables.items():
            if not name.isidentifier() or name in FUNCTIONS or name in CONSTANTS:
                raise SimulationError(f"{name!r} cannot be a variable name", where=calc.id,
                                      hint="use letters, digits and _, not a function name")
            if isinstance(ref, ParamRef):
                _read_param(model, ref, calc.id)  # checks the reference
                params[name] = ref
            else:
                results[name] = Variable.resolve(ref, model, calc.id)
        expressions = []
        for assign in calc.assign:
            _param_quantity(model, assign.unit, assign.param, assign.component, calc.id)
            expressions.append((assign, Expression(assign.expression, set(calc.variables),
                                                   calc.id)))
        return cls(calc, results, params, expressions)

    @property
    def reads_results(self) -> bool:
        return any(e.names & set(self.results) for _, e in self.expressions)

    def apply(self, model: FlowsheetModel, result: SolveResult | None) -> list[dict]:
        """Evaluate every assignment against ``model`` (and ``result``) and set it in ``model``;
        returns what was set."""
        where = self.model.id
        values: dict[str, float] = {}
        for name, ref in self.params.items():
            value = _read_param(model, ref, where)
            if value is None:
                raise SimulationError(f"variable {name}: {ref.param} of {ref.unit} is not set",
                                      where=where)
            values[name] = value
        if result is not None:
            for name, var in self.results.items():
                try:
                    values[name] = var.read(result)
                except SimulationError as exc:
                    raise exc.located(where)
        done = []
        for assign, expr in self.expressions:
            if result is None and expr.names & set(self.results):
                continue
            value = expr.evaluate(values)
            set_param(model, assign.unit, assign.param, assign.component, value)
            done.append({"unit": assign.unit, "param": assign.param,
                         "component": assign.component, "value": value})
        return done


# ---- solving -----------------------------------------------------------------------------

def solve_with_calculators(model: FlowsheetModel) -> SolveResult:
    """Solve ``model`` with its active calculator blocks applied; see the module docstring."""
    work = model.model_copy(deep=True, update={"calculators": None})
    try:
        calcs = [_Calculator.build(c, work) for c in model.calculators or () if c.active]
        for calc in calcs:  # assignments from parameters alone hold before the first solve
            calc.apply(work, None)
        loop = [c for c in calcs if c.reads_results]
        result = Solver(work).solve()
        passes = 1
        while loop and result.ok:
            if passes >= MAX_PASSES:
                raise SimulationError(
                    f"calculator blocks did not settle in {MAX_PASSES} flowsheet solves",
                    where=", ".join(c.model.id for c in loop),
                    hint="an assignment feeds back on itself too strongly; check the formulas",
                )
            before = [_assigned(work, c) for c in loop]
            for calc in loop:
                calc.apply(work, result)
            if all(_settled(b, _assigned(work, c)) for b, c in zip(before, loop)):
                break
            result = Solver(work).solve()
            passes += 1
    except SimulationError as exc:
        return SolveResult(ok=False, error=exc.to_dict())
    if result.ok:
        report = []
        for calc in calcs:
            values = calc.apply(work.model_copy(deep=True), result)
            report.append({"id": calc.model.id, "assigned": values, "passes": passes})
        result.calculators = report
    return result


def _assigned(model: FlowsheetModel, calc: _Calculator) -> list[float | None]:
    return [_read_param(model, a, calc.model.id) for a, _ in calc.expressions]


def _settled(before: list[float | None], after: list[float | None]) -> bool:
    return all(a is not None and b is not None and abs(a - b) <= TOLERANCE * max(abs(a), 1e-12)
               for a, b in zip(before, after))
