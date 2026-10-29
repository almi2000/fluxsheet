"""Command line: ``fluxsheet run file.json``, ``fluxsheet sensitivity file.json`` and
``fluxsheet serve``."""

from __future__ import annotations

import argparse
import json
import sys

from .flowsheet import FlowsheetModel
from .sensitivity import run_sensitivity, to_csv
from .solver import solve
from .units import UNIT_TYPES


def _fmt_table(result) -> str:
    rows = [("Stream", "T [°C]", "P [bar]", "Vap. frac.", "Flow [kmol/h]", "Flow [kg/h]")]
    for sid, s in result.streams.items():
        rows.append((
            sid,
            f"{s['T'] - 273.15:.2f}",
            f"{s['P'] / 1e5:.4g}",
            f"{s['vapor_fraction']:.4f}",
            f"{s['molar_flow'] * 3.6:.4g}",
            f"{s['mass_flow'] * 3600:.4g}",
        ))
    widths = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    lines = ["  ".join(c.rjust(w) if i else c.ljust(w) for i, (c, w) in enumerate(zip(r, widths)))
             for r in rows]
    lines.insert(1, "-" * len(lines[0]))
    return "\n".join(lines)


def _run(args) -> int:
    result = solve(FlowsheetModel.load(args.file))
    if args.json:
        print(json.dumps(result.to_dict(), indent=2))
        return 0 if result.ok else 1
    if not result.ok:
        err = result.error
        print(f"error{' in ' + err['where'] if err['where'] else ''}: {err['message']}",
              file=sys.stderr)
        if err["hint"]:
            print(f"  hint: {err['hint']}", file=sys.stderr)
        return 1
    print(f"Property method: {result.property_method}")
    if result.property_method_reason:
        print(f"  ({result.property_method_reason})")
    if result.tear_streams:
        print(f"Recycle: tear {', '.join(result.tear_streams)}, "
              f"converged in {result.iterations} iterations")
    print()
    print(_fmt_table(result))
    # Energy results each unit declares (duty, power, condenser duty, ...), in kW.
    lines = []
    for uid, u in result.units.items():
        cls = UNIT_TYPES.get(u.get("type"))
        for f in cls.result_fields if cls else ():
            if f.quantity == "power" and isinstance(u.get(f.name), (int, float)):
                lines.append(f"{uid}: {f.label.lower()} {u[f.name] / 1e3:.4g} kW")
    if lines:
        print()
        print("\n".join(lines))
    for spec in result.design_specs:
        print(f"design spec {spec['id']}: {spec['target']} = {spec['achieved']:.6g} with "
              f"{spec['vary']} = {spec['vary_value']:.6g} ({spec['iterations']} solves)")
    for w in result.warnings:
        print(f"warning: {w}")
    return 0


def _sensitivity(args) -> int:
    table = run_sensitivity(FlowsheetModel.load(args.file), args.study)
    if not table["ok"]:
        err = table["error"]
        print(f"error{' in ' + err['where'] if err['where'] else ''}: {err['message']}",
              file=sys.stderr)
        if err["hint"]:
            print(f"  hint: {err['hint']}", file=sys.stderr)
        return 1
    text = to_csv(table)
    if args.csv:
        with open(args.csv, "w", encoding="utf-8", newline="") as f:
            f.write(text)
        print(f"{len(table['rows'])} points written to {args.csv}")
    else:
        print(text, end="")
    for row in table["rows"]:
        if row["error"]:
            print(f"warning: point {row['inputs']} failed: {row['error']['message']}",
                  file=sys.stderr)
    return 0


def _serve(args) -> int:
    import uvicorn

    uvicorn.run("fluxsheet.server:app", host=args.host, port=args.port, reload=args.reload)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="fluxsheet", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="solve a flowsheet file and print the results")
    run.add_argument("file")
    run.add_argument("--json", action="store_true", help="print full results as JSON")
    run.set_defaults(func=_run)
    sens = sub.add_parser("sensitivity", help="run a sensitivity study and print it as CSV")
    sens.add_argument("file")
    sens.add_argument("--study", help="id of the study (default: the first)")
    sens.add_argument("--csv", help="write the table to this file instead of printing it")
    sens.set_defaults(func=_sensitivity)
    serve = sub.add_parser("serve", help="start the API server for the web editor")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8000)
    serve.add_argument("--reload", action="store_true")
    serve.set_defaults(func=_serve)
    args = parser.parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
