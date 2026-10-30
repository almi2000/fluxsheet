"""Command line: ``fluxsheet run file.json``, ``fluxsheet sensitivity file.json``,
``fluxsheet fit-binary`` and ``fluxsheet serve``."""

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
    # Binary parameter provenance, grouped by source so wide component lists stay readable.
    by_source: dict[str, list[str]] = {}
    for src in result.parameter_sources:
        pair = f"{src['pair']} ({src['method']})" if src.get("method") else src["pair"]
        by_source.setdefault(src["source"], []).append(pair)
    def rank(source: str) -> int:  # measured data first, then estimates, then nothing
        return {"ThermoML fit": 0, "ChemSep": 1}.get(source, 2 if source.startswith("UNIFAC") else 3)

    for source, pairs in sorted(by_source.items(), key=lambda kv: rank(kv[0])):
        shown = ", ".join(pairs) if len(pairs) <= 3 else f"{', '.join(pairs[:3])} and {len(pairs) - 3} more"
        print(f"  binary parameters, {source}: {shown}")
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


def _fit_binary(args) -> int:
    from pathlib import Path

    from .fitting import archive, pipeline, store
    from .fitting.regression import FitOptions
    from .thermo_package import _binary_splits

    def progress(msg):
        if not args.quiet:
            print(msg, file=sys.stderr)

    try:
        collected = pipeline.collect(args.component_a, args.component_b,
                                     cache=Path(args.cache) if args.cache else None,
                                     max_files=args.max_files, offline=args.offline,
                                     progress=progress)
    except (LookupError, OSError, archive.ArchiveError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    pair = collected.pair
    print(f"{pair.names[0]} (1) / {pair.names[1]} (2): {len(collected.datasets)} VLE data sets "
          f"in {collected.files} ThermoML files")
    for line in collected.skipped:
        print(f"  skipped {line}")
    options = FitOptions(model=args.model, alpha=args.alpha, max_pressure=args.max_pressure * 1e3,
                         require_consistency=not args.no_consistency)
    try:
        result = pipeline.fit(collected, options)
    except (ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print()
    for r in result.datasets:
        pt, at = r.point_test, r.area_test
        tests = (f"point {pt.value:.4f}" if pt.applicable else "point n/a") + (
            f", area {at.value:.1f}" if at.applicable else "")
        status = (f"used, {r.n_rejected} outliers" if r.used else f"not used: {r.reason}")
        print(f"  {r.dataset.doi} #{r.dataset.number} {r.dataset.label}, "
              f"{len(r.dataset.points)} pts [{tests}] {status}")
    entry = result.entry()
    print()
    print(f"{args.model.upper()}: " + ", ".join(f"{k} = {entry[k]:.6g}" for k in
                                                 ("a12", "a21", "b12", "b21"))
          + (f", alpha = {result.alpha}" if result.alpha is not None else "")
          + f"  (fitted: {', '.join(result.free)})")
    print(f"{entry['n_points']} points from {entry['n_datasets']} data sets, "
          f"{entry['n_rejected_points']} rejected; T {entry['T_range_K'][0]}-"
          f"{entry['T_range_K'][1]} K, P {entry['P_range_kPa'][0]}-{entry['P_range_kPa'][1]} kPa")
    rows = [("this fit", result.stats), *result.comparison.items()]
    keys = [k for k in ("AAD_T_K", "RMS_T_K", "AAD_P_pct", "RMS_P_pct", "AAD_y", "RMS_y")
            if k in result.stats]
    print("  " + "".ljust(16) + "".join(k.rjust(11) for k in keys))
    for label, stats in rows:
        print("  " + label.ljust(16) + "".join(
            (f"{stats[k]:.4g}" if stats.get(k) is not None else "-").rjust(11) for k in keys))
    if _binary_splits(*pair.cas):
        print("\nnote: UNIFAC predicts two liquid phases for this pair. The fit uses vapor-liquid "
              "data only; check its liquid-liquid predictions (mutual solubilities) before "
              "storing it.")
    if args.write:
        path = store.save_entry(args.model, entry,
                                Path(args.output) if args.output else None)
        print(f"\nwrote {path}")
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
    fit = sub.add_parser(
        "fit-binary", help="fit NRTL/UNIQUAC binary parameters to NIST ThermoML VLE data",
        description="Download the binary VLE data of two compounds from the NIST TRC ThermoML "
        "Archive (cached outside the repository), check their consistency, fit the model and "
        "report the deviations, compared with ChemSep's parameters on the same points.")
    fit.add_argument("component_a")
    fit.add_argument("component_b")
    fit.add_argument("--model", choices=["nrtl", "uniquac"], default="nrtl")
    fit.add_argument("--alpha", type=float, default=None,
                     help="NRTL non-randomness (default 0.3, or 0.2 if UNIFAC predicts a split)")
    fit.add_argument("--max-files", type=int, default=0,
                     help="read at most this many ThermoML files, those with most binary data "
                     "first (default 0: all)")
    fit.add_argument("--max-pressure", type=float, default=300.0, help="kPa (default 300)")
    fit.add_argument("--no-consistency", action="store_true",
                     help="also use data sets that fail the point consistency test")
    fit.add_argument("--cache", help="download cache (default: user cache directory)")
    fit.add_argument("--offline", action="store_true", help="use cached files only")
    fit.add_argument("--write", action="store_true",
                     help="store the result in the package's parameter file")
    fit.add_argument("--output", help="directory for --write instead of the package data")
    fit.add_argument("--quiet", action="store_true", help="no download progress")
    fit.set_defaults(func=_fit_binary)
    args = parser.parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
