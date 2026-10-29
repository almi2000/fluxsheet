# Architecture

```
 web editor (React, React Flow)
        │  JSON over HTTP (/api/…)
        ▼
 server.py (FastAPI) ──► solver.py ──► units/ ──► stream.py ──► thermo_package.py ──► thermo / chemicals
                              │                                        ▲
                              └──────────── flowsheet.py (file format) │ compounds.py (lookup, search)
```

## Python core (`src/fluxsheet`)

| Module | Responsibility |
|---|---|
| `flowsheet.py` | The file format (pydantic models). Units, streams as `UNIT.port` connections, solver options. |
| `quantities.py` | Converts `"80 degC"`-style values to SI. Everything inside the core is SI. |
| `compounds.py` | Compound lookup and prefix search over the `chemicals` database. |
| `thermo_package.py` | Property methods, the automatic method suggestion, and `flash()` returning a `PhaseState`. |
| `fitting/` | Binary parameter regression: ThermoML parser and downloader, consistency tests, NRTL/UNIQUAC least-squares fits, and the parameter files in `data/binary_parameters/` that `thermo_package.py` reads. |
| `stream.py` | `Stream`: component molar flows plus a `PhaseState`. |
| `units/` | Unit operations. Each declares ports and parameters as data (`Param`, `Port`); the same schema validates input and generates the editor's forms. |
| `solver.py` | Sequential-modular solver: builds the graph, tears recycles (DFS back edges), orders units topologically, iterates tears with Wegstein acceleration, checks the mass balance. |
| `errors.py` | `SimulationError` with `where` and `hint`. |
| `server.py`, `cli.py` | HTTP API for the editor; command-line interface. |

### Adding a unit operation

1. Subclass `UnitOperation` in `units/library.py`: set `type`, `label`, `category`,
   `description`, `inlets`, `outlets`, `params`, and any `exclusive_groups` /
   `required_groups` (alternative specifications such as *temperature or duty*).
2. Implement `calculate(pkg, inlets) -> {outlet_port: Stream}`; put reported values
   (duty, power, …) in `self.results`.
3. Register it in `UNIT_TYPES`, add an icon in `web/src/icons.tsx`, and write a test.

The editor picks up the new unit automatically from `/api/unit-types`.

## Web editor (`web/`)

- `App.tsx`: canvas state (React Flow nodes and edges), conversion to and from the file
  format, run/save/open.
- `UnitNode.tsx`: node rendering with one handle per port.
- `Inspector.tsx`: forms generated from the unit schema, stream results, flowsheet settings
  (compound search, property method with suggestion).
- `ResultsPanel.tsx`: stream table.

The server serves `web/dist` when it exists, so a single `fluxsheet serve` runs everything.
