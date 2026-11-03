# Architecture

```
 web editor (React, React Flow)
        │  JSON over HTTP (/api/…)
        ▼
 server.py (FastAPI) ──► solver.py ──► units/ ──► stream.py ──► thermo_package.py ──► thermo / chemicals
        │                  ▲   ▲                                       ▲
        │      specs.py ───┘   └─── calculator.py                     │ compounds.py (lookup, search)
        │      sensitivity.py                                          │
        └────► dynamics/ (simulator, unit models, controllers) ────────┘
                              flowsheet.py (file format)
```

## Python core (`src/fluxsheet`)

| Module | Responsibility |
|---|---|
| `flowsheet.py` | The file format (pydantic models). Units, streams as `UNIT.port` connections, solver options, and the optional sections below. |
| `quantities.py` | Converts `"80 degC"`-style values to SI. Everything inside the core is SI. |
| `compounds.py` | Compound lookup and prefix search over the `chemicals` database. |
| `thermo_package.py` | Property methods, the automatic method suggestion, `flash()` returning a `PhaseState`, single-phase properties without a flash. |
| `fitting/` | Binary parameter regression: ThermoML parser and downloader, consistency tests, NRTL/UNIQUAC least-squares fits, and the parameter files in `data/binary_parameters/` that `thermo_package.py` reads. |
| `property_sections.py` | Per-unit property method overrides: the package each unit uses, the T-P re-flash of streams crossing between methods, per-unit method suggestions. |
| `stream.py` | `Stream`: component molar flows plus a `PhaseState`. |
| `reactions.py` | Reaction sets: stoichiometry checks, kinetics, equilibrium constants. |
| `utilities.py` | Heating and cooling media; rates jackets and coils (U·A) together with a unit's energy balance. |
| `units/` | Unit operations, one module per family (`streams`, `mixers`, `heat`, `separators`, `columns`, `reactors`, `pressure*`, ...). Each declares ports and parameters as data (`Param`, `Port`); the same schema validates input and generates the editor's forms. |
| `columns/` | The equilibrium-stage solver (MESH equations by Newton, Naphtali-Sandholm) and the shortcut method. |
| `solver.py` | Sequential-modular solver: builds the graph, tears recycles (DFS back edges), orders units topologically, iterates tears with Wegstein or Broyden acceleration, checks the mass balance. `solve_flowsheet` is one solve with the calculator blocks. |
| `specs.py` | Design specs: an outer loop (secant, then Brent) around `solve_flowsheet`. Also `Variable` (a stream property or unit result) and `Handle` (a unit parameter), which sensitivity studies and calculators reuse. |
| `sensitivity.py` | Sensitivity studies over one or two parameters. |
| `calculator.py` | Calculator blocks: a safe arithmetic evaluator over `ast`, and the substitution loop when expressions read results. |
| `dynamics/` | Dynamic simulation, see below. |
| `limits.py` | Time and size limits of the public demo server (off locally). |
| `errors.py` | `SimulationError` with `where` and `hint`. |
| `server.py`, `cli.py` | HTTP API for the editor; command-line interface. |

### Dynamic simulation (`dynamics/`)

- `models.py`: the `"dynamics"` section of the file: end time, the units with a hold-up and
  their settings, controllers, events.
- `units.py`: dynamic models wrapping solved units. `Vessel` (Flash, MixingTank) and
  `StirredReactor` (CSTR) hold component amounts and enthalpy; `Lag` (Heater) holds its
  outlet enthalpy. Each gives its outlets from the state and its time derivatives. A
  "source" model's outlets depend on its state alone, so it cuts the flowsheet into pieces
  calculated in flow order. Single-phase contents skip the full flash with its stability
  test (`_ContentsFlash`) while they stay near a confirmed state.
- `controllers.py`: PID in percent of range; filtered measurement and integral are states,
  so the outputs never form an algebraic loop with the process.
- `simulator.py`: `Simulator` extends `Solver`: solves the steady state, initializes the
  models and controllers from it, and integrates with `scipy`'s BDF method. Each right-hand
  side applies events and controller outputs, evaluates the hold-ups, recalculates the
  algebraic units in order (cached by inlets and parameters, so finite-difference Jacobian
  columns only redo what a perturbed state reaches), and collects the derivatives. The
  Jacobian's sparsity comes from which hold-ups and controllers can reach which units.

To give another unit a dynamic model, subclass `DynamicModel` in `dynamics/units.py`
(`params`, `result_fields`, `initialize`, `outlets`, `derivatives`, `scale`) and register it
in `DYNAMIC_MODELS`. Start it from the unit's steady state, so a run without events stays put
(`tests/test_dynamics.py` checks this).

### Adding a unit operation

1. Subclass `UnitOperation` in a module under `units/` (or a new one, which
   `units/__init__.py` imports): set `type`, `label`, `category`, `description`, `inlets`,
   `outlets`, `params`, and any `exclusive_groups` / `required_groups` (alternative
   specifications such as *temperature or duty*).
2. Implement `calculate(pkg, inlets) -> {outlet_port: Stream}`; put reported values
   (duty, power, …) in `self.results` and declare them in `result_fields`.
3. Every subclass with a `type` is registered automatically. Add an equipment symbol in
   `web/src/icons/` and write a test.

The editor picks up the new unit automatically from `/api/unit-types`.

## Web editor (`web/`)

- `App.tsx`: canvas state (React Flow nodes and edges), conversion to and from the file
  format, run/save/open, undo/redo (`history.ts`), copy/paste (`clipboard.ts`), auto-layout
  (`layout.ts`).
- `UnitNode.tsx`, `StreamEdge.tsx`, `icons/`: nodes drawn as equipment, streams attached at
  their nozzles.
- `Inspector.tsx`: forms generated from the unit schema, stream results, flowsheet settings
  (compound search, property method with suggestion), with `ReactionsEditor.tsx`,
  `UtilitiesEditor.tsx`, `StudiesEditor.tsx` (design specs, sensitivity) and
  `DynamicsEditor.tsx` (hold-ups, controllers, events).
- `ResultsPanel.tsx`: stream table and tabs for phase diagrams (`Diagrams.tsx`),
  sensitivity studies (`SensitivityView.tsx`) and dynamic runs (`DynamicsView.tsx`).

The server serves `web/dist` when it exists, so a single `fluxsheet serve` runs everything.
