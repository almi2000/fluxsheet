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

### Plug-ins (`plugins.py`)

Another installed package can add unit operations, property methods and compounds without
changing Fluxsheet. It names them in entry points of its `pyproject.toml`:

| Entry-point group | Value | Registered in |
| --- | --- | --- |
| `fluxsheet.units` | a module of `UnitOperation` subclasses, a class, or a list | `units.UNIT_TYPES` |
| `fluxsheet.property_methods` | a `plugins.ActivityMethod` or `plugins.EquationOfState`, or a list | `thermo_package.METHODS` |
| `fluxsheet.components` | custom-compound definitions (as in a flowsheet's `custom_compounds`) | `thermo_package.PLUGIN_COMPOUNDS` |

A value may also be a function without arguments that returns one of these. Plug-ins load
when `fluxsheet.units` and `fluxsheet.thermo_package` are first imported, after Fluxsheet's
own entries.

- **No replacements.** A plug-in cannot replace a built-in or another plug-in. A unit
  type, method key or compound name that is already taken rejects that whole entry point. A
  plug-in compound cannot reuse a name the compound database knows. A flowsheet's own
  `custom_compounds` still override a plug-in compound of the same name.
- **Failures are reported.** An entry point that fails to import or register is skipped,
  and Fluxsheet keeps running. The reason appears in `fluxsheet plugins` (exit code 1),
  `GET /api/plugins` and the editor's palette. `FLUXSHEET_DISABLE_PLUGINS=1` loads no
  plug-ins.
- **Activity methods** provide only the excess Gibbs energy:
  `gibbs_excess(consts, indices)` returns a `thermo` `GibbsExcess` model. Fluxsheet adds the
  ideal gas (with carboxylic-acid dimers), vapor pressures, liquid enthalpies and the flash,
  as for NRTL. Henry components are not available with a plug-in method.
- **Equations of state** name a `thermo` mixture class such as `PRMIX`, the `kijs` and any
  extra keyword arguments. They get no volume translation. `thermo_package.is_eos` treats
  them like Peng-Robinson, for example for phase envelopes and the activity basis of
  reactive columns.
- **Energy balances.** Fluxsheet knows the energy terms of its own units by type. A
  plug-in unit declares its energy terms in `energy_results`, the result fields carrying
  heat or work into the unit (W), for example `("duty",)`. Without it the unit's energy
  balance is reported as not checkable.
- **The editor** lists plug-in units in their category of the palette with a "plug-in"
  tag, and plug-in methods with "(plug-in)". A unit without its own symbol in
  `web/src/icons/` is drawn as a plain box.

Plug-ins are imported into the Fluxsheet process (and the server) with its permissions;
only installed packages can register one. [`examples/plugin/`](../examples/plugin/) is a
complete one, with a unit, an activity method, a compound and two flowsheets that use them.
Dynamic models, costing and DWSIM import do not see plug-in units yet.

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
