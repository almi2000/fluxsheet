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

#### Pressure-driven dynamics (`dynamics/pressure.py`)

With `"pressure_driven": true` the flows between pressure nodes follow from their pressures
instead of the steady state. Nodes are vessels with a dynamic model (`PressureVessel`, used
for Flash and MixingTank in such a run), feeds and products. `Network` finds the *lines*
from a node's outlet through valves, pumps, compressors and heaters (one inlet, one outlet)
to the next node. A line is pressure-driven if it holds a valve with a flow coefficient
(`Cv` or `Kv` on the valve) or a pump with a curve; every other line keeps the flow-driven
rule, and a valve with a coefficient that sits on no line is reported in the warnings.

- `PressureVessel`: with gas in it at the start, its states are the component amounts and
  the internal energy U = H - P V; the contents fill the vessel, and `_UVFlash` finds T and
  P by Newton steps in (T, ln P) on T-P flashes (on one phase's properties while the
  contents stay one phase), converged far below the Jacobian's step. Without gas it keeps its
  pressure. Outlets on a pressure-driven line take the flow the line sets (`forced`); the
  liquid leaves at the gas pressure plus rho g level.
- `Network.update`, called in `_evaluate` after the controllers: the vessels' pressures,
  then the flow of each line: the valve equation (IEC 60534-2-1, liquid and gas with
  choking) at the pressures before and after the valve, with the pump curve and the other
  units' steady-state pressure changes in between, solved for the flow if the line has a
  pump. The line's flow goes to the vessel's outlet or scales the feed; the valve's outlet
  pressure is set before the sweep recalculates it.
- `Network.start`: openings not given are those that pass the steady-state flows; a line
  with a pump only gets the line resistance that makes the steady state balance. An
  undisturbed run therefore equals the flow-driven one.
- The Jacobian's sparsity adds, for each line, the node downstream and the controllers on
  the line to the source's dependencies (`Network.dependencies`).

Limits: flow runs downstream only (no reverse flow through a valve), one valve with a
coefficient and one pump curve per line, no flashing or cavitation in liquid valves, gamma
of a gas valve from the steady state, compressors keep their steady pressure rise, and
mixers, splitters and other units end a line (the flows there stay flow-driven).

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

Another installed package can add unit operations, property methods, compounds and dynamic
models without changing Fluxsheet. It names them in entry points of its `pyproject.toml`:

| Entry-point group | Value | Registered in |
| --- | --- | --- |
| `fluxsheet.units` | a module of `UnitOperation` subclasses, a class, or a list | `units.UNIT_TYPES` |
| `fluxsheet.property_methods` | a `plugins.ActivityMethod` or `plugins.EquationOfState`, or a list | `thermo_package.METHODS` |
| `fluxsheet.components` | custom-compound definitions (as in a flowsheet's `custom_compounds`) | `thermo_package.PLUGIN_COMPOUNDS` |
| `fluxsheet.dynamic_models` | a dict `{unit type: DynamicModel subclass}`, or a list of them | `dynamics.units.DYNAMIC_MODELS` |

A value may also be a function without arguments that returns one of these. Plug-ins load
when `fluxsheet.units`, `fluxsheet.thermo_package` and `fluxsheet.dynamics.units` are first
imported, after Fluxsheet's own entries; `plugins.report()` imports all three, so
`fluxsheet plugins` lists every group.

- **No replacements.** A plug-in cannot replace a built-in or another plug-in. A unit
  type, method key, compound name or unit type with a dynamic model that is already taken
  rejects that whole entry point. A plug-in compound cannot reuse a name the compound
  database knows. A flowsheet's own
  `custom_compounds` still override a plug-in compound of the same name.
- **Failures are reported.** An entry point that fails to import or register is skipped,
  and Fluxsheet keeps running. The reason appears in `fluxsheet plugins` (exit code 1),
  `GET /api/plugins` and the editor's palette. `FLUXSHEET_DISABLE_PLUGINS=1` loads no
  plug-ins.
- **Activity methods** provide only the excess Gibbs energy:
  `gibbs_excess(consts, indices)` returns a `thermo` `GibbsExcess` model. Fluxsheet adds the
  ideal gas (with carboxylic-acid dimers), vapor pressures, liquid enthalpies and the flash,
  as for NRTL. With `henry=True` light gases become Henry components as with NRTL
  (`henry.py`): `gibbs_excess` then gets the indices of the solvents only, possibly a single
  one, and Fluxsheet wraps the model it returns in `HenryDilutedModel`. Without it,
  `henry_components` is ignored with a warning and `gibbs_excess` gets every index. A model
  of another size than the indices given is an error.
- **Equations of state** name a `thermo` mixture class such as `PRMIX`, the `kijs` and any
  extra keyword arguments. They get no volume translation. `thermo_package.is_eos` treats
  them like Peng-Robinson, for example for phase envelopes and the activity basis of
  reactive columns.
- **Energy balances.** Fluxsheet knows the energy terms of its own units by type. A
  plug-in unit declares its energy terms in `energy_results`, the result fields carrying
  heat or work into the unit (W), for example `("duty",)`. Without it the unit's energy
  balance is reported as not checkable.
- **Dynamic models.** `fluxsheet.dynamic_models` maps a unit type, a plug-in's or a
  built-in one without a model (Pump, Mixer, ...), to a `DynamicModel` subclass (see "To give
  another unit a dynamic model" above). The unit type must be registered when the dynamic
  models load, so a plug-in's own unit type comes from its `fluxsheet.units` entry point.
  `BUILTIN_MODELS` reserves Fluxsheet's own types, the column's included; `PLUGIN_MODELS`
  records which plug-in gave which. The model is used in flow- and pressure-driven runs; in
  a pressure-driven one it ends a line as a non-vessel model (no pressure node). The editor
  gets its settings and results from `/api/unit-types` (`dynamic`) like the built-in ones.
- **Costing.** Fluxsheet costs its own units by type (`costing._HANDLERS`). A unit type
  without a handler may price itself with a class or static method `estimate_cost(c)`; `c` is
  a `costing.UnitCosting` with the unit's results, inlet and outlet streams, material,
  pressure, temperatures, U and parameters, and adds `equipment(...)` with a Turton
  `Correlation` (`VESSEL`, `EXCHANGER`, `PUMP`, ...) and a bare-module factor, `heat(...)`
  (utility and, by default, a shell-and-tube exchanger sized as for a Heater),
  `electricity(...)` and `note(...)`. The costing section's CEPCI, material, U, prices and
  hours apply as for built-in units. A failure in the method is reported as the unit's
  costing error. A unit without `estimate_cost` keeps the utility of its `duty` only
  ("equipment not costed"), one whose method adds nothing is listed as not costed. A class
  hook rather than another entry-point group: the cost belongs to the unit and needs nothing
  to register.
- **The editor** lists plug-in units in their category of the palette with a "plug-in"
  tag, and plug-in methods with "(plug-in)". A unit without its own symbol in
  `web/src/icons/` is drawn as a plain box.

Plug-ins are imported into the Fluxsheet process (and the server) with its permissions;
only installed packages can register one. [`examples/plugin/`](../examples/plugin/) is a
complete one, with a unit that prices itself and has a dynamic model, an activity method
with Henry components, a compound, and three flowsheets that use them. The DWSIM importer
does not know plug-in units.

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
