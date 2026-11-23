# Changelog

All notable changes to Fluxsheet. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/) (before 1.0, minor versions may change the file
format or API).

## [Unreleased]

### Added

- **Pressure-driven dynamics** (`"pressure_driven": true` under `dynamics`): flows between
  vessels, feeds and products follow from their pressures. Valves get a flow coefficient
  (`Cv` or `Kv`), an opening that controllers and events move, a linear, equal-percentage or
  quick-opening characteristic and the IEC 60534 liquid and gas equations with choked gas
  flow. Vessels with gas in them carry their pressure (the contents fill the volume), so it
  builds up behind a closed valve and relieves through an open one; the liquid outflow adds
  the static head of the level. Pumps with a curve deliver against the line. Valves start
  at the opening that passes the steady-state flow, so an undisturbed run stays put. The
  flow-driven mode stays the default. Example: `flash_drum_pressure_control.json`.
- **Cooled tubular reactor:** the plug flow reactor exchanges heat along its length with a
  utility (cooling water, steam, thermal oil) or with a process stream on new medium ports
  (`medium: "stream"`, ports `coolant_in` and `coolant_out`), co- or countercurrent. Give U
  with the length and diameter; the area per length defaults to the tube circumference
  (`area_per_length`). Countercurrent is solved as the boundary-value problem it is, by
  shooting on the duty with a bracketing root finder; a fluid utility without a flow takes
  the flow that brings it to its outlet temperature. Results: both temperature profiles,
  the heat through the wall, the medium's outlet state, and the hot-spot temperature and
  position. Example: `examples/phthalic_anhydride_tube.json`.

### Changed

- The plug flow reactor takes its heat capacity by a central difference, which moves
  temperature profiles by up to about a millikelvin. A PFR with a utility and U but neither
  temperature nor duty now rates the wall along the tube (before, nothing was exchanged)
  and needs the length and diameter for it.

### Fixed

- **Three-phase flash with vapor:** with a vapor and two liquids, the heavy-liquid fraction
  was the heavy liquid's share of the liquid instead of the mixture, so decanters with a gas
  phase put too much into their heavy outlet. Found by the DWSIM comparison below.
- **DWSIM import of real DWSIM files:** stream states were not read from files DWSIM
  writes, which store a phase's properties and compounds twice (a type name first, the data
  after). Such files lost every feed. Three-phase vessels now assign their liquids by the
  stored densities. On a cubic equation of state, which finds no second liquid in Fluxsheet,
  they split with UNIFAC, with a warning. Gibbs reactors minimize on the property method's
  fugacities, as DWSIM's do. `fluxsheet.importers.stored_streams` returns the stream states
  DWSIM saved. Three of DWSIM's own sample simulations (Cavett's problem, a three-phase
  separator, Gibbs and equilibrium reactors) are identified by SHA-256 in
  `tests/data/dwsim/real-samples.json`, but not redistributed because they are GPL-licensed.
  `scripts/fetch_dwsim_samples.py` downloads them, and the tests then compare every stream
  with DWSIM's result: within 1 K and 2-8 % in flow for Cavett's problem, 1.5 % for the
  separator and 5e-5 for the Gibbs reactor.
- **Reactive columns on a boilup ratio:** an equilibrium reactive column specified by its
  boilup ratio now converges by continuation from a column at moderate reflux, as one
  specified by its reflux ratio already did. A direct solve that makes no progress is
  abandoned after 20 iterations, which halves the solve time at high reflux (methyl acetate
  column at reflux 3-5: about 8 s instead of 14 s).

## [0.3.0] - 2026-11-23

The first published release; the prepared 0.2.0 was never published and is included here.
New since 0.2.0: dynamic simulation (vessels, reactors, heaters and the rigorous column,
with PID controllers and events), batch reactors, columns with side draws, pumparounds and
kinetic or equilibrium reactive stages, a Gibbs reactor on fugacities, equipment costing and
optimization, pump and compressor curves, logic blocks, Henry components and user-defined
compounds, Excel and HTML reports, DWSIM import and plug-ins for units, property methods and
compounds. Several results change through corrections listed under Fixed: liquid densities
of the cubic equations, heats of vaporization of the activity methods, carboxylic-acid vapor
and two-phase pipe flow.

### Added

- **Plug-ins:** other installed packages can add unit operations, property methods
  (activity models or cubic equations of state) and compounds through the entry points
  `fluxsheet.units`, `fluxsheet.property_methods` and `fluxsheet.components`. Plug-ins cannot
  replace anything built in. One that fails to load is skipped and reported by
  `fluxsheet plugins`, `GET /api/plugins` and the editor palette, which also marks plug-in
  units. Plug-in units declare their energy terms for the balance checks. A complete
  example is in `examples/plugin/`.

- **Feature reference campaign:** `scripts/validate_features.py` compares the Gibbs
  reactor on fugacities with the Larson and Dodge (1923) and Larson (1924) ammonia
  equilibria and the equilibrium methyl acetate column with the measured profile in
  Mekala (2019), against limits registered in `docs/validation-campaign-features.json`.
  Ammonia passes at 10-300 atm and fails at 600-1000 atm; the equilibrium column fails
  against the kinetically limited pilot column. `docs/model-validation.md` records the
  results and the verification of the dynamic column and batch costing.

- **Gibbs reactor on fugacities:** `fugacities: property_method` minimizes the Gibbs
  energy on fugacities from the property method (equation of state or activity model)
  with the vapor and liquid phases the mixture forms, instead of an ideal gas. Ammonia at
  300 atm comes within 0.1-1.5 points of Larson's (1924) data (the ideal gas is 4-5
  points low); condensing methanol and water drive the synthesis almost to completion.
  Example `ammonia_equilibrium_fugacities.json`.
- **Batch reactor costing:** the `BatchReactor` gets Turton costs as a jacketed agitated
  vessel (Table A.1, the CSTR's bare-module factor, scaled for material and pressure like
  the other vessels) holding its largest contents at 80 % fill, once per vessel the cycle
  needs: a vessel is occupied for charge, reaction and turnaround (the unit's down time,
  else the new costing setting `batch_turnaround`, default 1 h). Heating and cooling are
  priced apart per batch (the unit now reports `heating_per_batch` and
  `cooling_per_batch`) and per year from the batches per year (operating hours over the
  cycle time). Results `batches_per_year`, `vessels_costed`, `utility_cost_per_batch` and
  `cost_per_batch`; the batch reactor counts in `tac`, so studies, design specs and
  optimizations of batch flowsheets see it. Editor field for the turnaround.
- **DWSIM import:** `fluxsheet import-dwsim FILE [-o OUT.json]`, `POST /api/import/dwsim`
  and `fluxsheet.importers.import_dwsim` turn a DWSIM simulation (`.dwxml`, zipped
  `.dwxmz`) into a flowsheet: compounds, property method, feeds, mixers, splitters,
  heaters and coolers, pumps, compressors, expanders, valves, separators, heat exchangers,
  shortcut columns, compound separators, conversion and Gibbs reactors with their specs,
  the connections and the layout. Unsupported parts (energy streams, logical blocks,
  other units and property packages, DWSIM's interaction parameters) come back as
  warnings; an unsupported unit is left out and its streams become feeds and products at
  their stored state. Safe XML parsing (defusedxml), size and zip-bomb limits.
- **Dynamic distillation column:** the rigorous column under `dynamics.units`, started
  from its steady MESH solution. Trays hold liquid over the weir (Francis), reflux drum and
  sump hold minutes of liquid; reflux, distillate and bottoms flows and the reboiler duty
  are settings controllers move. The condenser condenses all vapor (the default for a
  total condenser) or runs on its duty and vents what it leaves. Results include drum and
  sump levels, a control-stage temperature and every stage's temperature. At rest the
  model ends on the steady column with the same flows. Example
  `benzene_toluene_column_control.json`. Side draws, pumparounds, stage duties, reactive
  stages and Murphree efficiencies are not modelled dynamically yet.
- **Result overlays in the editor:** a toolbar choice of what stream labels show after a
  run (id, temperature, T and P, plus molar flow, plus the main component); remembered per
  browser. Unit notes and unit balance checks that do not close from the last run appear
  as a warning marker (with a count) on the unit.
- **Optimization around the flowsheet:** objective and constraints as calculator
  expressions, bounded unit parameters as decision variables, SLSQP on scaled variables
  with a COBYLA fallback when trials fail. CLI `fluxsheet optimize`, `POST /api/optimize`,
  editor settings and results tab. Example `ethyl_acetate_optimization.json`.
- **Equipment costing:** Turton purchased and bare-module costs (CEPCI-scaled) sized from
  the solution for columns, drums, exchangers, pumps, compressors and stirred reactors;
  utilities priced per GJ; total annualized cost (`tac`) readable by specs, sensitivity
  studies and optimizations.
- **Column side streams:** liquid or vapor side draws, pumparounds (duty or return
  temperature), side heaters and coolers, and kinetic reactive stages on the rigorous
  column. Example `btx_side_draw_pumparound.json`.
- **Chemical equilibrium on distillation stages.** Selected stages solve liquid reaction
  equilibrium with a temperature-dependent `keq` on a mole-fraction or activity basis.
  Reaction extents enter the stage material and energy balances. Example
  `methyl_acetate_equilibrium_column.json`.
- **Performance curves:** pump head and efficiency curves with affinity laws and NPSH
  available (warning below NPSH required); compressor curves (polytropic or isentropic)
  with fan laws and surge/stonewall warnings.
- **Logic blocks:** Mult, Dupl and Selector.
- **Henry components:** light and supercritical gases in NRTL, UNIQUAC, UNIFAC and Wilson
  dissolve by Henry's law instead of an extrapolated vapor pressure; the components and the
  source of each constant are reported.
- **User-defined compounds:** a flowsheet can define its own compounds (formula, SMILES or
  Joback groups, known constants); missing properties are estimated by Joback and marked as
  estimates. Editor form in the components section.
- **Batch reactor:** one batch cycle integrated with stop criteria (time, conversion,
  composition, temperature) at constant pressure (optionally vented at the bubble point)
  or constant volume; isothermal, adiabatic, duty, temperature program or utility
  jacket. The product leaves as the cycle-average stream, with time profiles in the
  results.
- **Reports with runs:** the report also holds the costing summary and the editor's latest
  sensitivity, optimization and (in Excel) dynamic runs; `POST /api/report` takes
  `{model, runs}`, the CLI `--with-runs`.
- **Editor:** copy and paste carries user-defined compounds, and renaming one updates every
  reference; autosave no longer lets one browser tab silently overwrite another's newer
  changes; "Use optimum" is disabled after the flowsheet changed; a warning for connected
  column side outlets without a side draw.
- **Reports:** Excel workbook (summary, stream table, unit results, profiles, design specs)
  and a printable HTML report, from the results panel, `fluxsheet report` or
  `POST /api/report`. Excel needs the `report` extra (`openpyxl`), included in `server`.
- **Editor autosave:** the open flowsheet is kept in the browser and restored after a
  reload; the Open menu lists recently opened and saved flowsheets.

- **Dynamic simulation:** from the steady state, flash drums and mixing tanks (level, hold-up),
  CSTRs (kinetics and heat of reaction) and heaters (thermal lag) are integrated over time
  with a stiff BDF method and a sparse Jacobian. PID controllers in percent of range with
  sensor lag and anti-windup; steps and ramps on any unit parameter or setpoint. Editor:
  Dynamics section in the flowsheet settings, Dynamics tab with charts and CSV export; CLI
  `fluxsheet dynamic`; API `POST /api/dynamic`. Example `esterification_control.json`.
- **Calculator blocks:** set unit parameters from arithmetic expressions over stream
  properties, unit results and other parameters; parsed safely, iterated with the solve
  when they read results.
- **Broyden recycle convergence** (`"solver": {"acceleration": "broyden"}`) next to Wegstein.

- **Equipment symbols:** units are drawn as shaded equipment (drums with liquid levels, a
  column with condenser, reflux drum and reboiler, pump, stirred tank, shell and tube
  exchanger, ...) with streams attached at their nozzles.
- **Several feeds per unit:** most units take any number of streams on their inlet, mixed
  adiabatically before they enter. A distillation column keeps its feeds apart and takes a
  feed stage per stream.
- Dropping a stream on a unit's body connects it to the nearest free port.
- **Utilities and jacket heat transfer:** cooling water, steam (IAPWS-95), thermal oil and
  other media defined per flowsheet. Heaters, flash drums, decanters and reactors report
  their utility use and the UA or area a duty needs, or take U and the area of a jacket or
  coil and find the duty together with their energy balance. Heated and cooled units show
  a red or blue jacket.
- **More unit operations:** thin-film evaporator (plug flow along the heated wall,
  differential evaporation), mixing tank, turbine/expander, multistage compressor with
  intercoolers and condensate drain, pipe (friction, height change, heat loss) and
  liquid-liquid extraction column. Example `glycol_evaporator.json`.
- **Product rates relative to the feed:** a distillation column takes a distillate or
  bottoms to feed ratio, which follows the feed inside a recycle loop.

### Fixed

- **Chemical equilibrium near complete conversion:** the equilibrium solver no longer
  stops short when a species falls to around 1e-9 of the feed, where the Gibbs energy
  changes by less than its rounding; it accepts Newton steps that halve the gradient.
- **Convergence of reactive distillation.** If the direct Newton solve fails, columns with
  chemical equilibrium stages can continue from a converged solution at lower reflux to
  the specified reflux ratio with adaptive steps.
- **Controllers on bounded parameters:** a PI on a pump's speed no longer fails with a
  division by zero when it saturates at an output range starting at 0. Controller output
  ranges are narrowed to the manipulated parameter's own bounds (with a warning when that
  matters), and controller and event values are clamped to them.
- A compressed supercritical gas (hydrogen at 70 bar or more, nitrogen at 340 bar) is
  identified as vapor, not liquid, while superheated vapors stay vapor.
- A bubble or dew point that thermo converges to the other boundary is refused with an
  error instead of returned as a plausible wrong temperature.
- A distillation column whose feed carries dissolved hydrogen or methane starts without a
  bubble point of the feed; when it fails, the hint names a total condenser or a distillate
  too small for these gases as the likely cause.

- A pressure-enthalpy flash that thermo cannot converge (e.g. an almost pure component in a
  mixture model) is solved by temperature with T-P flashes instead of failing.
- A long flowsheet name no longer pushes the editor wider than the window.
- **Heat of vaporization of the activity-coefficient methods** (NRTL, UNIQUAC, UNIFAC,
  Wilson, ideal): the liquid enthalpy now uses the compounds' heat of vaporization
  correlations. It came from the vapor pressure curve with an ideal vapor
  (Clausius-Clapeyron), which is 4 to 6 % too high at 1 atm, 11 to 20 % at 5 bar and 20 to
  40 % at 10 bar. Condenser, reboiler and vaporizer duties of these methods change
  accordingly; phase equilibria do not. Dissolved gases keep their heat of solution from the
  Henry constant.
- **Liquid density of Peng-Robinson and Soave-Redlich-Kwong:** the volumes they report
  carry a Peneloux volume translation that matches each compound's liquid volume at 25 °C
  (or Tr = 0.7 for light compounds). The plain equations were 15 % (water, Peng-Robinson)
  and 8 to 24 % (Soave) off, and with them pump power, velocities, residence times and
  liquid concentrations in rate laws. Phase equilibria and enthalpies are unchanged.
- **Pipe:** gas and flashing flow include the pressure drop of acceleration and the kinetic
  energy; the outlet Mach number is reported and a flow that would exceed the speed of
  sound stops with a choking error instead of a too small pressure drop. Liquid results
  are unchanged.
- **Multistage compressor:** with a pressure drop per cooler the last stage now ends at the
  specified outlet pressure (it ended lower by the drops times the stage ratio).
- **Carboxylic acid dimers:** the activity-coefficient methods describe the vapor-phase
  dimerization of acetic, formic and propionic acid (chemical theory, published
  dimerization constants). The acids' heats of vaporization now equal the measured ones
  (they were 40 to 65 % too high), vapor densities include the dimers, and liquid-phase
  heats of reaction follow the formation enthalpies: the esterification of acetic acid with
  ethanol comes out at -5.2 instead of -13.8 kJ/mol (formation data: -4.5). Acetic
  acid/water no longer shows a minimum-boiling azeotrope near pure water, and the formic
  acid/water azeotrope lies at 382.3 K instead of 377.5 K (measured 380.5 K). The binary
  regression uses the same vapor, so fitted acid parameters stay consistent. The
  esterification examples change accordingly: their reaction is now almost thermoneutral,
  so the batch example cools to 300 K and needs 81 instead of 19 minutes to 60 %
  conversion. Other carboxylic acids, and acids with an equation of state, are flagged with a
  warning.
- An isentropic flash of a pure liquid with an activity-coefficient method no longer fails.
- An isentropic flash of a mixture into the two-phase region with an activity-coefficient
  method is solved by temperature instead of failing.
- Dynamic runs build their own finite-difference Jacobian with one step for the scaled
  states. The integrator's differences shrank to 1e-17 for states at zero (a controller's
  integral), where rounding noise became spurious couplings: a run could stall at an exact
  steady state or fail at an event. Vessel and reactor temperatures no longer depend on
  the calls before.
- A vented batch that sits on its bubble point without heat no longer switches between
  boiling and subcooled.

## [0.2.0] - 2026-10-30

First packaged release: `pip install "fluxsheet[server]"` installs the simulator, the web
editor and the examples.

### Added

- **Columns:** shortcut design (Fenske, Underwood, Gilliland/Molokanov, Kirkbride); rigorous
  equilibrium-stage distillation solved by Newton's method (Naphtali-Sandholm) with total,
  partial or no condenser, kettle reboiler, stage pressure drop and Murphree efficiency;
  absorber/stripper on the same stage model; stage profiles in the editor.
- **Reactions and reactors:** flowsheet-level reactions with element-balance checks, power-law
  kinetics and equilibrium constants; stoichiometric, yield, equilibrium, Gibbs, CSTR and plug
  flow reactors; energy balances on formation-enthalpy basis.
- **Thermodynamics:** NRTL, UNIQUAC and Wilson with ChemSep parameters; vapor-liquid-liquid
  equilibrium and a decanter; T-xy, P-xy and phase-envelope diagrams; property methods per
  unit, with streams re-flashed where methods meet.
- **Parameters fitted to NIST ThermoML data** for eight binary pairs (NRTL and UNIQUAC), with
  fit statistics and source DOIs; `fluxsheet fit-binary` fits further pairs. Results show where
  every binary parameter set comes from.
- **Heat exchanger** (two-sided, LMTD/UA, minimum approach), **design specs** that adjust a
  parameter to reach a target, and **sensitivity studies** with chart and CSV export
  (`fluxsheet sensitivity`).
- **Editor:** undo/redo, copy/paste/duplicate, auto-layout, keyboard shortcut help, unsaved
  changes warning, generic result fields and profile charts.
- Examples: methanol synthesis loop, steam reforming, ethyl acetate, benzene/toluene column,
  propane absorber, butanol decanter, feed/effluent exchanger and more.
- `FLUXSHEET_DEMO=1` limits for public servers (time budget per calculation, one calculation
  at a time, request and flowsheet size), and a Dockerfile for the online demo.

### Fixed

- Units downstream of a recycle loop are calculated once after the loop converges instead of
  in every iteration (a column behind a reactor loop failed on unconverged feeds).

## [0.1.0] - 2026-10-29

Initial prototype: compound search over ~76,000 compounds, Peng-Robinson, SRK, UNIFAC and
ideal property methods with an automatic suggestion, feeds, mixers, splitters, heaters, flash
drums, valves, pumps and compressors, recycle tearing with Wegstein acceleration, readable JSON
flowsheets with units, CLI, HTTP API and a React Flow editor.

[Unreleased]: https://github.com/almi2000/fluxsheet/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/almi2000/fluxsheet/releases/tag/v0.3.0
[0.2.0]: https://github.com/almi2000/fluxsheet/blob/v0.3.0/CHANGELOG.md#020---2026-10-30
