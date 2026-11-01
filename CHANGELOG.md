# Changelog

All notable changes to Fluxsheet. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/) (before 1.0, minor versions may change the file
format or API).

## [Unreleased]

### Added

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

### Fixed

- A pressure-enthalpy flash that thermo cannot converge (e.g. an almost pure component in a
  mixture model) is solved by temperature with T-P flashes instead of failing.
- A long flowsheet name no longer pushes the editor wider than the window.

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

[0.2.0]: https://github.com/almi2000/fluxsheet/releases/tag/v0.2.0
