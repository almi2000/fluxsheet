# Fluxsheet

**An open, friendly chemical process simulator.** Draw a flowsheet in the browser, pick your
compounds from a database of ~76,000, and get rigorous phase-equilibrium results, with
error messages that tell you what to do next.

![Fluxsheet editor: an NGL recovery process with a gas recycle](docs/screenshot.jpg)

> **Status: early prototype (v0.1).** Steady-state, sequential-modular, vapor–liquid. Good for
> learning and conceptual design; not yet a replacement for commercial simulators.

## Why another simulator?

Commercial process simulators are powerful but hard to learn, expensive, and closed. Open
alternatives exist, but most are either libraries for programmers or desktop apps with dated
interfaces. Fluxsheet aims for both: a modern editor *and* a clean Python core.

| | |
|---|---|
| **Guides you through property methods** | Picks a thermodynamic model for your components and explains why (e.g. *"water and ethanol are polar, so the liquid is non-ideal, use UNIFAC"*), and warns when your manual choice looks wrong. |
| **Explains every error** | Each error names the unit or stream, says what went wrong, and suggests a fix: *"Valve V1: outlet pressure is above inlet pressure. Hint: a valve can only reduce pressure; use a pump or compressor."* |
| **Values with units** | Type `80 degC`, `1.5 bar`, `10 kmol/h` anywhere, in the editor or in files. |
| **Files you can diff** | Flowsheets are plain, readable JSON. Put them in Git, review changes in pull requests. |
| **Recycles just work** | Tear streams are found automatically and converged with Wegstein acceleration. |
| **Scriptable** | The editor is a thin client over a Python package you can use directly. |

## Quick start

Requires Python ≥ 3.10 and, for the editor, Node.js ≥ 18.

```bash
git clone https://github.com/almi2000/fluxsheet.git
cd fluxsheet
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e ".[server,dev]"
```

**Run a flowsheet from the command line**

```bash
fluxsheet run examples/ethanol_flash.json
```

```
Property method: unifac
  (water, ethanol are polar, so the liquid is likely non-ideal. ...)

Stream  T [°C]  P [bar]  Vap. frac.  Flow [kmol/h]  Flow [kg/h]
---------------------------------------------------------------
S1       25.00      1.2      0.0000            100         2363
S2       92.00      1.1      0.5052            100         2363
S3       92.00      1.1      1.0000          50.52         1396
S4       92.00      1.1      0.0000          49.48        966.4

HEATER: duty 733.5 kW
```

**Start the web editor**

```bash
cd web && npm install && npm run build && cd ..
fluxsheet serve                # open http://127.0.0.1:8000
```

For editor development, run `fluxsheet serve` and `npm run dev` (in `web/`) side by side,
then open http://localhost:5173.

**Use it from Python**

```python
from fluxsheet import FlowsheetModel, solve

result = solve(FlowsheetModel.load("examples/gas_recycle.json"))
print(result.streams["NGL-OUT"]["mass_flow"] * 3600, "kg/h")
```

## A flowsheet file

```json
{
  "format": "fluxsheet/1",
  "name": "Ethanol-water flash",
  "components": ["water", "ethanol"],
  "property_method": "auto",
  "units": [
    {"id": "FEED", "type": "Feed",
     "params": {"T": "25 degC", "P": "1.2 bar",
                "flows": {"water": "80 kmol/h", "ethanol": "20 kmol/h"}}},
    {"id": "HEATER", "type": "Heater", "params": {"T": "92 degC", "dP": "0.1 bar"}},
    {"id": "DRUM", "type": "Flash"}
  ],
  "streams": [
    {"id": "S1", "from": "FEED.out", "to": "HEATER.in"},
    {"id": "S2", "from": "HEATER.out", "to": "DRUM.in"},
    {"id": "S3", "from": "DRUM.vapor"},
    {"id": "S4", "from": "DRUM.liquid"}
  ]
}
```

Plain numbers are SI (K, Pa, mol/s, W). A stream without `to` leaves the process.

## What's included

**Unit operations:** Feed, Product, Mixer, Splitter, Heater/Cooler, Flash drum, Component
splitter, Valve, Pump, Compressor.

**Property methods:**

| Method | Use for |
|---|---|
| Peng-Robinson | Hydrocarbons, light gases, high pressure (ChemSep binary parameters) |
| Soave-Redlich-Kwong | Same territory, common in gas processing |
| UNIFAC (Dortmund) | Polar, non-ideal liquids at low/moderate pressure (water, alcohols, ketones…) |
| Ideal (Raoult) | Similar molecules at low pressure; single compounds |

Thermodynamics and compound data come from the excellent MIT-licensed
[`thermo`](https://github.com/CalebBell/thermo) and
[`chemicals`](https://github.com/CalebBell/chemicals) libraries by Caleb Bell. See
[docs/data-sources.md](docs/data-sources.md) for where the data comes from and how we plan to
extend it.

## Roadmap

- [ ] Distillation: shortcut column (Fenske–Underwood–Gilliland), then rigorous MESH column
- [ ] Reactors: stoichiometric, equilibrium (Gibbs), CSTR/PFR with kinetics
- [ ] NRTL / UNIQUAC with regressed binary parameters; parameter fitting from NIST ThermoML data
- [ ] Vapor–liquid–liquid equilibrium (decanters, three-phase flash)
- [ ] Heat exchanger (two-sided), design specs, sensitivity studies
- [ ] Phase envelopes, T-xy / P-xy diagrams in the editor
- [ ] Undo/redo, copy/paste, auto-layout in the editor
- [ ] CAPE-OPEN interoperability

## Project layout

```
src/fluxsheet/   Python core: compounds, thermo, units, solver, CLI, HTTP API
web/                 React + React Flow editor
examples/            Example flowsheets (also loadable in the editor)
tests/               pytest suite
docs/                Architecture and data-source notes
```

More in [docs/architecture.md](docs/architecture.md). Contributions are welcome, see
[CONTRIBUTING.md](CONTRIBUTING.md).

## License

[MIT](LICENSE). Fluxsheet is an independent project, not affiliated with or derived from
any commercial simulator.
