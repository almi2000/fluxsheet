# Example Fluxsheet plug-in

A small package that adds four things to Fluxsheet through entry points:

| What | Entry point | Where |
| --- | --- | --- |
| `DewPointCooler` unit: cools a vapor to its dew point plus a margin, and prices its exchanger (`estimate_cost`) | `fluxsheet.units` | `src/fluxsheet_example_plugin/units.py` |
| `regular-solution` property method (Hildebrand), with light gases as Henry components (`henry=True`) | `fluxsheet.property_methods` | `src/fluxsheet_example_plugin/thermo.py` |
| `light-naphtha-cut` pseudo-component (Joback estimates) | `fluxsheet.components` | `src/fluxsheet_example_plugin/compounds.py` |
| `CoolerLag`, the dew-point cooler's dynamic model: a first-order lag of its outlet | `fluxsheet.dynamic_models` | `src/fluxsheet_example_plugin/dynamics.py` |

Install it next to Fluxsheet and check that it loaded:

```
pip install -e examples/plugin
fluxsheet plugins
```

Then run the flowsheets in this folder, or open them in the editor (`fluxsheet serve`):

```
fluxsheet run examples/plugin/naphtha_cooler.json
fluxsheet dynamic examples/plugin/naphtha_cooler.json
fluxsheet run examples/plugin/regular_solution_flash.json
fluxsheet run examples/plugin/hydrogen_separator.json
```

The naphtha cooler has a costing section, which the cooler's `estimate_cost` answers with
its exchanger and cooling water, and a dynamic run in which the superheat steps from 5 to
15 K and the cooler's outlet follows with a 40 s lag. In the hydrogen separator, hydrogen
and methane dissolve in the hydrocarbon liquid by Henry's law, so the regular-solution model
needs no solubility parameters for them.

All need the plug-in installed. Without it, Fluxsheet reports the unknown unit type, the
unknown property method or the unknown compound. See "Plug-ins" in
[docs/architecture.md](../../docs/architecture.md) for the rules.

The regular-solution model uses the solubility parameters and liquid volumes of the compound
database. It therefore rejects `light-naphtha-cut`, which has neither, with an error from the
plug-in.
