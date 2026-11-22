# Example Fluxsheet plug-in

A small package that adds three things to Fluxsheet through entry points:

| What | Entry point | Where |
| --- | --- | --- |
| `DewPointCooler` unit: cools a vapor to its dew point plus a margin | `fluxsheet.units` | `src/fluxsheet_example_plugin/units.py` |
| `regular-solution` property method (Hildebrand) | `fluxsheet.property_methods` | `src/fluxsheet_example_plugin/thermo.py` |
| `light-naphtha-cut` pseudo-component (Joback estimates) | `fluxsheet.components` | `src/fluxsheet_example_plugin/compounds.py` |

Install it next to Fluxsheet and check that it loaded:

```
pip install -e examples/plugin
fluxsheet plugins
```

Then run the two flowsheets in this folder, or open them in the editor (`fluxsheet serve`):

```
fluxsheet run examples/plugin/naphtha_cooler.json
fluxsheet run examples/plugin/regular_solution_flash.json
```

Both need the plug-in installed. Without it, Fluxsheet reports the unknown unit type, the
unknown property method or the unknown compound. See "Plug-ins" in
[docs/architecture.md](../../docs/architecture.md) for the rules.

The regular-solution model uses the solubility parameters and liquid volumes of the compound
database. It therefore rejects `light-naphtha-cut`, which has neither, with an error from the
plug-in.
