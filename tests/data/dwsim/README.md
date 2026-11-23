# DWSIM sample files

Small simulations in DWSIM's `.dwxml` layout, for `tests/test_import_dwsim.py`. They were
written for Fluxsheet (not exported from DWSIM, and no DWSIM file was copied), following
the element names DWSIM's source code writes (`DWSIM.FlowsheetBase` SaveToXML, the
unit operations' and streams' SaveData). Only the parts the importer reads are present;
real DWSIM files also carry every compound's property data, settings and results.

The stored stream states (temperature, pressure, flows, composition), which DWSIM writes
after a run, were computed with Fluxsheet, not with DWSIM.

- `heater_flash.dwxml`: water/ethanol, NRTL. Feed, heater (outlet temperature) with an
  energy stream, flash drum, a text annotation.
- `gas_recycle.dwxml`: light hydrocarbons, Peng-Robinson. Mixer, cooler, flash, splitter,
  compressor and a Recycle block closing the loop, plus an Adjust block.
- `shortcut_column.dwxml`: benzene/toluene, Raoult's law. Pump, shortcut column, a pipe
  (no Fluxsheet counterpart) and a valve.
- `hydrogenation.dwxml`: benzene hydrogenation, Peng-Robinson. Feed/effluent heat exchanger
  (hot side on connector 1), heater, conversion reactor with its reaction set, cooler,
  flash and a compound separator. The tests also zip it into a `.dwxmz`.

## Real DWSIM files

`real-samples.json` identifies three of the sample simulations that ship with DWSIM (its
repository at a fixed commit, by SHA-256). They are GPL-licensed and are not copied into
this repository. `python scripts/fetch_dwsim_samples.py --cache-dir PATH` downloads and
verifies them; `FLUXSHEET_DWSIM_SAMPLES=PATH pytest tests/test_import_dwsim.py` then imports
each one, solves it and compares every stream with the state DWSIM stored. Without the
variable those tests are skipped.
