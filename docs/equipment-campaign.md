# Equipment campaign, manifest version 2

Run the current qualification register without measurements or network access:

```powershell
.\.venv\Scripts\python.exe scripts/validate_models.py --manifest docs/equipment-campaign.json --offline --cache-dir "$env:TEMP/fluxsheet-validation-cache" --output-dir "$env:TEMP/fluxsheet-equipment-output"
```

The expected exit code is **2**. Seven source assessments remain incomplete.
`validation.json` includes `equipment_sources` with admission, independence,
evidence, source hashes where available, and missing inputs. `validation.md`
renders the same records and data-request checklists. These are qualification
results, not numerical comparisons: `passed` is null and screening is
`not_checkable`. A missing hash remains null, never a fabricated value.

Version 2 retains the version-1 `sources` array for existing thermodynamic
comparisons and adds `equipment_sources`. Both may coexist; incomplete equipment
records force exit 2 even when other comparisons pass or fail. Existing comparisons
and errors remain in the report. Source identifiers must be unique across both
arrays. Version 1 rejects equipment fields rather than silently ignoring them.

Each equipment record requires `id`, `family` (dynamics/kinetics/column), `url`,
`admission`, `evidence`, and `independence` (separate/shared_fit_data/unconfirmed).
For `incomplete` or `incompatible`, supply `reason` and `missing_inputs`.
`sha256` is optional or null for open sources. Recorded hashes identify the original
sources. Download mode applies to the legacy NIST sources. Equipment files are supplied
through the external cache.

## Running a stirred vessel comparison

Use the command above with your manifest. For an executable source set `family`
to `dynamics` and `admission` to `qualified`. Set `data_kind` to `measured` or
`synthetic`, and supply `parameter_provenance` explaining the origin of the model
parameters and whether they were derived from the comparison runs. Document the
source of the initial conditions, cylindrical geometry, constant pressures and
net heat delivered to the liquid. Net heat may be based on independently established
losses if the source measures gross power. The model assumes negligible wall storage.
For measured sources supply the original source SHA-256 and at least two distinct
runs. Executable kinetics and column cases remain pending.

Each entry in the source's `runs` list has these fields:

| Field | Format |
| --- | --- |
| `id` | Unique experiment identifier within the source |
| `flowsheet`, `measurements` | Objects with cache-relative `path` and lowercase `sha256`. The two files are a Fluxsheet JSON flowsheet and a UTF-8, comma-delimited CSV. |
| `time` | `{"column": "time", "unit": "s"}`. The common time axis starts at zero and strictly increases. |
| `feeds` | Object keyed by each Feed ID. Each entry has a `temperature` and `flow` mapping. |
| `duty`, `liquid_flow` | Input mappings for net liquid heat duty and imposed outlet flow |
| Input mapping | For example `{"column": "inlet_flow", "unit": "kmol/h", "interpolation": "step"}`. Use `step` for values held from their timestamp and `linear` for ramps between samples. Flows are molar. |
| `outputs` | One or both named mappings `temperature` and `level`, each with `column`, `unit` and `uncertainty`. Use a nonnegative absolute uncertainty in the channel unit or `null` when unknown. |
| `initial` | `temperature_K`, `temperature_tolerance_K`, `inventory_quantity` (`level`, `holdup`, or `liquid_volume`), corresponding `inventory_si` and `inventory_tolerance_si`, and source `evidence` for the values and tolerances |
| `baseline_s`, `comparison_s` | Inclusive `[start, end]` intervals in seconds. Include at least three baseline samples before the comparison. |
| `end_start_s` | Start of the registered settled interval within the comparison window. Include at least two samples from here through the comparison end. |

The flowsheet contains built-in `water`, direct Feed connections to one
`MixingTank`, and a liquid outlet. Set the feed and vessel pressures explicitly.
Register the tank under `dynamics.units` with `diameter`, `height` and
`initial_level` (fraction of height). The adapter uses the CSV to set feed
temperatures and flows, fixed liquid outflow and tank duty. Use an empty event
and controller list. The tank must stay liquid, within its height, and at or
above 1 cm, where the simulator can deliver the specified outflow.

The first CSV row supplies the initial boundary conditions for the stationary
solve. The adapter checks the resulting initial temperature and inventory
against the registered tolerances. It checks the signal against five times the
detrended baseline noise and compares results at the original CSV timestamps.
The maximum prediction difference between 501 and 1001 output points with
relative tolerances of `1e-6` and `1e-7` must stay within 1% of the measured step
amplitude. The fine result must meet both 10% limits for RMSE and end-window
error. Output errors and uncertainties are reported in SI.

Numerical results appear under `equipment_sources[].run_results` in JSON, with
the run definition, file hashes, original comparison times, predictions,
residuals and uncertainties. Markdown summarizes the limits and the largest
errors. A complete passing synthetic or reproduction campaign can return code
`0` while experimental validation stays open. The independent status applies
only to passing measured runs with documented separate provenance and coverage
of temperature and level. A failed threshold returns `1`. Unavailable or
non-checkable cases return `2` and retain completed comparisons.

The stirred-vessel adapter has passed the synthetic CLI checks for a heat input
step and an inlet temperature step. The independent experimental comparison of
temperature and level remains open pending suitable measured runs. The Dynamics
tab shows this status when a stirred tank is configured as a dynamic unit.

See [source details](equipment-reference-qualification.md) for the current evidence
and [model validation](model-validation.md) for fixed engineering limits.
