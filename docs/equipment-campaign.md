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
`admission` (incomplete/incompatible), `evidence`, `reason`, `missing_inputs`, and
`independence` (separate/shared_fit_data/unconfirmed). `sha256` is optional or null;
a supplied value must be a lowercase SHA-256 digest. Recorded hashes identify
previously inspected archives; this qualification-only command does not fetch or
reverify those archives. Download mode still applies only to legacy NIST sources.

**Current implementation boundary:** qualified apparatus simulation adapters and
measurement-file schemas are not yet registered. An admission of `qualified`
therefore fails explicitly, rather than returning a false pass. The numerical
helpers in `scripts/equipment_metrics.py` remain separately tested building blocks.
No parameter fitting or experimental model calculation is performed by this register.

See [source details](equipment-reference-qualification.md) for the current evidence
and [model validation](model-validation.md) for fixed engineering limits.
