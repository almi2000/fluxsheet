"""Validate source qualification records without implying experimental success."""

import re

FAMILIES = {"dynamics", "kinetics", "column"}
INDEPENDENCE = {"separate", "shared_fit_data", "unconfirmed"}


def qualify(records, occupied_ids=()):
    if not isinstance(records, list):
        raise TypeError("equipment_sources must be a list")
    seen = set(occupied_ids)
    results = []
    for record in records:
        if not isinstance(record, dict):
            raise TypeError("equipment source must be an object")
        name = record.get("id")
        if not isinstance(name, str) or not name.strip() or name in seen:
            raise ValueError("missing or duplicate equipment source id")
        seen.add(name)
        if record.get("family") not in FAMILIES:
            raise ValueError(f"{name}: unsupported equipment family")
        if record.get("independence") not in INDEPENDENCE:
            raise ValueError(f"{name}: missing independence classification")
        admission = record.get("admission")
        for key in ("url", "evidence"):
            if not isinstance(record.get(key), str) or not record[key].strip():
                raise ValueError(f"{name}: missing {key}")
        sha = record.get("sha256")
        if sha is not None and (not isinstance(sha, str) or not re.fullmatch("[0-9a-f]{64}", sha)):
            raise ValueError(f"{name}: invalid SHA-256")
        if admission == "qualified" and record["family"] == "dynamics":
            _check_runs(record)
            results.append({**record, "screening": "pending", "passed": None,
                            "experimental_validation": "open"})
        else:
            if admission not in {"incomplete", "incompatible"}:
                raise ValueError(f"{name}: unsupported admission or no simulator adapter")
            if not isinstance(record.get("reason"), str) or not record["reason"].strip():
                raise ValueError(f"{name}: missing reason")
            missing = record.get("missing_inputs")
            if not isinstance(missing, list) or not missing or any(
                    not isinstance(s, str) or not s.strip() for s in missing):
                raise ValueError(f"{name}: missing input requirements")
            results.append({**record, "screening": "not_checkable", "passed": None,
                            "experimental_validation": "open"})
    return results


def _check_runs(record):
    name = record["id"]
    kind = record.get("data_kind")
    if kind not in {"measured", "synthetic"}:
        raise ValueError(f"{name}: data_kind must be measured or synthetic")
    if not isinstance(record.get("parameter_provenance"), str) or not record["parameter_provenance"].strip():
        raise ValueError(f"{name}: missing parameter provenance")
    if kind == "measured" and not record.get("sha256"):
        raise ValueError(f"{name}: missing original source SHA-256")
    if record.get("missing_inputs"):
        raise ValueError(f"{name}: qualified source has unresolved inputs")
    runs = record.get("runs")
    minimum = 2 if kind == "measured" else 1
    if not isinstance(runs, list) or len(runs) < minimum:
        raise ValueError(f"{name}: requires at least {minimum} registered runs")
    ids = set()
    measurement_hashes = set()
    for run in runs:
        run_id = run.get("id") if isinstance(run, dict) else None
        if not isinstance(run_id, str) or not run_id.strip() or run_id.strip() in ids:
            raise ValueError(f"{name}: missing or duplicate run id")
        ids.add(run_id.strip())
        reference = run.get("measurements")
        sha = reference.get("sha256") if isinstance(reference, dict) else None
        if kind == "measured" and isinstance(sha, str):
            if sha in measurement_hashes:
                raise ValueError(f"{name}: measured runs must use distinct measurement files")
            measurement_hashes.add(sha)


def markdown(records):
    rows = ["## Equipment source qualification", ""]
    for record in records:
        rows.extend([f"### {record['id']} ({record['family']})", "",
                     f"Source: {record['url']}", "",
                     (f"Admission: {record['admission']}; screening: {record['screening']}; "
                      f"independence: {record['independence']}; experimental validation: "
                      f"{record['experimental_validation']}."), "",
                     f"Evidence: {record['evidence']}", "",
                     f"SHA-256: {record.get('sha256') or 'not established'}", ""])
        if record["admission"] == "qualified":
            rows.extend([f"Data: {record['data_kind']}", "",
                         f"Parameter provenance: {record['parameter_provenance']}", ""])
        else:
            rows.extend([record["reason"], "", "Missing inputs:", ""])
            rows.extend(f"- [ ] {item}" for item in record["missing_inputs"])
        rows.append("")
    return rows
