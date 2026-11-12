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
        if record.get("admission") not in {"incomplete", "incompatible"}:
            raise ValueError(f"{name}: no qualified simulator adapter is registered yet")
        for key in ("url", "evidence", "reason"):
            if not isinstance(record.get(key), str) or not record[key].strip():
                raise ValueError(f"{name}: missing {key}")
        missing = record.get("missing_inputs")
        if not isinstance(missing, list) or not missing or any(
                not isinstance(s, str) or not s.strip() for s in missing):
            raise ValueError(f"{name}: missing input requirements")
        sha = record.get("sha256")
        if sha is not None and (not isinstance(sha, str) or not re.fullmatch("[0-9a-f]{64}", sha)):
            raise ValueError(f"{name}: invalid SHA-256")
        results.append({**record, "screening": "not_checkable", "passed": None,
                        "experimental_validation": "open"})
    return results


def markdown(records):
    rows = ["## Equipment source qualification", ""]
    for record in records:
        rows.extend([f"### {record['id']} ({record['family']})", "",
                     f"Source: {record['url']}", "",
                     (f"Admission: {record['admission']}; screening: not_checkable; "
                      f"independence: {record['independence']}; experimental validation: open."),
                     "", record["reason"], "", f"Evidence: {record['evidence']}", "",
                     f"SHA-256: {record.get('sha256') or 'not established'}", "",
                     "Missing inputs:", ""])
        rows.extend(f"- [ ] {item}" for item in record["missing_inputs"])
        rows.append("")
    return rows
