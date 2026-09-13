#!/usr/bin/env python3
"""Audit label-based control candidates; never authorize peak calling."""
import argparse
import csv
import io
import json
import re
import xml.etree.ElementTree as ET
from collections import Counter
from decimal import Decimal
from pathlib import Path
from fetch_chipseq_ena_xml import digest, validate_xml, write_atomic


def table(path, key):
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if any(not row[key] for row in rows):
        raise ValueError(f"Empty identifier: {path}")
    result = {row[key]: row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f"Duplicate identifiers: {path}")
    return result


def extract(texts, assay, patterns, study):
    hits = []
    for rule in patterns:
        pattern = rule["pattern"].replace("{assay}", re.escape(assay))
        match = re.fullmatch(pattern, texts[rule["field"]], flags=re.I)
        if not match:
            continue
        value, unit = match["value"], rule["unit"]
        if unit == "label":
            value = value.upper()
            mapping = study.get("interpretation", {}).get("time_label_to_days", {})
            if value in mapping:
                value, unit = str(mapping[value]), "days"
        if unit != "label":
            value = format(Decimal(value).normalize(), "f")
        hits.append({
            "rule_id": rule["id"], "source_field": rule["field"],
            "source_text": texts[rule["field"]],
            "context": " ".join(match["context"].casefold().split()),
            "time_value": value, "time_unit": unit,
        })
    return hits


def candidates(row, controls):
    fields = ("study_accession", "condition_rule", "label_context",
              "time_value", "time_unit", "library_layout", "instrument_platform")
    if row["condition_status"] != "resolved_from_labels":
        return [], "condition_unresolved"
    if not all(row[field] for field in fields):
        return [], "matching_fields_missing"
    matches = [c for c in controls
               if c["condition_status"] == "resolved_from_labels"
               and all(c[field] == row[field] for field in fields)]
    status = ("no_candidate" if not matches else
              "single_candidate_needs_review" if len(matches) == 1 else
              "multiple_candidates_needs_review")
    return matches, status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("runs", "snapshot", "rules", "study-evidence", "outdir"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    rules = json.loads(args.rules.read_text())
    studies = json.loads(args.study_evidence.read_text())
    if rules["schema_version"] != 1 or studies["schema_version"] != 1:
        raise ValueError("Unsupported schema")
    labeldir = args.snapshot / "annotations"
    paths = {
        "runs": args.runs, "rules": args.rules,
        "study_evidence": args.study_evidence,
        "inventory": args.snapshot / "xml_inventory.tsv",
        "annotations": labeldir / "chipseq_target_annotations.tsv",
        "label_summary": labeldir / "chipseq_label_summary.json",
        "script": Path(__file__),
        "validator": Path(__file__).with_name("fetch_chipseq_ena_xml.py"),
    }
    previous = json.loads(paths["label_summary"].read_text())
    for key in ("runs", "inventory"):
        if digest(paths[key].read_bytes()) != previous[key + "_sha256"]:
            raise ValueError(f"Label provenance mismatch: {key}")
    runs = table(args.runs, "run_accession")
    labels = table(paths["annotations"], "run_accession")
    inventory = table(paths["inventory"], "accession")
    if set(runs) != set(labels) or len(labels) != previous["run_count"]:
        raise ValueError("Run catalog mismatch")
    for column, key in (("library_role", "role_counts"),
                        ("declared_target", "target_counts"),
                        ("label_status", "label_status_counts")):
        counts = dict(Counter(r[column] for r in labels.values() if r[column]))
        if counts != previous[key]:
            raise ValueError(f"Label summary mismatch: {key}")

    def record(acc, kind):
        item = inventory[acc]
        data = (args.snapshot / "records" / acc / "record.xml").read_bytes()
        if (item["status"] != "verified" or item["record_type"] != kind
                or digest(data) != item["sha256"]):
            raise ValueError(f"XML inventory mismatch: {acc}")
        validate_xml(data, kind, acc)
        root = ET.fromstring(data)
        return root if root.tag == kind.upper() else root.find(kind.upper())

    conditions = []
    for acc, row in sorted(labels.items()):
        raw = runs[acc]
        for field in ("study_accession", "sample_accession", "experiment_accession"):
            if row[field] != raw[field]:
                raise ValueError(f"Run identity mismatch: {acc}/{field}")
        sample = record(row["sample_accession"], "sample")
        experiment = record(row["experiment_accession"], "experiment")
        study = studies["studies"].get(row["study_accession"], {})
        assay = "Input" if row["library_role"] == "input" else row["declared_target"]
        texts = {"sample_alias": sample.get("alias", ""),
                 "experiment_title": experiment.findtext("TITLE") or ""}
        hits = extract(texts, assay, rules["patterns"], study) if assay else []
        if row["label_status"] != "resolved_from_labels":
            hits = []
        for hit in hits:
            source = (row["sample_accession"] if hit["source_field"] == "sample_alias"
                      else row["experiment_accession"])
            hit.update(source_accession=source, source_url=inventory[source]["source_url"],
                       source_sha256=inventory[source]["sha256"])
        hit = hits[0] if len(hits) == 1 else {}
        flagged = any(acc in case.get("run_accessions", [])
                      for case in study.get("unresolved_cases", []))
        conditions.append({
            **{k: row[k] for k in ("run_accession", "study_accession",
               "sample_accession", "experiment_accession", "library_role", "declared_target")},
            "library_layout": raw["library_layout"],
            "instrument_platform": raw["instrument_platform"],
            "condition_status": ("resolved_from_labels" if len(hits) == 1 else
                                 "multiple_pattern_matches" if hits else "unresolved"),
            "condition_rule": hit.get("rule_id", ""),
            "label_context": hit.get("context", ""),
            "time_value": hit.get("time_value", ""),
            "time_unit": hit.get("time_unit", ""),
            "documented_issue": str(flagged).lower(),
            "evidence_json": json.dumps(hits, sort_keys=True),
        })

    controls = [r for r in conditions if r["library_role"] == "input"]
    audit = []
    for row in conditions:
        if row["library_role"] != "ip":
            continue
        matches, status = candidates(row, controls)
        audit.append({
            "ip_run_accession": row["run_accession"],
            "study_accession": row["study_accession"],
            "candidate_status": status,
            "candidate_count": len(matches),
            "candidate_input_runs": json.dumps([m["run_accession"] for m in matches]),
            "candidates_with_documented_issue": json.dumps([
                m["run_accession"] for m in matches if m["documented_issue"] == "true"]),
            "pairing_approved": "false",
        })
    if not conditions or not audit:
        raise ValueError("No conditions or IP records to audit")
    args.outdir.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for name, rows in (("chipseq_conditions.tsv", conditions),
                       ("chipseq_control_candidates.tsv", audit)):
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(buffer, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        data = buffer.getvalue().encode()
        write_atomic(args.outdir / name, data)
        outputs[name] = digest(data)
    summary = {
        "schema_version": 1, "run_count": len(conditions), "ip_count": len(audit),
        "condition_status_counts": dict(Counter(r["condition_status"] for r in conditions)),
        "candidate_status_counts": dict(Counter(r["candidate_status"] for r in audit)),
        "approved_pairings": 0,
        "input_sha256": {key: digest(path.read_bytes()) for key, path in paths.items()},
        "output_sha256": outputs,
    }
    text = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    write_atomic(args.outdir / "control_candidate_summary.json", text.encode())
    print(text)


if __name__ == "__main__":
    main()
