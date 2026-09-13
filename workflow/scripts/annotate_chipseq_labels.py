#!/usr/bin/env python3
"""Extract declared histone targets and control roles from verified ENA XML."""
import argparse
import csv
import io
import json
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from fetch_chipseq_ena_xml import digest, validate_xml, write_atomic


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def extract(texts, rules):
    hits = []
    for source in rules["source_fields"]:
        text = texts[source]
        for match in re.finditer(rules["histone_pattern"], text, re.I):
            target = (
                match["histone"].upper() + match["residue"].upper()
                + match["position"] + match["modification"].lower()
            )
            hits.append(
                (source, "histone_label", "target", target, match.group())
            )
        for rule in rules["control_patterns"]:
            for match in re.finditer(rule["pattern"], text, re.I):
                hits.append((
                    source, rule["rule_id"], "control",
                    rule["role"], match.group()
                ))

    targets = sorted({h[3] for h in hits if h[2] == "target"})
    controls = sorted({h[3] for h in hits if h[2] == "control"})

    target, role, status = "", "unresolved", "unresolved"
    if len(targets) > 1 or len(controls) > 1 or (targets and controls):
        status = "conflict"
    elif controls:
        role, status = controls[0], "resolved_from_labels"
    elif targets:
        target, role, status = targets[0], "ip", "resolved_from_labels"

    return target, role, status, targets, controls, hits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("runs", "inventory", "snapshot", "rules", "outdir"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()

    rules = json.loads(args.rules.read_text())
    if rules["schema_version"] != 1:
        raise ValueError("Unsupported rule schema")

    runs = read_tsv(args.runs)
    inventory_rows = read_tsv(args.inventory)
    inventory = {r["accession"]: r for r in inventory_rows}

    if len(inventory) != len(inventory_rows):
        raise ValueError("Duplicate inventory accessions")
    if not runs or len({r["run_accession"] for r in runs}) != len(runs):
        raise ValueError("Empty run table or duplicate runs")

    cache = {}

    def record(accession, kind):
        if accession not in cache:
            item = inventory[accession]
            if item["record_type"] != kind or item["status"] != "verified":
                raise ValueError(f"Invalid inventory record: {accession}")
            data = (
                args.snapshot / "records" / accession / "record.xml"
            ).read_bytes()
            if digest(data) != item["sha256"]:
                raise ValueError(f"XML hash mismatch: {accession}")

            validate_xml(data, kind, accession)
            root = ET.fromstring(data)
            cache[accession] = (
                root if root.tag == kind.upper() else root.find(kind.upper())
            )
        return cache[accession]

    annotations, evidence = [], []

    for row in sorted(
        runs, key=lambda r: (r["study_accession"], r["run_accession"])
    ):
        sample = record(row["sample_accession"], "sample")
        experiment = record(row["experiment_accession"], "experiment")

        texts = {
            "sample_alias": sample.get("alias", ""),
            "sample_title": sample.findtext("TITLE") or "",
            "sample_description": sample.findtext("DESCRIPTION") or "",
            "experiment_title": experiment.findtext("TITLE") or "",
            "library_name": experiment.findtext(
                "./DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME"
            ) or "",
        }

        target, role, status, targets, controls, hits = extract(texts, rules)

        annotations.append({
            "run_accession": row["run_accession"],
            "study_accession": row["study_accession"],
            "experiment_accession": row["experiment_accession"],
            "sample_accession": row["sample_accession"],
            "declared_target": target,
            "library_role": role,
            "label_status": status,
            "target_candidates": json.dumps(targets),
            "control_role_candidates": json.dumps(controls),
        })

        for source, rule_id, kind, value, matched in hits:
            accession = (
                row["sample_accession"] if source.startswith("sample_")
                else row["experiment_accession"]
            )
            item = inventory[accession]
            evidence.append({
                "run_accession": row["run_accession"],
                "source_field": source,
                "source_accession": accession,
                "source_url": item["source_url"],
                "source_sha256": item["sha256"],
                "rule_id": rule_id,
                "extracted_property": kind,
                "extracted_value": value,
                "matched_text": matched,
                "source_text": texts[source],
            })

    args.outdir.mkdir(parents=True, exist_ok=True)
    evidence_fields = [
        "run_accession", "source_field", "source_accession", "source_url",
        "source_sha256", "rule_id", "extracted_property", "extracted_value",
        "matched_text", "source_text"
    ]

    for name, rows, fields in (
        ("chipseq_target_annotations.tsv", annotations, list(annotations[0])),
        ("chipseq_label_evidence.tsv", evidence, evidence_fields),
    ):
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(
            buffer, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
        write_atomic(args.outdir / name, buffer.getvalue().encode())

    summary = {
        "run_count": len(annotations),
        "role_counts": dict(Counter(
            r["library_role"] for r in annotations
        )),
        "label_status_counts": dict(Counter(
            r["label_status"] for r in annotations
        )),
        "target_counts": dict(Counter(
            r["declared_target"] for r in annotations if r["declared_target"]
        )),
        "rules_sha256": digest(args.rules.read_bytes()),
        "runs_sha256": digest(args.runs.read_bytes()),
        "inventory_sha256": digest(args.inventory.read_bytes()),
        "script_sha256": digest(Path(__file__).read_bytes()),
        "control_pairing": "not_performed",
    }
    serialized = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    write_atomic(
        args.outdir / "chipseq_label_summary.json", serialized.encode()
    )
    print(serialized)


if __name__ == "__main__":
    main()
