#!/usr/bin/env python3

from pathlib import Path
import csv
import sys

ROOT = Path(__file__).resolve().parents[1]

STUDIES = ROOT / "sources" / "studies.tsv"
OBS = ROOT / "curated" / "benchmark_observations.tsv"

TARGET = "GCF_003668045.3"

ALLOWED_ROLES = {
    "positive",
    "negative",
    "support_only",
}

ALLOWED_CLASSES = {
    "targeted_validated",
    "screen_validated",
    "observed_stable",
    "observed_unstable",
    "predicted_candidate",
    "orthogonal_support",
}

ALLOWED_MAPPING = {
    "",
    "native_target_reference",
    "mapped_unique",
    "mapped_ambiguous",
    "unresolved",
}

def read_tsv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

errors = []

studies = read_tsv(STUDIES)
observations = read_tsv(OBS)

study_ids = [x["study_id"] for x in studies]

if len(study_ids) != len(set(study_ids)):
    errors.append("Duplicate study_id in studies.tsv")

known_studies = set(study_ids)

for i, row in enumerate(observations, start=2):
    oid = row["observation_id"] or f"line_{i}"

    if row["study_id"] not in known_studies:
        errors.append(f"{oid}: unknown study_id {row['study_id']}")

    if row["benchmark_role"] and row["benchmark_role"] not in ALLOWED_ROLES:
        errors.append(
            f"{oid}: invalid benchmark_role={row['benchmark_role']}"
        )

    if row["evidence_class"] and row["evidence_class"] not in ALLOWED_CLASSES:
        errors.append(
            f"{oid}: invalid evidence_class={row['evidence_class']}"
        )

    if row["mapping_status"] not in ALLOWED_MAPPING:
        errors.append(
            f"{oid}: invalid mapping_status={row['mapping_status']}"
        )

    if row["evidence_class"] in {
        "predicted_candidate",
        "orthogonal_support",
    }:
        if row["benchmark_role"] not in {"", "support_only"}:
            errors.append(
                f"{oid}: prediction-only evidence cannot be gold positive/negative"
            )

    target_fields = [
        row["target_seqname"],
        row["target_start"],
        row["target_end"],
    ]

    if any(target_fields):
        if not all(target_fields):
            errors.append(
                f"{oid}: incomplete target coordinates"
            )

        if row["target_assembly"] != TARGET:
            errors.append(
                f"{oid}: target assembly must be {TARGET}"
            )

        try:
            start = int(row["target_start"])
            end = int(row["target_end"])
            if start < 1 or end < start:
                errors.append(
                    f"{oid}: invalid 1-based closed coordinates"
                )
        except ValueError:
            errors.append(
                f"{oid}: target coordinates must be integers"
            )

if errors:
    print("BENCHMARK VALIDATION: FAIL")
    for error in errors:
        print(f" - {error}")
    sys.exit(1)

print("BENCHMARK VALIDATION: PASS")
print(f"Studies registered: {len(studies)}")
print(f"Experimental observations curated: {len(observations)}")
print(f"Canonical target reference: {TARGET}")
