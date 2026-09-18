#!/usr/bin/env python3
"""Build the auditable ATAC production manifest from the complete catalogue.

The output contains every ATAC FASTQ in the full manifest.  Eligibility is
recomputed with the shared selection gate; validation subset files are never
consulted.  Structural uncertainty is represented as a closed decision.
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from selection_gate import catalogue, decision, fingerprint
from file_hash import sha256_file

EXPECTED_ROLES = {"SINGLE": {"SINGLE"}, "PAIRED": {"R1", "R2"}}
REQUIRED = {
    "run_accession", "study_accession", "sample_accession", "omics",
    "library_layout", "processing_layout", "fastq_structure", "fastq_role",
    "filename", "fastq_ftp", "fastq_md5", "fastq_bytes",
}


def build_rows(manifest_path, repo):
    with Path(manifest_path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED - set(reader.fieldnames or ())
        if missing:
            raise ValueError("Full FASTQ manifest lacks columns: " + ",".join(sorted(missing)))
        rows = [dict(row) for row in reader if row.get("omics") == "ATAC-seq"]
    if not rows:
        raise ValueError("Full FASTQ manifest contains no ATAC-seq records")

    groups = defaultdict(list)
    for row in rows:
        groups[row["run_accession"]].append(row)

    output = []
    policy_records = catalogue(repo)
    for run, group in sorted(groups.items()):
        ok, selection_reason = decision(run, repo=repo, omics="ATAC-seq")
        structures = {r["fastq_structure"] for r in group if r["fastq_structure"]}
        layouts = {r["processing_layout"] for r in group if r["processing_layout"]}
        library_layouts = {r["library_layout"] for r in group if r["library_layout"]}
        roles = {r["fastq_role"] for r in group if r["fastq_role"]}
        identities = {(r["study_accession"], r["sample_accession"]) for r in group}
        problems = []
        if len(structures) != 1 or next(iter(structures), "") not in EXPECTED_ROLES:
            problems.append("unsupported_or_ambiguous_fastq_structure")
        structure = next(iter(structures), "UNRESOLVED")
        if roles != EXPECTED_ROLES.get(structure, set()):
            problems.append("inconsistent_fastq_roles")
        if len(layouts) != 1 or next(iter(layouts), "") != structure:
            problems.append("inconsistent_processing_layout")
        if len(library_layouts) != 1 or len(identities) != 1:
            problems.append("ambiguous_run_metadata")
        if any(not r[k].strip() for r in group for k in REQUIRED):
            problems.append("missing_metadata")
        if len({r["filename"] for r in group}) != len(group):
            problems.append("duplicate_fastq_identity")

        eligible = ok and not problems
        reason = "eligible" if eligible else ";".join(
            ([selection_reason] if not ok else []) + problems
        )
        for row in group:
            row.update(
                selection_eligible="true" if eligible else "false",
                selection_decision="eligible" if eligible else "excluded_or_review",
                selection_reason=reason,
                selection_policy=policy_records.get(run, {}).get("selection_policy", "unknown"),
                selection_policy_source="workflow/scripts/selection_policy.py",
            )
            output.append(row)
    return output


def write_manifest(rows, output):
    fields = list(rows[0])
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--provenance", required=True)
    parser.add_argument("--repo", default=Path(__file__).resolve().parents[2], type=Path)
    args = parser.parse_args()
    rows = build_rows(args.manifest, args.repo)
    write_manifest(rows, args.output)
    runs = sorted({r["run_accession"] for r in rows})
    eligible = sorted({r["run_accession"] for r in rows if r["selection_eligible"] == "true"})
    provenance = {
        "schema_version": 1,
        "source_manifest": str(Path(args.manifest)),
        "source_manifest_sha256": sha256_file(args.manifest),
        "selection_inputs": fingerprint(args.repo),
        "atac_runs": runs,
        "eligible_runs": eligible,
        "fastq_records": len(rows),
    }
    out = Path(args.provenance)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
