"""Strict readers for the ATAC production manifest and replicate groups."""
import csv
from pathlib import Path

EXPECTED_ROLES = {"SINGLE": {"SINGLE"}, "PAIRED": {"R1", "R2"}}


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def eligible_run(manifest, run, require_single=False):
    rows = [r for r in read_tsv(manifest) if r.get("run_accession") == run]
    if not rows:
        raise ValueError(f"ATAC run {run} is absent from the production manifest")
    decisions = {(r.get("selection_eligible"), r.get("selection_decision"), r.get("selection_reason")) for r in rows}
    if len(decisions) != 1 or next(iter(decisions))[0] != "true":
        reason = ";".join(sorted({r.get("selection_reason", "missing_reason") for r in rows}))
        raise ValueError(f"ATAC run {run} is not production eligible: {reason}")
    structures = {r.get("fastq_structure") for r in rows}
    roles = {r.get("fastq_role") for r in rows}
    if len(structures) != 1:
        raise ValueError(f"Ambiguous ATAC structure for {run}")
    structure = next(iter(structures))
    if structure not in EXPECTED_ROLES or roles != EXPECTED_ROLES[structure]:
        raise ValueError(f"Invalid ATAC FASTQ roles for {run}: {structure}/{sorted(roles)}")
    if require_single and structure != "SINGLE":
        raise ValueError(f"ATAC downstream is validated only for SINGLE; {run} is {structure}")
    return rows


def eligible_runs(manifest, require_single=False):
    runs = sorted({r["run_accession"] for r in read_tsv(manifest) if r.get("selection_eligible") == "true"})
    for run in runs:
        eligible_run(manifest, run, require_single=require_single)
    if not runs:
        raise ValueError("No eligible ATAC production runs")
    return runs


def role_path(manifest, run, role, fastp_dir):
    rows = eligible_run(manifest, run)
    match = [r for r in rows if r["fastq_role"] == role]
    if len(match) != 1:
        raise ValueError(f"Expected one {role} FASTQ for {run}")
    structure = match[0]["fastq_structure"]
    base = match[0]["filename"].removesuffix(".fastq.gz")
    return f"{fastp_dir}/{run}/{structure}/{base}.fastq.gz"


def load_groups(path, manifest):
    rows = read_tsv(path)
    required = {"replicate_group", "run_accession", "compatibility", "replicate_type"}
    if not rows or required - set(rows[0]):
        raise ValueError("Invalid ATAC replicate-group table")
    groups = {}
    for row in rows:
        group = row["replicate_group"].strip()
        run = row["run_accession"].strip()
        if not group or not run or row["compatibility"] != "compatible":
            raise ValueError("Incomplete/incompatible ATAC replicate group")
        eligible_run(manifest, run, require_single=True)
        groups.setdefault(group, []).append(run)
    for group, runs in groups.items():
        if len(runs) != len(set(runs)) or len(runs) < 2:
            raise ValueError(f"Replicate group {group} needs at least two unique eligible runs")
    return groups
