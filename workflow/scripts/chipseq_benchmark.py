#!/usr/bin/env python3
"""Fail-closed support utilities for the ChIP-seq P10/P50/P90 benchmark."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BENCHMARK_DIR = ROOT / "benchmarks/chipseq/2026-09-18"
AUTHORITATIVE_SOURCE = ROOT / "config/samples.tsv"
SELECTION = BENCHMARK_DIR / "benchmark_selection.tsv"
MANIFEST = BENCHMARK_DIR / "chipseq_benchmark_manifest.tsv"
EXPECTED = {
    "P10": "ERR868176",
    "P50": "ERR868152",
    "P90": "SRR20770294",
}
MANIFEST_FIELDS = [
    "benchmark_class",
    "run_accession",
    "study_accession",
    "library_layout",
    "instrument_platform",
    "filename",
    "fastq_ftp",
    "fastq_bytes",
    "fastq_md5",
    "authoritative_source",
    "authoritative_source_sha256",
    "source_git_commit",
]
SACCT_FIELDS = [
    "JobID",
    "JobName",
    "State",
    "Elapsed",
    "TotalCPU",
    "AllocCPUS",
    "MaxRSS",
    "AveRSS",
    "ExitCode",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(path)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Missing or duplicate TSV columns: {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"TSV contains no rows: {path}")
    return rows


def tsv_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    import io

    output = io.StringIO()
    writer = csv.DictWriter(
        output,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def git_value(root: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), *args], text=True
    ).strip()


def load_selection(path: Path = SELECTION) -> list[dict[str, str]]:
    rows = read_tsv(path)
    if set(rows[0]) != {"benchmark_class", "run_accession"}:
        raise ValueError("Benchmark selection must contain exactly two columns")
    observed = {row["benchmark_class"]: row["run_accession"] for row in rows}
    if len(rows) != 3 or len(observed) != 3 or observed != EXPECTED:
        raise ValueError(
            f"Benchmark selection must be exactly {EXPECTED}; observed {observed}"
        )
    if len({row["run_accession"] for row in rows}) != 3:
        raise ValueError("Benchmark run accessions must be unique")
    return rows


def _single(value: str, label: str, run: str) -> str:
    values = [item.strip() for item in value.split(";") if item.strip()]
    if len(values) != 1:
        raise ValueError(f"Expected one {label} for {run}; found {len(values)}")
    return values[0]


def generate_manifest_rows(
    source: Path = AUTHORITATIVE_SOURCE,
    selection: Path = SELECTION,
    root: Path = ROOT,
    source_commit: str | None = None,
) -> list[dict[str, str]]:
    selected = load_selection(selection)
    source_rows = read_tsv(source)
    required = {
        "run_accession", "study_accession", "library_layout",
        "instrument_platform", "fastq_ftp", "fastq_md5", "fastq_bytes", "omics",
    }
    if not required.issubset(source_rows[0]):
        raise ValueError(f"Authoritative metadata lacks columns: {sorted(required - set(source_rows[0]))}")
    source_digest = sha256(source)
    source_relative = source.relative_to(root)
    commit = source_commit or git_value(
        root, "log", "-1", "--format=%H", "--", str(source_relative)
    )
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError(f"Malformed authoritative source commit: {commit!r}")
    result = []
    for choice in selected:
        run = choice["run_accession"]
        matches = [row for row in source_rows if row["run_accession"] == run]
        if len(matches) != 1:
            raise ValueError(f"Expected one authoritative row for {run}; found {len(matches)}")
        row = matches[0]
        if row["omics"] != "ChIP-seq":
            raise ValueError(f"Benchmark run is not ChIP-seq: {run}")
        if row["library_layout"] != "SINGLE" or row["instrument_platform"] != "ILLUMINA":
            raise ValueError(f"Unsupported benchmark layout/platform: {run}")
        url = _single(row["fastq_ftp"], "FASTQ URL", run)
        md5 = _single(row["fastq_md5"], "FASTQ MD5", run).lower()
        byte_text = _single(row["fastq_bytes"], "FASTQ byte value", run)
        if not re.fullmatch(r"[0-9a-f]{32}", md5):
            raise ValueError(f"Malformed MD5 for {run}")
        if not byte_text.isdigit() or int(byte_text) <= 0:
            raise ValueError(f"Malformed byte value for {run}")
        if not url.startswith("ftp.sra.ebi.ac.uk/"):
            raise ValueError(f"Unexpected authoritative FASTQ host for {run}: {url}")
        filename = url.rsplit("/", 1)[-1]
        if filename != f"{run}.fastq.gz":
            raise ValueError(f"Unexpected SINGLE FASTQ filename for {run}: {filename}")
        result.append({
            "benchmark_class": choice["benchmark_class"],
            "run_accession": run,
            "study_accession": row["study_accession"],
            "library_layout": row["library_layout"],
            "instrument_platform": row["instrument_platform"],
            "filename": filename,
            "fastq_ftp": url,
            "fastq_bytes": byte_text,
            "fastq_md5": md5,
            "authoritative_source": str(source_relative),
            "authoritative_source_sha256": source_digest,
            "source_git_commit": commit,
        })
    return result


def write_manifest(output: Path, source: Path, selection: Path, root: Path) -> None:
    rows = generate_manifest_rows(source, selection, root)
    atomic_text(output, tsv_text(MANIFEST_FIELDS, rows))


def validate_manifest(
    manifest: Path = MANIFEST,
    source: Path = AUTHORITATIVE_SOURCE,
    selection: Path = SELECTION,
    root: Path = ROOT,
) -> list[dict[str, str]]:
    rows = read_tsv(manifest)
    if list(rows[0]) != MANIFEST_FIELDS:
        raise ValueError("Unexpected benchmark manifest columns or order")
    commits = {row["source_git_commit"] for row in rows}
    if len(commits) != 1:
        raise ValueError("Benchmark manifest has inconsistent source commits")
    manifest_commit = next(iter(commits))
    if (root / ".git").exists():
        source_relative = source.relative_to(root)
        observed_commit = git_value(
            root, "log", "-1", "--format=%H", "--", str(source_relative)
        )
        if observed_commit != manifest_commit:
            raise ValueError("Benchmark manifest authoritative-source commit differs")
    expected = generate_manifest_rows(
        source, selection, root, source_commit=manifest_commit
    )
    if rows != expected:
        raise ValueError("Committed benchmark manifest differs from deterministic authoritative derivation")
    return rows


def manifest_item(rows: list[dict[str, str]], benchmark_class: str) -> dict[str, str]:
    if benchmark_class not in EXPECTED:
        raise ValueError(f"Unsafe benchmark class: {benchmark_class!r}")
    matches = [row for row in rows if row["benchmark_class"] == benchmark_class]
    if len(matches) != 1 or matches[0]["run_accession"] != EXPECTED[benchmark_class]:
        raise ValueError(f"Missing, duplicate, or mismatched benchmark class: {benchmark_class}")
    return matches[0]


def existing_ancestor(path: Path) -> Path:
    candidate = path.resolve()
    while not candidate.exists():
        if candidate.parent == candidate:
            raise ValueError(f"No existing ancestor for {path}")
        candidate = candidate.parent
    if not candidate.is_dir():
        candidate = candidate.parent
    return candidate


def available_bytes(path: Path) -> int:
    return shutil.disk_usage(existing_ancestor(path)).free


def storage_preflight(
    row: dict[str, str], output_root: Path, scratch_root: Path
) -> dict[str, object]:
    run_dir = output_root / f"{row['benchmark_class']}_{row['run_accession']}"
    if run_dir.exists():
        raise ValueError(f"Benchmark output collision: {run_dir}")
    known = int(row["fastq_bytes"])
    output_parent = existing_ancestor(output_root)
    scratch_parent = existing_ancestor(scratch_root)
    for label, parent in (("output", output_parent), ("scratch", scratch_parent)):
        if not os.access(parent, os.W_OK | os.X_OK):
            raise ValueError(f"Required {label} directory is not writable: {parent}")
    output_free = shutil.disk_usage(output_parent).free
    scratch_free = shutil.disk_usage(scratch_parent).free
    if output_free < known:
        raise ValueError("Output storage cannot hold even the known compressed input")
    if scratch_free < known:
        raise ValueError("Scratch cannot hold even the known compressed input")
    return {
        "benchmark_class": row["benchmark_class"],
        "run_accession": row["run_accession"],
        "known_compressed_input_bytes": known,
        "output_available_bytes": output_free,
        "scratch_available_bytes": scratch_free,
        "output_existing_ancestor": str(output_parent),
        "scratch_existing_ancestor": str(scratch_parent),
        "downstream_expansion_bytes": "unknown_pending_empirical_benchmark",
        "output_directory": str(run_dir),
        "collision": False,
        "minimum_known_storage_satisfied": True,
    }


def parse_sacct(input_path: Path, output_path: Path) -> None:
    lines = [line for line in input_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise ValueError("sacct input is empty")
    reader = csv.DictReader(lines, delimiter="|")
    if reader.fieldnames is None:
        raise ValueError("sacct input has no header")
    missing_columns = [field for field in SACCT_FIELDS if field not in reader.fieldnames]
    if missing_columns:
        raise ValueError(f"sacct input lacks requested columns: {missing_columns}")
    rows = []
    for raw in reader:
        if raw.get(None):
            raise ValueError("sacct row has more fields than its header")
        rows.append({field: (raw.get(field, "").strip() or "not_available") for field in SACCT_FIELDS})
    if not rows:
        raise ValueError("sacct input has no accounting rows")
    atomic_text(output_path, tsv_text(SACCT_FIELDS, rows))


def completion_record(stage: str, outputs: list[Path]) -> dict[str, object]:
    if not outputs:
        raise ValueError("At least one stage output is required")
    records = []
    for output in outputs:
        if not output.is_file() or output.stat().st_size <= 0:
            raise ValueError(f"Missing or empty stage output: {output}")
        records.append({"path": str(output), "bytes": output.stat().st_size, "sha256": sha256(output)})
    return {"schema_version": 1, "stage": stage, "status": "complete", "outputs": records}


def write_completion(stage: str, marker: Path, outputs: list[Path]) -> None:
    atomic_text(marker, json.dumps(completion_record(stage, outputs), indent=2, sort_keys=True) + "\n")


def verify_completion(stage: str, marker: Path) -> None:
    obj = json.loads(marker.read_text(encoding="utf-8"))
    if obj.get("schema_version") != 1 or obj.get("stage") != stage or obj.get("status") != "complete":
        raise ValueError(f"Invalid completion marker: {marker}")
    outputs = obj.get("outputs")
    if not isinstance(outputs, list) or not outputs:
        raise ValueError(f"Completion marker contains no outputs: {marker}")
    for item in outputs:
        path = Path(item["path"])
        if not path.is_file() or path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
            raise ValueError(f"Completed stage output changed: {path}")


def make_alignment_plan(run: str, processed_fastq: Path, fastp_json: Path, output: Path) -> None:
    if run not in EXPECTED.values():
        raise ValueError(f"Unsafe benchmark accession: {run}")
    report = json.loads(fastp_json.read_text(encoding="utf-8"))
    reads = report.get("summary", {}).get("after_filtering", {}).get("total_reads")
    if not isinstance(reads, int) or reads <= 0:
        raise ValueError("fastp report has no positive after-filtering read count")
    obj = {
        "schema_version": 1,
        "library_layout": "SINGLE",
        "instrument_platform": "ILLUMINA",
        "runs": [{
            "run_accession": run,
            "role": "ip",
            "destination": f"inputs/{run}.fastq.gz",
            "library_layout": "SINGLE",
            "instrument_platform": "ILLUMINA",
            "reads": reads,
            "processed_fastq": str(processed_fastq),
            "processed_fastq_sha256": sha256(processed_fastq),
        }],
    }
    atomic_text(output, json.dumps(obj, indent=2, sort_keys=True) + "\n")


def make_filtering_plan(run: str, bam: Path, qc_path: Path, config: Path, output: Path) -> None:
    if run not in EXPECTED.values():
        raise ValueError(f"Unsafe benchmark accession: {run}")
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    if qc.get("run_accession") != run:
        raise ValueError("Alignment QC accession mismatch")
    obj = {
        "schema_version": 1,
        "library_layout": "SINGLE",
        "instrument_platform": "ILLUMINA",
        "source_root": ".",
        "configuration_sha256": sha256(config),
        "upstream_raw_bam_cleanup_authorized": False,
        "runs": [{
            "run_accession": run,
            "role": "ip",
            "source_relative": str(bam),
            "sha256": sha256(bam),
            "bytes": bam.stat().st_size,
            "input_reads": qc["primary_reads"],
            "mapped_reads": qc["mapped_reads"],
            "nuclear_mapq_ge_threshold": qc["nuclear_mapq_ge_threshold"],
        }],
    }
    atomic_text(output, json.dumps(obj, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    generate = sub.add_parser("generate-manifest")
    generate.add_argument("--source", type=Path, default=AUTHORITATIVE_SOURCE)
    generate.add_argument("--selection", type=Path, default=SELECTION)
    generate.add_argument("--output", type=Path, default=MANIFEST)
    check = sub.add_parser("check-manifest")
    check.add_argument("--manifest", type=Path, default=MANIFEST)
    preflight = sub.add_parser("preflight")
    preflight.add_argument("--benchmark-class", required=True)
    preflight.add_argument("--manifest", type=Path, default=MANIFEST)
    preflight.add_argument("--output-root", type=Path, required=True)
    preflight.add_argument("--scratch-root", type=Path, required=True)
    sacct = sub.add_parser("parse-sacct")
    sacct.add_argument("--input", type=Path, required=True)
    sacct.add_argument("--output", type=Path, required=True)
    complete = sub.add_parser("complete-stage")
    complete.add_argument("--stage", required=True)
    complete.add_argument("--marker", type=Path, required=True)
    complete.add_argument("outputs", nargs="+", type=Path)
    verify = sub.add_parser("verify-stage")
    verify.add_argument("--stage", required=True)
    verify.add_argument("--marker", type=Path, required=True)
    align = sub.add_parser("make-alignment-plan")
    align.add_argument("--run", required=True)
    align.add_argument("--processed-fastq", type=Path, required=True)
    align.add_argument("--fastp-json", type=Path, required=True)
    align.add_argument("--output", type=Path, required=True)
    filtering = sub.add_parser("make-filtering-plan")
    filtering.add_argument("--run", required=True)
    filtering.add_argument("--bam", type=Path, required=True)
    filtering.add_argument("--qc", type=Path, required=True)
    filtering.add_argument("--config", type=Path, default=Path("config/chipseq_filtering.json"))
    filtering.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.action == "generate-manifest":
        write_manifest(args.output, args.source, args.selection, ROOT)
    elif args.action == "check-manifest":
        rows = validate_manifest(args.manifest)
        print(f"CHIPSEQ_BENCHMARK_MANIFEST=PASS runs={len(rows)}")
    elif args.action == "preflight":
        row = manifest_item(validate_manifest(args.manifest), args.benchmark_class)
        print(json.dumps(storage_preflight(row, args.output_root, args.scratch_root), indent=2, sort_keys=True))
    elif args.action == "parse-sacct":
        parse_sacct(args.input, args.output)
    elif args.action == "complete-stage":
        write_completion(args.stage, args.marker, args.outputs)
    elif args.action == "verify-stage":
        verify_completion(args.stage, args.marker)
    elif args.action == "make-alignment-plan":
        make_alignment_plan(args.run, args.processed_fastq, args.fastp_json, args.output)
    elif args.action == "make-filtering-plan":
        make_filtering_plan(args.run, args.bam, args.qc, args.config, args.output)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, json.JSONDecodeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2) from error
