#!/usr/bin/env python3
"""Fail-closed support for one-at-a-time, scratch-resident ChIP production."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import os
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PEAK_PLAN = ROOT / "config/chipseq_peak_calling_plan.tsv"
ANALYSIS_PLAN = ROOT / "snapshots/chipseq/incremental_planning_validation_001/chipseq_analysis_plan.tsv"
PROCESSING_PLAN = ROOT / "snapshots/chipseq/runtime_planning_validation_001/chipseq_processing_runs.tsv"
SAMPLES = ROOT / "config/samples.tsv"
PEAK_SUMMARY = ROOT / "config/chipseq_peak_calling_plan_summary.json"
PEAK_POLICY = ROOT / "config/chipseq_peak_calling_policy.json"
FILTER_CONFIG = ROOT / "config/chipseq_filtering.json"
MTDNA = "NC_007936.1"
ACCESSION = re.compile(r"[DES]RR[0-9]+")
MANIFEST_FIELDS = [
    "analysis_id", "study", "ip", "control", "target", "peak_mode",
    "fragment_size_policy", "fragment_size_bp", "start_utc", "end_utc",
    "slurm_job", "final_status", "persistent_artifact", "sha256", "size_bytes",
]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


DYNAMIC = load_module("chipseq_peak_dynamic_production", ROOT / "workflow/scripts/chipseq_peak_calling_dynamic.py")
BENCHMARK = load_module("chipseq_benchmark_production", ROOT / "workflow/scripts/chipseq_benchmark.py")


def sha256(path: Path | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: Path | str, text: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    if path.exists() or temporary.exists():
        raise ValueError(f"Refusing output collision: {path}")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def atomic_json(path: Path | str, value: object) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def table(path: Path | str) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Missing or duplicate TSV columns: {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty TSV: {path}")
    return rows


def one(rows: list[dict[str, str]], key: str, value: str, label: str) -> dict[str, str]:
    matches = [row for row in rows if row.get(key) == value]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {label} for {value}; found {len(matches)}")
    return matches[0]


def _single(value: str, label: str, run: str) -> str:
    values = [item.strip() for item in value.split(";") if item.strip()]
    if len(values) != 1:
        raise ValueError(f"Expected one {label} for {run}; found {len(values)}")
    return values[0]


def validated_plans() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    peaks = DYNAMIC.read_peak_plan(PEAK_PLAN)
    ready = [row for row in table(ANALYSIS_PLAN) if row.get("analysis_status") == "ready"]
    peak_pairs = {
        row["analysis_id"]: (row["ip_run_accession"], row["control_run_accession"])
        for row in peaks
    }
    ready_pairs = {
        row["analysis_id"]: (row["ip_run_accession"], row["control_run_accession"])
        for row in ready
    }
    if len(peaks) != 18 or len(ready) != 18 or peak_pairs != ready_pairs:
        raise ValueError("Peak plan does not exactly inherit all authoritative ready IP/Input pairs")
    summary = json.loads(PEAK_SUMMARY.read_text(encoding="utf-8"))
    hashes = summary.get("input_sha256", {})
    DYNAMIC.verify_hash(ANALYSIS_PLAN, hashes.get("analysis_plan", ""), "analysis plan")
    DYNAMIC.verify_hash(PEAK_POLICY, hashes.get("policy", ""), "peak policy")
    DYNAMIC.verify_hash(PEAK_PLAN, summary.get("peak_calling_plan_sha256", ""), "peak plan")
    return peaks, ready


def sample(run: str, expected_study: str) -> dict[str, object]:
    if not ACCESSION.fullmatch(run):
        raise ValueError(f"Unsafe run accession: {run!r}")
    row = one(table(SAMPLES), "run_accession", run, "authoritative sample row")
    if (row.get("study_accession") != expected_study or row.get("omics") != "ChIP-seq"
            or row.get("library_layout") != "SINGLE"
            or row.get("instrument_platform") != "ILLUMINA"):
        raise ValueError(f"Unsupported or mismatched authoritative metadata: {run}")
    url = _single(row.get("fastq_ftp", ""), "FASTQ URL", run)
    md5 = _single(row.get("fastq_md5", ""), "FASTQ MD5", run).lower()
    size = _single(row.get("fastq_bytes", ""), "FASTQ size", run)
    if (not url.startswith("ftp.sra.ebi.ac.uk/") or not re.fullmatch(r"[0-9a-f]{32}", md5)
            or not size.isdigit() or int(size) <= 0
            or url.rsplit("/", 1)[-1] != f"{run}.fastq.gz"):
        raise ValueError(f"Malformed ENA FASTQ evidence: {run}")
    return {"run_accession": run, "filename": f"{run}.fastq.gz", "fastq_ftp": url,
            "fastq_md5": md5, "fastq_bytes": int(size)}


def select_analysis(analysis_id: str, expected_control: str | None = None) -> dict[str, object]:
    peaks, ready = validated_plans()
    row = one(peaks, "analysis_id", analysis_id, "eligible peak analysis")
    authoritative = one(ready, "analysis_id", analysis_id, "authoritative ready analysis")
    if expected_control and row["control_run_accession"] != expected_control:
        raise ValueError("Explicit control does not match the authoritative pairing")
    if any(row[key] != authoritative[key] for key in
           ("analysis_id", "study_accession", "ip_run_accession", "control_run_accession", "declared_target")):
        raise ValueError("Peak row differs from the authoritative biological pairing layer")
    processing = table(PROCESSING_PLAN)
    for run, role in ((row["ip_run_accession"], "ip"), (row["control_run_accession"], "input")):
        item = one(processing, "run_accession", run, "processing-plan run")
        if item.get("library_role") != role or item.get("study_accession") != row["study_accession"]:
            raise ValueError(f"Processing role/study mismatch: {run}")
    return {"analysis": row, "samples": {
        "ip": sample(row["ip_run_accession"], row["study_accession"]),
        "input": sample(row["control_run_accession"], row["study_accession"]),
    }}


def preflight(analysis_id: str, control: str | None, reference_root: Path,
              scratch_root: Path, production_root: Path) -> dict[str, object]:
    selected = select_analysis(analysis_id, control)
    if not scratch_root.is_absolute() or not scratch_root.is_dir():
        raise ValueError("TMPDIR/scratch must be an existing absolute directory")
    probe = scratch_root / f".chipseq-production-write-probe-{os.getpid()}"
    try:
        probe.mkdir()
        probe.rmdir()
    except OSError as exc:
        raise ValueError(f"TMPDIR/scratch is not writable: {scratch_root}") from exc
    sizes = [item["fastq_bytes"] for item in selected["samples"].values()]
    required = max(8 * 1024**3, 4 * max(sizes) + 2 * 1024**3)
    free = shutil.disk_usage(scratch_root).free
    if free < required:
        raise ValueError(f"Scratch free bytes {free} are below structural reserve {required}")
    reference = BENCHMARK.reference_inventory(reference_root)
    analysis_dir = production_root / "analyses" / analysis_id
    if analysis_dir.exists():
        raise ValueError(f"Production analysis output collision: {analysis_dir}")
    result = dict(selected)
    result.update({
        "schema_version": 1,
        "scope": "single_analysis_lean_chipseq_production",
        "scratch_root": str(scratch_root.resolve()),
        "scratch_free_bytes": free,
        "scratch_structural_reserve_bytes": required,
        "production_root": str(production_root.resolve()),
        "reference": reference,
        "input_sha256": {
            "analysis_plan": sha256(ANALYSIS_PLAN), "peak_plan": sha256(PEAK_PLAN),
            "peak_plan_summary": sha256(PEAK_SUMMARY), "peak_policy": sha256(PEAK_POLICY),
            "processing_plan": sha256(PROCESSING_PLAN), "samples": sha256(SAMPLES),
            "filtering_policy": sha256(FILTER_CONFIG),
        },
    })
    return result


def fastq_manifest(plan: dict[str, object], output: Path) -> None:
    fields = ["run_accession", "filename", "fastq_ftp", "fastq_bytes", "fastq_md5"]
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for role in ("input", "ip"):
        row = plan["samples"][role]
        writer.writerow({key: row[key] for key in fields})
    atomic_text(output, buffer.getvalue())


def alignment_plan(plan: dict[str, object], run: str, role: str, fastq: Path,
                   fastp_json: Path, output: Path) -> None:
    expected = plan["analysis"]["ip_run_accession" if role == "ip" else "control_run_accession"]
    if run != expected or role not in {"ip", "input"}:
        raise ValueError("Run/role differs from selected authoritative analysis")
    report = json.loads(fastp_json.read_text(encoding="utf-8"))
    reads = report.get("summary", {}).get("after_filtering", {}).get("total_reads")
    if type(reads) is not int or reads <= 0 or not fastq.is_file() or fastq.stat().st_size <= 0:
        raise ValueError("Processed FASTQ or fastp read count is invalid")
    atomic_json(output, {"schema_version": 1, "library_layout": "SINGLE",
        "instrument_platform": "ILLUMINA", "runs": [{"run_accession": run, "role": role,
        "destination": f"inputs/{run}.fastq.gz", "library_layout": "SINGLE",
        "instrument_platform": "ILLUMINA", "reads": reads,
        "processed_fastq": str(fastq), "processed_fastq_sha256": sha256(fastq)}]})


def filtering_plan(plan: dict[str, object], run: str, role: str, bam: Path,
                   qc_path: Path, output: Path) -> None:
    expected = plan["analysis"]["ip_run_accession" if role == "ip" else "control_run_accession"]
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    if run != expected or qc.get("run_accession") != run or qc.get("role") != role:
        raise ValueError("Alignment QC differs from selected authoritative run/role")
    atomic_json(output, {"schema_version": 1, "library_layout": "SINGLE",
        "instrument_platform": "ILLUMINA", "source_root": str(bam.parent),
        "configuration_sha256": sha256(FILTER_CONFIG), "upstream_raw_bam_cleanup_authorized": True,
        "runs": [{"run_accession": run, "role": role, "source_relative": bam.name,
        "sha256": sha256(bam), "bytes": bam.stat().st_size, "input_reads": qc["primary_reads"],
        "mapped_reads": qc["mapped_reads"],
        "nuclear_mapq_ge_threshold": qc["nuclear_mapq_ge_threshold"]}]})


def _validated_filtered(run: str, role: str, directory: Path) -> dict[str, object]:
    bam, csi = directory / "filtered.bam", directory / "filtered.bam.csi"
    validation_path, qc_path = directory / "filtered_validation.json", directory / "filtering_qc.json"
    for path in (bam, csi, validation_path, qc_path):
        if not path.is_file() or path.stat().st_size <= 0:
            raise ValueError(f"Missing filtered-run artifact: {path}")
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    if (validation.get("run_accession") != run or validation.get("full_read_verified") is not True
            or qc.get("run_accession") != run or qc.get("role") != role
            or validation.get("retained_reads") != qc.get("retained_reads")
            or validation.get("bam_sha256") != sha256(bam)
            or validation.get("index_sha256") != sha256(csi)):
        raise ValueError(f"Filtered-run validation mismatch: {run}")
    return {"run_accession": run, "role": role, "bam_path": str(bam.resolve()),
        "csi_path": str(csi.resolve()), "bam_sha256": validation["bam_sha256"],
        "csi_sha256": validation["index_sha256"], "retained_reads": validation["retained_reads"],
        "validation_sha256": sha256(validation_path), "qc_sha256": sha256(qc_path)}


def persist_control(plan: dict[str, object], source: Path, destination: Path) -> None:
    run = plan["analysis"]["control_run_accession"]
    record = _validated_filtered(run, "input", source)
    if destination.exists() or destination.with_name(destination.name + ".part").exists():
        raise ValueError(f"Control destination collision: {destination}")
    temporary = destination.with_name(destination.name + ".part")
    temporary.mkdir(parents=True)
    try:
        for name in ("filtered.bam", "filtered.bam.csi", "filtered_validation.json", "filtering_qc.json",
                     "filtering_flow.tsv", "picard_metrics.txt", "alignment_qc.json", "fastp.json", "fastp.html"):
            source_path = source / name
            if source_path.is_file():
                shutil.copy2(source_path, temporary / name)
        copied = _validated_filtered(run, "input", temporary)
        atomic_json(temporary / "control_manifest.json", {"schema_version": 1,
            "status": "validated_temporary_control", "analysis_plan_sha256": plan["input_sha256"]["analysis_plan"],
            "peak_plan_sha256": plan["input_sha256"]["peak_plan"], "run": copied})
        temporary.replace(destination)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def load_control(plan: dict[str, object], destination: Path) -> dict[str, object]:
    run = plan["analysis"]["control_run_accession"]
    record = _validated_filtered(run, "input", destination)
    manifest = json.loads((destination / "control_manifest.json").read_text(encoding="utf-8"))
    if (manifest.get("status") != "validated_temporary_control"
            or manifest.get("analysis_plan_sha256") != plan["input_sha256"]["analysis_plan"]
            or manifest.get("peak_plan_sha256") != plan["input_sha256"]["peak_plan"]
            or manifest.get("run", {}).get("bam_sha256") != record["bam_sha256"]):
        raise ValueError("Persisted control manifest is invalid or stale")
    return record


def make_runtime(plan: dict[str, object], ip_dir: Path, control_dir: Path, output: Path) -> None:
    row = plan["analysis"]
    ip = _validated_filtered(row["ip_run_accession"], "ip", ip_dir)
    control = load_control(plan, control_dir)
    reference_root = Path(plan["reference"]["root"])
    fai = reference_root / "genome_plus_mt.fa.fai"
    provenance_path = reference_root / "reference_provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    reference = DYNAMIC.validate_reference(fai, provenance, sha256(fai))
    reference["source_path"] = str(fai.resolve())
    reference["provenance_sha256"] = sha256(provenance_path)
    runtime = {"schema_version": 1, "scope": "dynamic_chipseq_peak_calling_runtime",
        "analysis_count": 1, "run_count": 2, "analyses": [row],
        "runs": {row["ip_run_accession"]: ip, row["control_run_accession"]: control},
        "reference": reference, "input_sha256": dict(plan["input_sha256"]),
        "pairing_source": str(ANALYSIS_PLAN),
        "pairing_policy": "exact_inheritance_no_downstream_repairing",
        "full_bam_and_csi_hashes_verified": True}
    atomic_json(output, runtime)


def compression_record(inputs: list[Path], outputs: list[Path], destination: Path) -> None:
    if len(inputs) != len(outputs) or not inputs:
        raise ValueError("Signal compression requires paired input/output paths")
    rows = []
    for source, compressed in zip(inputs, outputs):
        if not source.is_file() or not compressed.is_file():
            raise ValueError("Missing signal compression input/output")
        rows.append({"source_name": source.name, "source_bytes": source.stat().st_size,
                     "source_sha256": sha256(source), "gzip_name": compressed.name,
                     "gzip_bytes": compressed.stat().st_size, "gzip_sha256": sha256(compressed)})
    atomic_json(destination, {"schema_version": 1, "method": "gzip_lossless_no_timestamp", "files": rows})


def fragment_record(plan: dict[str, object], parameters_path: Path, output: Path) -> None:
    row = plan["analysis"]
    parameters = json.loads(parameters_path.read_text(encoding="utf-8"))
    if parameters.get("analysis_id") != row["analysis_id"]:
        raise ValueError("Fragment record parameter/analysis mismatch")
    atomic_json(output, {"schema_version": 1, "analysis_id": row["analysis_id"],
        "study_accession": row["study_accession"], "fragment_size_policy": row["fragment_size_policy"],
        "fragment_size_bp": parameters.get("fragment_size_bp"),
        "fragment_size_source": parameters.get("fragment_size_source"),
        "parameters_sha256": sha256(parameters_path)})


def finalize(plan: dict[str, object], artifact_root: Path, output: Path,
             start_utc: str, slurm_job: str) -> None:
    row = plan["analysis"]
    primary = artifact_root / "peaks" / f"{row['analysis_id']}_peaks.{ 'narrowPeak' if row['peak_mode'] == 'narrow' else 'broadPeak'}"
    required_files = (primary, artifact_root / "peak_qc.json", artifact_root / "provenance.json",
                      artifact_root / "parameters.json", artifact_root / "macs3.log")
    for required in required_files:
        if not required.is_file():
            raise ValueError(f"Missing compact final analysis output: {required}")
    for required in required_files[1:4]:
        if required.stat().st_size <= 0:
            raise ValueError(f"Empty compact final analysis output: {required}")
    parameters = json.loads((artifact_root / "parameters.json").read_text(encoding="utf-8"))
    fragment_size = parameters.get("fragment_size_bp")
    if type(fragment_size) is not int or fragment_size <= 0:
        raise ValueError("Final parameters lack a positive fragment size")
    files = sorted(path for path in artifact_root.rglob("*") if path.is_file())
    end = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for path in files:
        writer.writerow({"analysis_id": row["analysis_id"], "study": row["study_accession"],
            "ip": row["ip_run_accession"], "control": row["control_run_accession"],
            "target": row["declared_target"], "peak_mode": row["peak_mode"],
            "fragment_size_policy": row["fragment_size_policy"],
            "fragment_size_bp": fragment_size,
            "start_utc": start_utc, "end_utc": end, "slurm_job": slurm_job,
            "final_status": "completed", "persistent_artifact": path.relative_to(output).as_posix(),
            "sha256": sha256(path), "size_bytes": path.stat().st_size})
    manifest = output / "production_run_manifest.tsv"
    atomic_text(manifest, buffer.getvalue())
    atomic_json(output / "completed.json", {"schema_version": 1, "status": "completed",
        "analysis_id": row["analysis_id"], "manifest_sha256": sha256(manifest), "end_utc": end})


def verify_analysis(root: Path, analysis_id: str) -> None:
    directory = root / "analyses" / analysis_id
    if not (directory / "completed.json").is_file() or not (directory / "production_run_manifest.tsv").is_file():
        raise ValueError(f"Expected analysis is not complete: {analysis_id}")
    marker = json.loads((directory / "completed.json").read_text(encoding="utf-8"))
    manifest = directory / "production_run_manifest.tsv"
    if marker.get("status") != "completed" or marker.get("analysis_id") != analysis_id \
            or marker.get("manifest_sha256") != sha256(manifest):
        raise ValueError(f"Invalid completion record: {analysis_id}")
    rows = table(manifest)
    if not rows or any(row["analysis_id"] != analysis_id or row["final_status"] != "completed" for row in rows):
        raise ValueError(f"Invalid production manifest identity/status: {analysis_id}")
    for row in rows:
        path = directory / row["persistent_artifact"]
        if (not path.is_file() or path.stat().st_size != int(row["size_bytes"])
                or sha256(path) != row["sha256"]):
            raise ValueError(f"Persistent artifact validation failed: {path}")


def control_cleanup(root: Path, control: str, remove: bool) -> None:
    peaks, _ = validated_plans()
    expected = sorted(row["analysis_id"] for row in peaks if row["control_run_accession"] == control)
    if not expected:
        raise ValueError(f"Control is not referenced by an eligible analysis: {control}")
    for analysis_id in expected:
        verify_analysis(root, analysis_id)
    control_dir = root / "shared_controls" / control
    plan = select_analysis(expected[0], control)
    synthetic = {"analysis": plan["analysis"], "input_sha256": {
        "analysis_plan": sha256(ANALYSIS_PLAN), "peak_plan": sha256(PEAK_PLAN)}}
    load_control(synthetic, control_dir)
    print(f"SAFE_TO_REMOVE_CONTROL_BAM={control}")
    print(f"COMPLETED_EXPECTED_ANALYSES={','.join(expected)}")
    if remove:
        for name in ("filtered.bam", "filtered.bam.csi"):
            (control_dir / name).unlink()
        atomic_json(control_dir / "control_bam_removed.json", {"schema_version": 1,
            "control": control, "verified_analyses": expected,
            "removed_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")})
        print(f"CONTROL_BAM_REMOVED={control}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("preflight")
    p.add_argument("--analysis", required=True); p.add_argument("--control")
    p.add_argument("--reference-root", type=Path, required=True)
    p.add_argument("--scratch-root", type=Path, required=True)
    p.add_argument("--production-root", type=Path, required=True); p.add_argument("--output", type=Path)
    p = sub.add_parser("fastq-manifest"); p.add_argument("--plan", type=Path, required=True); p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("alignment-plan"); p.add_argument("--plan", type=Path, required=True); p.add_argument("--run", required=True); p.add_argument("--role", choices=("ip", "input"), required=True); p.add_argument("--fastq", type=Path, required=True); p.add_argument("--fastp-json", type=Path, required=True); p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("filtering-plan"); p.add_argument("--plan", type=Path, required=True); p.add_argument("--run", required=True); p.add_argument("--role", choices=("ip", "input"), required=True); p.add_argument("--bam", type=Path, required=True); p.add_argument("--qc", type=Path, required=True); p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("persist-control"); p.add_argument("--plan", type=Path, required=True); p.add_argument("--source", type=Path, required=True); p.add_argument("--destination", type=Path, required=True)
    p = sub.add_parser("load-control"); p.add_argument("--plan", type=Path, required=True); p.add_argument("--destination", type=Path, required=True)
    p = sub.add_parser("make-runtime"); p.add_argument("--plan", type=Path, required=True); p.add_argument("--ip-dir", type=Path, required=True); p.add_argument("--control-dir", type=Path, required=True); p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("compression-record"); p.add_argument("--input", action="append", type=Path, required=True); p.add_argument("--gzip", action="append", type=Path, required=True); p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("fragment-record"); p.add_argument("--plan", type=Path, required=True); p.add_argument("--parameters", type=Path, required=True); p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("finalize"); p.add_argument("--plan", type=Path, required=True); p.add_argument("--artifacts", type=Path, required=True); p.add_argument("--output", type=Path, required=True); p.add_argument("--start-utc", required=True); p.add_argument("--slurm-job", required=True)
    p = sub.add_parser("control-cleanup"); p.add_argument("--production-root", type=Path, required=True); p.add_argument("--control", required=True); p.add_argument("--remove", action="store_true")
    args = parser.parse_args()
    if args.action == "preflight":
        value = preflight(args.analysis, args.control, args.reference_root, args.scratch_root, args.production_root)
        text = json.dumps(value, indent=2, sort_keys=True) + "\n"
        atomic_text(args.output, text) if args.output else print(text, end="")
    elif args.action == "fastq-manifest": fastq_manifest(json.loads(args.plan.read_text()), args.output)
    elif args.action == "alignment-plan": alignment_plan(json.loads(args.plan.read_text()), args.run, args.role, args.fastq, args.fastp_json, args.output)
    elif args.action == "filtering-plan": filtering_plan(json.loads(args.plan.read_text()), args.run, args.role, args.bam, args.qc, args.output)
    elif args.action == "persist-control": persist_control(json.loads(args.plan.read_text()), args.source, args.destination)
    elif args.action == "load-control": print(json.dumps(load_control(json.loads(args.plan.read_text()), args.destination), indent=2, sort_keys=True))
    elif args.action == "make-runtime": make_runtime(json.loads(args.plan.read_text()), args.ip_dir, args.control_dir, args.output)
    elif args.action == "compression-record": compression_record(args.input, args.gzip, args.output)
    elif args.action == "fragment-record": fragment_record(json.loads(args.plan.read_text()), args.parameters, args.output)
    elif args.action == "finalize": finalize(json.loads(args.plan.read_text()), args.artifacts, args.output, args.start_utc, args.slurm_job)
    elif args.action == "control-cleanup": control_cleanup(args.production_root, args.control, args.remove)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        raise SystemExit(2) from error
