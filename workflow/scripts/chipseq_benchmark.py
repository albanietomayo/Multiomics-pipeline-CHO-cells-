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
from typing import Iterable


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
REFERENCE_IDENTITY = {
    "assembly_name": "CriGri-PICRH-1.0",
    "refseq_accession": "GCF_003668045.3",
    "genbank_accession": "GCA_003668045.2",
    "species": "Cricetulus griseus",
    "taxid": "10029",
    "configured_annotation_release": "104",
    "mitochondrial_accession": "NC_007936.1",
}
EXPECTED_NUCLEAR_SPAN_BP = 2_366_634_374
EXPECTED_MITOCHONDRIAL_LENGTH_BP = 16_284
REFERENCE_SOURCE_SHA256 = {
    "genome_plus_mt.fa.gz": "f5e1effa4d063b9005eeaa0f245e7a09e84d970d6efaeb80f442cfbbf6216923",
    "genome_plus_mt.fa.fai": "2dfb28d82459be6cbddb6e4be79e2c328c436654cc37f2e703c7abb9c50a0d57",
    "reference_provenance.json": "3646c60b547d946814704af383464807dda890cfd371ff78f24edf8307fd582f",
}
REFERENCE_SOURCE_BYTES = {
    "genome_plus_mt.fa.gz": 882_813_762,
    "genome_plus_mt.fa.fai": 24_707,
    "reference_provenance.json": 710,
}
REFERENCE_CONTENT_SHA256 = {
    "nuclear": "5c81f08eafe8f5051704f692cc21a970f58b24b2bc1bedef82cb17ec7a2b4478",
    "mitochondrial": "b7ccf6b1c6981c2b4a0c9e57245bd6712a861e685571db56d274c04fc30ed37c",
    "mapping": "dfd445e136c4bc9c11f9616d251b86732d1aab0cbdbf9d1ababb8093f7430e94",
}
INDEX_COMPONENTS = [
    f"bowtie2_index/genome_plus_mt.{suffix}.bt2l"
    for suffix in ("1", "2", "3", "4", "rev.1", "rev.2")
]
REFERENCE_FILES = [*REFERENCE_SOURCE_SHA256, *INDEX_COMPONENTS]
REFERENCE_INVENTORY = "reference_inventory.json"
REFERENCE_BUILD = {
    "builder": "bowtie2-build",
    "version": "2.5.5",
    "mode": "large-index",
    "input": "TMPDIR decompression of genome_plus_mt.fa.gz",
    "command": (
        "bowtie2-build --large-index --threads 8 "
        "${TMPDIR}/genome_plus_mt.fa ${TMPDIR}/bowtie2_index/genome_plus_mt"
    ),
}


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


def _fai_summary(path: Path, mitochondrial_accession: str) -> dict[str, int]:
    sequence_count = 0
    nuclear_span = 0
    mitochondrial_length = None
    seen: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5 or not fields[1].isdigit() or int(fields[1]) <= 0:
                raise ValueError(f"Malformed FAI row: {path}")
            name, length = fields[0], int(fields[1])
            if not name or name in seen:
                raise ValueError(f"Missing or duplicate FAI sequence: {path}")
            seen.add(name)
            sequence_count += 1
            if name == mitochondrial_accession:
                mitochondrial_length = length
            else:
                nuclear_span += length
    if sequence_count < 2 or mitochondrial_length is None or nuclear_span <= 0:
        raise ValueError("Combined-reference FAI lacks nuclear or mitochondrial span")
    return {
        "sequence_count": sequence_count,
        "nuclear_span_bp": nuclear_span,
        "mitochondrial_length_bp": mitochondrial_length,
    }


def _reference_document(root: Path) -> dict[str, object]:
    paths = {relative: root / relative for relative in REFERENCE_FILES}
    for relative, path in paths.items():
        if not path.is_file() or path.stat().st_size <= 0:
            raise ValueError(f"Missing or empty shared reference resource: {relative}")

    for relative, expected in REFERENCE_SOURCE_BYTES.items():
        if paths[relative].stat().st_size != expected:
            raise ValueError(f"Validated shared reference source size mismatch: {relative}")

    observed_hashes = {relative: sha256(path) for relative, path in paths.items()}
    for relative, expected in REFERENCE_SOURCE_SHA256.items():
        if observed_hashes[relative] != expected:
            raise ValueError(f"Validated shared reference source hash mismatch: {relative}")

    provenance = json.loads(paths["reference_provenance.json"].read_text(encoding="utf-8"))
    provenance_identity = {
        "nuclear_accession": REFERENCE_IDENTITY["refseq_accession"],
        "mitochondrial_accession": REFERENCE_IDENTITY["mitochondrial_accession"],
        "configured_annotation_release": int(
            REFERENCE_IDENTITY["configured_annotation_release"]
        ),
        "annotation_release_independently_verified": False,
    }
    for field, expected in provenance_identity.items():
        if provenance.get(field) != expected:
            raise ValueError(f"Shared ChIP reference provenance mismatch for {field}")
    provenance_hashes = {
        "nuclear_sha256": REFERENCE_CONTENT_SHA256["nuclear"],
        "mitochondrial_sha256": REFERENCE_CONTENT_SHA256["mitochondrial"],
        "mapping_sha256": REFERENCE_CONTENT_SHA256["mapping"],
    }
    for field, expected in provenance_hashes.items():
        if provenance.get(field) != expected:
            raise ValueError(f"Shared ChIP reference provenance hash mismatch for {field}")

    fai = _fai_summary(
        paths["genome_plus_mt.fa.fai"],
        REFERENCE_IDENTITY["mitochondrial_accession"],
    )
    if fai["nuclear_span_bp"] != EXPECTED_NUCLEAR_SPAN_BP:
        raise ValueError("Shared ChIP reference nuclear span disagrees with the verified FAI")
    if fai["mitochondrial_length_bp"] != EXPECTED_MITOCHONDRIAL_LENGTH_BP:
        raise ValueError("Shared ChIP reference mitochondrial length disagrees with the verified FAI")
    if provenance.get("mitochondrial_length") != fai["mitochondrial_length_bp"]:
        raise ValueError("Shared ChIP reference provenance mitochondrial length disagrees with FAI")

    return {
        "schema_version": 2,
        "identity": dict(REFERENCE_IDENTITY),
        "fai": fai,
        "source_content_sha256": dict(REFERENCE_CONTENT_SHA256),
        "files": [
            {
                "relative_path": relative,
                "bytes": paths[relative].stat().st_size,
                "sha256": observed_hashes[relative],
            }
            for relative in REFERENCE_FILES
        ],
        "required_index_components": list(INDEX_COMPONENTS),
        "build": dict(REFERENCE_BUILD),
    }


def write_reference_inventory(root: Path, output: Path | None = None) -> Path:
    """Create the immutable inventory after one successful one-time index build."""
    if not root.is_absolute():
        raise ValueError("Shared reference root must be an explicit absolute path")
    root = root.resolve()
    if not root.is_dir():
        raise ValueError(f"Shared reference root is not a directory: {root}")
    destination = output or root / REFERENCE_INVENTORY
    if destination.exists() or destination.with_name(destination.name + ".part").exists():
        raise ValueError(f"Reference inventory output collision: {destination}")
    atomic_text(
        destination,
        json.dumps(_reference_document(root), indent=2, sort_keys=True) + "\n",
    )
    return destination


def reference_inventory(root: Path) -> dict[str, object]:
    """Validate one explicit, prebuilt and inventory-pinned shared reference."""
    if not root.is_absolute():
        raise ValueError("Shared reference root must be an explicit absolute path")
    root = root.resolve()
    if not root.is_dir():
        raise ValueError(f"Shared reference root is not a directory: {root}")
    inventory_path = root / REFERENCE_INVENTORY
    if not inventory_path.is_file() or inventory_path.stat().st_size <= 0:
        raise ValueError(f"Missing or empty shared reference resource: {REFERENCE_INVENTORY}")
    expected = json.loads(inventory_path.read_text(encoding="utf-8"))
    observed = _reference_document(root)
    if expected != observed:
        raise ValueError("Shared reference differs from its immutable reference inventory")
    result = dict(observed)
    result["root"] = str(root)
    result["inventory_file"] = {
        "relative_path": REFERENCE_INVENTORY,
        "execution_path": str(inventory_path),
        "bytes": inventory_path.stat().st_size,
        "sha256": sha256(inventory_path),
    }
    result["files"] = [
        {**item, "execution_path": str(root / item["relative_path"])}
        for item in observed["files"]
    ]
    return result


def validate_reference_inventory(root: Path, expected_path: Path) -> dict[str, object]:
    expected = json.loads(expected_path.read_text(encoding="utf-8"))
    if isinstance(expected, dict) and "shared_reference" in expected:
        expected = expected["shared_reference"]
    observed = reference_inventory(root)
    if observed != expected:
        raise ValueError("Shared reference differs from the immutable submission inventory")
    return observed


def artifact_record(
    stage: str,
    transient: Iterable[Path],
    persistent: Iterable[Path],
) -> dict[str, object]:
    records = []
    for retention, outputs in (("transient", transient), ("persistent", persistent)):
        for output in outputs:
            if not output.is_file() or output.stat().st_size <= 0:
                raise ValueError(f"Missing or empty stage output: {output}")
            records.append({
                "path": str(output),
                "bytes": output.stat().st_size,
                "sha256": sha256(output),
                "retention": retention,
                "regenerable": retention == "transient",
            })
    if not records:
        raise ValueError("At least one stage output is required")
    return {"schema_version": 2, "stage": stage, "status": "complete", "artifacts": records}


def write_artifact_record(
    stage: str,
    marker: Path,
    transient: Iterable[Path],
    persistent: Iterable[Path],
) -> None:
    atomic_text(
        marker,
        json.dumps(
            artifact_record(stage, transient, persistent),
            indent=2,
            sort_keys=True,
        ) + "\n",
    )


def persist_files(destination: Path, mappings: list[str]) -> list[Path]:
    """Copy an allowlisted set once, verify it, and refuse destination collisions."""
    if not mappings:
        raise ValueError("At least one persistence mapping is required")
    destination.mkdir(parents=True, exist_ok=True)
    persisted = []
    seen: set[Path] = set()
    for mapping in mappings:
        if "=" not in mapping:
            raise ValueError(f"Malformed persistence mapping: {mapping!r}")
        relative_text, source_text = mapping.split("=", 1)
        relative = Path(relative_text)
        if relative.is_absolute() or not relative.parts or ".." in relative.parts:
            raise ValueError(f"Unsafe persistent relative path: {relative_text!r}")
        source = Path(source_text)
        target = destination / relative
        temporary = target.with_name(target.name + ".part")
        if target in seen or target.exists() or temporary.exists():
            raise ValueError(f"Persistent output collision: {target}")
        seen.add(target)
        if not source.is_file() or source.stat().st_size <= 0:
            raise ValueError(f"Missing or empty persistence source: {source}")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, temporary)
        if (
            temporary.stat().st_size != source.stat().st_size
            or sha256(temporary) != sha256(source)
        ):
            raise ValueError(f"Persistent copy validation failed: {source}")
        temporary.replace(target)
        persisted.append(target)
    return persisted

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


def package_provenance(packages: object, package_name: str) -> str:
    if not isinstance(packages, list):
        raise ValueError("Conda package metadata must be a JSON list")
    matches = [
        package for package in packages
        if isinstance(package, dict) and package.get("name") == package_name
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one Conda package named {package_name!r}; "
            f"found {len(matches)}"
        )
    package = matches[0]
    fields = {
        "name": package.get("name"),
        "version": package.get("version"),
        "build": package.get("build_string"),
        "channel": package.get("channel"),
    }
    invalid = [
        field for field, value in fields.items()
        if not isinstance(value, str) or not value.strip()
    ]
    if invalid:
        raise ValueError(
            f"Conda metadata for {package_name!r} has missing or malformed fields: "
            f"{', '.join(invalid)}"
        )
    return "".join(f"{field}: {value}\n" for field, value in fields.items())


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
    package = sub.add_parser("package-provenance")
    package.add_argument("--package", required=True)
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
    reference = sub.add_parser("check-reference")
    reference.add_argument("--shared-reference-root", type=Path, required=True)
    reference.add_argument("--expected", type=Path)
    reference.add_argument("--output", type=Path)
    create_reference = sub.add_parser("create-reference-inventory")
    create_reference.add_argument("--shared-reference-root", type=Path, required=True)
    create_reference.add_argument("--output", type=Path)
    record = sub.add_parser("record-stage")
    record.add_argument("--stage", required=True)
    record.add_argument("--marker", type=Path, required=True)
    record.add_argument("--transient", action="append", type=Path, default=[])
    record.add_argument("--persistent", action="append", type=Path, default=[])
    persist = sub.add_parser("persist-files")
    persist.add_argument("--destination", type=Path, required=True)
    persist.add_argument("mappings", nargs="+")
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
    elif args.action == "package-provenance":
        print(package_provenance(json.load(sys.stdin), args.package), end="")
    elif args.action == "complete-stage":
        write_completion(args.stage, args.marker, args.outputs)
    elif args.action == "verify-stage":
        verify_completion(args.stage, args.marker)
    elif args.action == "make-alignment-plan":
        make_alignment_plan(args.run, args.processed_fastq, args.fastp_json, args.output)
    elif args.action == "make-filtering-plan":
        make_filtering_plan(args.run, args.bam, args.qc, args.config, args.output)
    elif args.action == "check-reference":
        inventory = (
            validate_reference_inventory(args.shared_reference_root, args.expected)
            if args.expected else reference_inventory(args.shared_reference_root)
        )
        content = json.dumps(inventory, indent=2, sort_keys=True) + "\n"
        if args.output:
            atomic_text(args.output, content)
        else:
            print(content, end="")
    elif args.action == "create-reference-inventory":
        output = write_reference_inventory(args.shared_reference_root, args.output)
        print(f"REFERENCE_INVENTORY={output}")
    elif args.action == "record-stage":
        write_artifact_record(args.stage, args.marker, args.transient, args.persistent)
    elif args.action == "persist-files":
        for path in persist_files(args.destination, args.mappings):
            print(path)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, json.JSONDecodeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2) from error
