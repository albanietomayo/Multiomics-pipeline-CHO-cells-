#!/usr/bin/env python3
"""Fail-closed support for one-time shared ChIP reference provisioning."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SUPPORT_PATH = ROOT / "workflow/scripts/chipseq_benchmark.py"
SPEC = importlib.util.spec_from_file_location("chipseq_benchmark_support", SUPPORT_PATH)
SUPPORT = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(SUPPORT)

EXPECTED_BRANCH = "chipseq-incremental-generalization"
BASELINE_COMMIT = "e66a41a1204196cc174472290fef50cb844de29c"
SOURCE_ROOT = Path(
    "/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_chipseq/results/chipseq/"
    "alignment/slurm/submission__niyvyr2/job_10297460/outputs/reference"
)
CANONICAL_ROOT = Path(
    "/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_benchmark/resources/"
    "reference/CriGri-PICRH-1.0/chipseq"
)
QUOTA_ROOT = Path("/cephyr/users/mayoa/Vera")
QUOTA_BYTES = 30 * 1024**3
P10_RESERVE_BYTES = 1_234_803_098
EXPECTED_MAPPING_SHA256 = SUPPORT.REFERENCE_CONTENT_SHA256["mapping"]
EXPECTED_SEQUENCE_COUNT = 648
EXPECTED_TOTAL_BP = 2_366_650_658


def hash_decompressed(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_source(source_root: Path, verify_decompressed: bool = True) -> dict[str, object]:
    if not source_root.is_absolute() or not source_root.is_dir():
        raise ValueError(f"Historical source root is not an absolute directory: {source_root}")
    files = []
    for name in SUPPORT.REFERENCE_SOURCE_SHA256:
        path = source_root / name
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"Historical source is missing or not a regular file: {path}")
        observed_bytes = path.stat().st_size
        if observed_bytes != SUPPORT.REFERENCE_SOURCE_BYTES[name]:
            raise ValueError(f"Historical source size mismatch: {name}")
        observed_hash = SUPPORT.sha256(path)
        if observed_hash != SUPPORT.REFERENCE_SOURCE_SHA256[name]:
            raise ValueError(f"Historical source SHA256 mismatch: {name}")
        files.append({"relative_path": name, "bytes": observed_bytes, "sha256": observed_hash})

    fai = SUPPORT._fai_summary(
        source_root / "genome_plus_mt.fa.fai",
        SUPPORT.REFERENCE_IDENTITY["mitochondrial_accession"],
    )
    if fai != {
        "sequence_count": EXPECTED_SEQUENCE_COUNT,
        "nuclear_span_bp": SUPPORT.EXPECTED_NUCLEAR_SPAN_BP,
        "mitochondrial_length_bp": SUPPORT.EXPECTED_MITOCHONDRIAL_LENGTH_BP,
    }:
        raise ValueError(f"Historical FAI identity mismatch: {fai}")
    if sum((fai["nuclear_span_bp"], fai["mitochondrial_length_bp"])) != EXPECTED_TOTAL_BP:
        raise ValueError("Historical FAI total span mismatch")

    provenance = json.loads(
        (source_root / "reference_provenance.json").read_text(encoding="utf-8")
    )
    expected_provenance = {
        "nuclear_accession": SUPPORT.REFERENCE_IDENTITY["refseq_accession"],
        "mitochondrial_accession": SUPPORT.REFERENCE_IDENTITY["mitochondrial_accession"],
        "configured_annotation_release": int(
            SUPPORT.REFERENCE_IDENTITY["configured_annotation_release"]
        ),
        "annotation_release_independently_verified": False,
        "mapping_sha256": EXPECTED_MAPPING_SHA256,
        "mitochondrial_length": SUPPORT.EXPECTED_MITOCHONDRIAL_LENGTH_BP,
    }
    for field, expected in expected_provenance.items():
        if provenance.get(field) != expected:
            raise ValueError(f"Historical provenance mismatch for {field}")

    decompressed_sha256 = None
    if verify_decompressed:
        decompressed_sha256 = hash_decompressed(source_root / "genome_plus_mt.fa.gz")
        if decompressed_sha256 != EXPECTED_MAPPING_SHA256:
            raise ValueError("Decompressed combined FASTA SHA256 mismatch")
    return {
        "status": "PASS",
        "source_root": str(source_root),
        "files": files,
        "fai": {**fai, "total_bp": EXPECTED_TOTAL_BP},
        "decompressed_sha256": decompressed_sha256,
    }


def existing_ancestor(path: Path) -> Path:
    candidate = path
    while not candidate.exists():
        if candidate.parent == candidate:
            raise ValueError(f"No existing ancestor for {path}")
        candidate = candidate.parent
    return candidate


def validate_destination(source_root: Path, destination: Path) -> dict[str, object]:
    if not destination.is_absolute():
        raise ValueError("Canonical destination must be an explicit absolute path")
    if destination.exists() or destination.is_symlink():
        try:
            SUPPORT.reference_inventory(destination)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            raise ValueError(f"Partial or invalid destination exists: {destination}") from error
        raise ValueError(f"Complete destination already exists: {destination}")
    if destination.parent.is_dir():
        stale = sorted(destination.parent.glob(".chipseq.provisioning.*.staging"))
        if stale:
            raise ValueError(f"Partial provisioning staging directory exists: {stale[0]}")
    ancestor = existing_ancestor(destination)
    source_devices = {(source_root / name).stat().st_dev for name in SUPPORT.REFERENCE_SOURCE_SHA256}
    if len(source_devices) != 1 or ancestor.stat().st_dev not in source_devices:
        raise ValueError("Historical sources and canonical destination are not on one filesystem")
    return {
        "destination": str(destination),
        "existing_ancestor": str(ancestor),
        "device": ancestor.stat().st_dev,
        "collision": False,
        "source_reuse_mode": "HARDLINK",
        "hardlink_supported_assessment": "same POSIX CephFS; actual ln is fail-closed",
    }


def current_used_bytes(quota_root: Path) -> int:
    try:
        raw = os.getxattr(quota_root, "ceph.dir.rbytes")
        value = int(raw.decode("ascii"))
    except (AttributeError, OSError, UnicodeError, ValueError) as error:
        raise ValueError(f"Cannot read Ceph recursive-byte usage for {quota_root}") from error
    if value < 0:
        raise ValueError("Negative quota usage is invalid")
    return value


def storage_report(
    quota_bytes: int,
    used_bytes: int,
    generated_index_bytes: int,
    inventory_bytes: int,
    reserve_bytes: int = P10_RESERVE_BYTES,
) -> dict[str, int | str | bool]:
    values = (quota_bytes, used_bytes, generated_index_bytes, inventory_bytes, reserve_bytes)
    if any(not isinstance(value, int) or value < 0 for value in values):
        raise ValueError("Storage quantities must be nonnegative integer bytes")
    if used_bytes > quota_bytes:
        raise ValueError("Current usage exceeds configured quota")
    additional = generated_index_bytes + inventory_bytes
    current_free = quota_bytes - used_bytes
    expected_free = current_free - additional
    after_reserve = expected_free - reserve_bytes
    report: dict[str, int | str | bool] = {
        "quota_source": "configured_30_GiB_home_quota",
        "current_quota_bytes": quota_bytes,
        "current_used_bytes": used_bytes,
        "current_free_bytes": current_free,
        "generated_index_bytes": generated_index_bytes,
        "hardlinked_source_additional_bytes": 0,
        "inventory_bytes": inventory_bytes,
        "additional_persistent_bytes_required": additional,
        "expected_free_bytes_after_publication": expected_free,
        "p10_reserve_bytes": reserve_bytes,
        "expected_free_bytes_after_publication_and_p10_reserve": after_reserve,
        "fits_reference_and_p10_reserve": after_reserve >= 0,
        "additional_operational_cushion_bytes": "not_defined_review_manually",
    }
    if after_reserve < 0:
        raise ValueError("Insufficient quota for the reference plus defined P10 reserve")
    return report


def verify_current_reserve(
    quota_root: Path = QUOTA_ROOT,
    quota_bytes: int = QUOTA_BYTES,
    reserve_bytes: int = P10_RESERVE_BYTES,
) -> dict[str, int | bool]:
    used = current_used_bytes(quota_root)
    free = quota_bytes - used
    if free < reserve_bytes:
        raise ValueError("Current free quota no longer preserves the defined P10 reserve")
    return {
        "current_quota_bytes": quota_bytes,
        "current_used_bytes": used,
        "current_free_bytes": free,
        "p10_reserve_bytes": reserve_bytes,
        "p10_reserve_preserved": True,
    }


def component_report(root: Path, quota_root: Path = QUOTA_ROOT) -> dict[str, object]:
    inventory = SUPPORT.reference_inventory(root)
    index_paths = [root / relative for relative in SUPPORT.INDEX_COMPONENTS]
    generated_index_bytes = sum(path.stat().st_size for path in index_paths)
    report = storage_report(
        QUOTA_BYTES,
        current_used_bytes(quota_root),
        generated_index_bytes,
        (root / SUPPORT.REFERENCE_INVENTORY).stat().st_size,
    )
    return {
        "status": "VALIDATED_BEFORE_PUBLICATION",
        "canonical_destination": str(CANONICAL_ROOT),
        "source_reuse_mode": "HARDLINK",
        "components": inventory["files"],
        "inventory_file": inventory["inventory_file"],
        "storage": report,
    }


def git(*args: str) -> str:
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()


def repository_preflight(allow_dirty: bool = False) -> dict[str, object]:
    branch = git("branch", "--show-current")
    if branch != EXPECTED_BRANCH:
        raise ValueError(f"Expected branch {EXPECTED_BRANCH}; observed {branch}")
    head = git("rev-parse", "HEAD")
    subprocess.run(
        ["git", "-C", str(ROOT), "merge-base", "--is-ancestor", BASELINE_COMMIT, head],
        check=True,
    )
    status = git("status", "--short")
    if status and not allow_dirty:
        raise ValueError("Provisioning source worktree is dirty")
    return {
        "branch": branch,
        "head": head,
        "baseline_commit": BASELINE_COMMIT,
        "baseline_is_ancestor": True,
        "worktree_clean": not bool(status),
        "status": status.splitlines(),
    }


def preflight(allow_dirty: bool = False) -> dict[str, object]:
    repository = repository_preflight(allow_dirty)
    source = validate_source(SOURCE_ROOT)
    destination = validate_destination(SOURCE_ROOT, CANONICAL_ROOT)
    used = current_used_bytes(QUOTA_ROOT)
    return {
        "mode": "check",
        "repository": repository,
        "source": source,
        "destination": destination,
        "quota": {
            "quota_source": "configured_30_GiB_home_quota",
            "current_quota_bytes": QUOTA_BYTES,
            "current_used_bytes": used,
            "current_free_bytes": QUOTA_BYTES - used,
            "generated_index_bytes": "unknown_until_TMPDIR_build",
            "additional_persistent_bytes_required": "unknown_until_TMPDIR_build",
            "p10_reserve_bytes": P10_RESERVE_BYTES,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    source = sub.add_parser("validate-source")
    source.add_argument("--source-root", type=Path, default=SOURCE_ROOT)
    source.add_argument("--output", type=Path)
    destination = sub.add_parser("validate-destination")
    destination.add_argument("--source-root", type=Path, default=SOURCE_ROOT)
    destination.add_argument("--destination", type=Path, default=CANONICAL_ROOT)
    components = sub.add_parser("component-report")
    components.add_argument("--root", type=Path, required=True)
    components.add_argument("--output", type=Path)
    sub.add_parser("verify-current-reserve")
    args = parser.parse_args()

    if args.action == "validate-source":
        result = validate_source(args.source_root)
    elif args.action == "validate-destination":
        result = validate_destination(args.source_root, args.destination)
    elif args.action == "component-report":
        result = component_report(args.root)
    else:
        result = verify_current_reserve()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if getattr(args, "output", None):
        SUPPORT.atomic_text(args.output, text)
    else:
        print(text, end="")


if __name__ == "__main__":
    try:
        main()
    except (
        OSError,
        ValueError,
        KeyError,
        json.JSONDecodeError,
        subprocess.CalledProcessError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2) from error
