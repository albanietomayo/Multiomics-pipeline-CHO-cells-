#!/usr/bin/env python3
"""Validate or explicitly submit one ChIP-seq percentile benchmark run."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SUPPORT_PATH = ROOT / "workflow/scripts/chipseq_benchmark.py"
SPEC = importlib.util.spec_from_file_location("chipseq_benchmark_support", SUPPORT_PATH)
SUPPORT = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(SUPPORT)

EXPECTED_BRANCH = "chipseq-incremental-generalization"
SOURCE_FILES = [
    "environment.yml",
    "config/config.yaml",
    "config/samples.tsv",
    "config/chipseq_alignment.json",
    "config/chipseq_filtering.json",
    "workflow/envs/reference.yaml",
    "workflow/envs/preprocessing.yaml",
    "workflow/envs/qc.yaml",
    "workflow/envs/chipseq_alignment.yaml",
    "workflow/envs/chipseq_duplicates.yaml",
    "workflow/envs/chipseq_filtering.yaml",
    "workflow/scripts/download_fastq.py",
    "workflow/scripts/fetch_reference_genome.py",
    "workflow/scripts/chipseq_alignment_support.py",
    "workflow/scripts/chipseq_filtering_support.py",
    "workflow/scripts/chipseq_benchmark.py",
    "workflow/slurm/chipseq_benchmark_p10_p50_p90.sbatch",
    "benchmarks/chipseq/2026-09-18/benchmark_selection.tsv",
    "benchmarks/chipseq/2026-09-18/chipseq_benchmark_manifest.tsv",
    "benchmarks/chipseq/2026-09-18/benchmark_execution_plan.tsv",
    "benchmarks/chipseq/2026-09-18/benchmark_metrics.tsv",
    "benchmarks/chipseq/2026-09-18/benchmark_summary.tsv",
    "benchmarks/chipseq/2026-09-18/README.md",
    "benchmarks/chipseq/2026-09-18/SHA256SUMS",
]


def git(*args: str) -> str:
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()


def sha256(path: Path) -> str:
    return SUPPORT.sha256(path)


def source_inventory() -> list[dict[str, object]]:
    inventory = []
    for relative in SOURCE_FILES:
        path = ROOT / relative
        if not path.is_file():
            raise ValueError(f"Required benchmark source is missing: {relative}")
        inventory.append({"path": relative, "bytes": path.stat().st_size, "sha256": sha256(path)})
    return inventory


def local_validation(args: argparse.Namespace) -> dict[str, object]:
    branch = git("branch", "--show-current")
    if branch != EXPECTED_BRANCH:
        raise ValueError(f"Expected branch {EXPECTED_BRANCH}; observed {branch}")
    status = git("status", "--short")
    if status and not (args.development_dirty_check and not args.submit):
        raise ValueError("Benchmark source worktree is dirty; commit reviewed changes before submission")
    rows = SUPPORT.validate_manifest()
    row = SUPPORT.manifest_item(rows, args.benchmark_class)
    if args.run_accession is not None and args.run_accession != row["run_accession"]:
        raise ValueError("Explicit run accession does not match the selected percentile")
    storage = SUPPORT.storage_preflight(row, args.output_root, args.scratch_root)
    return {
        "mode": "submit" if args.submit else "check",
        "development_dirty_check": bool(args.development_dirty_check),
        "branch": branch,
        "source_git_commit": git("rev-parse", "HEAD"),
        "source_worktree_clean": not bool(status),
        "source_status": status.splitlines(),
        "manifest_sha256": sha256(SUPPORT.MANIFEST),
        "benchmark": row,
        "storage": storage,
        "source_files": source_inventory(),
        "slurm_requests": {
            "status": "provisional_bootstrap_inherited_from_validated_chipseq_workers",
            "cpus": 8,
            "memory": "32G",
            "walltime": "08:00:00",
        },
    }


def write_json(path: Path, value: object) -> None:
    SUPPORT.atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def prepare_submission(plan: dict[str, object], output_root: Path) -> Path:
    row = plan["benchmark"]
    final = output_root / f"{row['benchmark_class']}_{row['run_accession']}"
    output_root.mkdir(parents=True, exist_ok=True)
    if final.exists():
        raise ValueError(f"Benchmark output collision: {final}")
    temporary = Path(tempfile.mkdtemp(prefix=".chipseq_benchmark_", dir=output_root))
    try:
        project = temporary / "project"
        for relative in SOURCE_FILES:
            destination = project / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / relative, destination)
        write_json(temporary / "submission_plan.json", plan)
        lines = [f"{item['sha256']}  {item['path']}\n" for item in plan["source_files"]]
        SUPPORT.atomic_text(temporary / "source.sha256", "".join(lines))
        temporary.replace(final)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return final


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--benchmark-class", required=True, choices=tuple(SUPPORT.EXPECTED))
    result.add_argument("--run-accession", help="Optional redundant fail-closed accession check")
    result.add_argument("--output-root", type=Path, default=ROOT / "results/chipseq/benchmark/runs")
    result.add_argument("--scratch-root", type=Path, default=Path(os.environ.get("TMPDIR", "/tmp")))
    result.add_argument("--submit", action="store_true", help="Explicitly create and submit; default is check only")
    result.add_argument(
        "--development-dirty-check",
        action="store_true",
        help="Allow dirty source only in non-submitting check mode",
    )
    return result


def main() -> None:
    args = parser().parse_args()
    if args.submit and args.development_dirty_check:
        raise ValueError("Dirty-source override is forbidden for submission")
    plan = local_validation(args)
    if not args.submit:
        print(json.dumps(plan, indent=2, sort_keys=True))
        print("CHECK_ONLY_NO_SBATCH=PASS")
        return
    submission = prepare_submission(plan, args.output_root)
    project = submission / "project"
    row = plan["benchmark"]
    command = [
        "sbatch",
        "--parsable",
        "--export=ALL,"
        f"CHIP_BENCHMARK_CLASS={row['benchmark_class']},"
        f"CHIP_BENCHMARK_ACCESSION={row['run_accession']}",
        "workflow/slurm/chipseq_benchmark_p10_p50_p90.sbatch",
    ]
    job_id = subprocess.check_output(command, cwd=project, text=True).strip()
    if not job_id or not job_id.split(";", 1)[0].isdigit():
        raise ValueError(f"Unexpected sbatch response: {job_id!r}")
    SUPPORT.atomic_text(submission / "job_id.txt", job_id + "\n")
    print(f"SUBMITTED={job_id}")
    print(f"SUBMISSION_DIRECTORY={submission}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2) from error
