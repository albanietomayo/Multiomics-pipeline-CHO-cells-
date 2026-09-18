#!/usr/bin/env python3
"""Prepare the dynamic peak runtime only from a completed filtering cohort."""

import argparse
import importlib.util
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/chipseq_peak_calling_dynamic.json"
SUPPORT = ROOT / "workflow/scripts/chipseq_peak_calling_dynamic.py"


def load_support():
    spec = importlib.util.spec_from_file_location("dynamic_peak", SUPPORT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def completed_filtering_job(pointer):
    pointer = Path(pointer)
    if not pointer.is_file():
        raise ValueError(
            "No real productive dynamic filtering pointer exists; peak calling remains fail-closed"
        )
    submission = Path(pointer.read_text(encoding="utf-8").strip())
    if not submission.is_absolute() or not re.fullmatch(r"submission_[A-Za-z0-9_]+", submission.name):
        raise ValueError("Filtering latest-submission pointer is unsafe")
    if submission.parent.resolve() != pointer.parent.resolve() or not submission.is_dir():
        raise ValueError("Filtering latest-submission pointer escapes its expected directory")
    job_id_path = submission / "job_id.txt"
    if not job_id_path.is_file():
        raise ValueError("Filtering submission has no productive job_id.txt")
    job_id = job_id_path.read_text(encoding="utf-8").strip()
    if not re.fullmatch(r"[0-9]+", job_id):
        raise ValueError("Filtering submission has an invalid job ID")
    job = submission / f"job_{job_id}"
    if not job.is_dir():
        raise ValueError("Completed filtering job directory is missing")
    return job


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="validate without writing runtime JSON")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.check and args.output:
        raise ValueError("--check cannot be combined with --output")
    if not args.check and args.output is None:
        raise ValueError("Use --check or provide --output")
    cfg = json.loads(CONFIG.read_text(encoding="utf-8"))
    if cfg.get("schema_version") != 1 or cfg.get("scope") != "dynamic_chipseq_peak_calling_runtime_preflight":
        raise ValueError("Unsupported dynamic peak-calling preflight configuration")
    support = load_support()
    job = completed_filtering_job(ROOT / cfg["filtering_latest_submission"])
    runtime = support.build_runtime_manifest(
        ROOT / cfg["peak_plan"],
        ROOT / cfg["peak_plan_summary"],
        ROOT / cfg["peak_policy"],
        ROOT / cfg["analysis_plan"],
        job,
        verify_large_files=True,
    )
    if runtime["analysis_count"] != cfg["expected_analysis_count"]:
        raise ValueError("Unexpected runtime analysis count")
    if runtime["run_count"] != cfg["expected_run_count"]:
        raise ValueError("Unexpected runtime referenced-run count")
    if args.output:
        output = args.output if args.output.is_absolute() else ROOT / args.output
        support.validate_output_destination(output)
        support.atomic_json(output, runtime)
        print(f"[OK] Wrote verified dynamic peak runtime: {output}")
    else:
        print("[OK] Dynamic peak runtime preflight passed; no files written")


if __name__ == "__main__":
    main()
