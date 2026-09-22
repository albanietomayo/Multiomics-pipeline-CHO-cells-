#!/usr/bin/env python3
"""Validate or explicitly submit one serialized lean ChIP-seq analysis."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
import re
import secrets
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXPECTED_BRANCH = "chipseq-lean-production"
SUPPORT_PATH = ROOT / "workflow/scripts/chipseq_production.py"
spec = importlib.util.spec_from_file_location("chipseq_production_submit_support", SUPPORT_PATH)
SUPPORT = importlib.util.module_from_spec(spec); assert spec.loader is not None; spec.loader.exec_module(SUPPORT)
SOURCE_FILES = [
    "environment.yml", "AGENTS.md", "config/samples.tsv", "config/chipseq_alignment.json",
    "config/chipseq_filtering.json", "config/chipseq_peak_calling_policy.json",
    "config/chipseq_peak_calling_plan.tsv", "config/chipseq_peak_calling_plan_summary.json",
    "snapshots/chipseq/incremental_planning_validation_001/chipseq_analysis_plan.tsv",
    "snapshots/chipseq/runtime_planning_validation_001/chipseq_processing_runs.tsv",
    "workflow/envs/chipseq_production.yaml", "workflow/scripts/download_fastq.py",
    "workflow/scripts/chipseq_alignment_support.py", "workflow/scripts/chipseq_filtering_support.py",
    "workflow/scripts/chipseq_peak_calling_dynamic.py", "workflow/scripts/chipseq_benchmark.py",
    "workflow/scripts/chipseq_production.py", "workflow/slurm/chipseq_production.sbatch",
    "workflow/slurm/submit_chipseq_production.py", "workflow/tests/test_chipseq_production.py",
    "docs/chipseq_lean_production.md",
]

def git(*args): return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()
def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""): digest.update(block)
    return digest.hexdigest()

def validate_sources():
    for name in SOURCE_FILES:
        if not (ROOT / name).is_file(): raise ValueError(f"Missing production source: {name}")
    subprocess.run(["bash", "-n", str(ROOT / "workflow/slurm/chipseq_production.sbatch")], check=True)
    subprocess.run(["python3", "-m", "py_compile", str(SUPPORT_PATH),
                    str(ROOT / "workflow/slurm/submit_chipseq_production.py")], check=True)

def local_plan(args):
    if git("branch", "--show-current") != EXPECTED_BRANCH:
        raise ValueError(f"Expected branch {EXPECTED_BRANCH}")
    status = git("status", "--short")
    if status and (args.submit or not args.development_dirty_check):
        raise ValueError("Production source must be clean (dirty override is check-only)")
    validate_sources()
    value = SUPPORT.preflight(args.analysis_id, args.control_run, args.shared_reference_root,
                              args.scratch_root, args.production_root)
    value["source_git_commit"] = git("rev-parse", "HEAD")
    value["source_worktree_clean"] = not bool(status)
    value["source_files"] = [{"path": name, "bytes": (ROOT/name).stat().st_size,
                              "sha256": sha(ROOT/name)} for name in SOURCE_FILES]
    return value

def prepare(plan, root):
    analysis = plan["analysis"]["analysis_id"]
    final = (root / "analyses" / analysis).resolve()
    final.parent.mkdir(parents=True, exist_ok=True)
    if final.exists(): raise ValueError(f"Production analysis output collision: {final}")
    temporary = Path(tempfile.mkdtemp(prefix=f".{analysis}_", dir=final.parent))
    try:
        project = temporary / "project"
        for item in plan["source_files"]:
            target = project / item["path"]; target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / item["path"], target)
        SUPPORT.atomic_text(temporary / "source.sha256", "".join(
            f"{item['sha256']}  {item['path']}\n" for item in plan["source_files"]))
        SUPPORT.atomic_json(temporary / "submission_plan.json", plan)
        temporary.replace(final)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True); raise
    return final

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-id", default="SRR20770297")
    parser.add_argument("--control-run", default="SRR20770287")
    parser.add_argument("--shared-reference-root", required=True, type=Path)
    parser.add_argument("--production-root", type=Path, default=ROOT / "results/chipseq/production")
    parser.add_argument("--scratch-root", type=Path, default=Path(os.environ.get("TMPDIR", "/tmp")))
    parser.add_argument("--development-dirty-check", action="store_true")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    if args.submit and args.development_dirty_check: raise ValueError("Dirty submission is forbidden")
    plan = local_plan(args)
    if not args.submit:
        print(json.dumps(plan, indent=2, sort_keys=True)); print("CHECK_ONLY_NO_SBATCH=PASS"); return
    token = secrets.token_hex(16)
    slot = SUPPORT.reserve_production_slot(args.production_root, args.analysis_id, token)
    output = None
    try:
        output = prepare(plan, args.production_root)
    except BaseException:
        SUPPORT.release_production_slot(args.production_root, slot, args.analysis_id, token)
        raise
    command = ["sbatch", "--parsable", "--export=ALL," +
        f"CHIP_ANALYSIS_ID={args.analysis_id},CHIP_CONTROL_RUN={args.control_run}," +
        f"CHIP_OUTPUT_DIR={output},CHIP_PRODUCTION_ROOT={args.production_root.resolve()}," +
        f"CHIP_PRODUCTION_SLOT={slot},CHIP_PRODUCTION_SLOT_TOKEN={token}",
        "workflow/slurm/chipseq_production.sbatch"]
    try:
        result = subprocess.check_output(command, cwd=output / "project", text=True).strip()
    except BaseException:
        shutil.rmtree(output)
        SUPPORT.release_production_slot(args.production_root, slot, args.analysis_id, token)
        raise
    if not re.fullmatch(r"[0-9]+(?:;[^\n]+)?", result): raise ValueError(f"Unexpected sbatch response: {result!r}")
    SUPPORT.atomic_text(output / "job_id.txt", result + "\n")
    print(f"SUBMITTED={result}"); print(f"PRODUCTION_SLOT={slot}"); print(f"OUTPUT_DIRECTORY={output}")

if __name__ == "__main__":
    try: main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr); raise SystemExit(2) from error
