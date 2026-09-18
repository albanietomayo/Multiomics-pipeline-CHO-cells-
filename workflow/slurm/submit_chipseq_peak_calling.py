#!/usr/bin/env python3
"""Validate or explicitly submit the dynamic Phase 6B peak-calling workflow."""

import argparse
import hashlib
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXPECTED_BRANCH = "chipseq-incremental-generalization"
FILES = [
    "environment.yml",
    "AGENTS.md",
    "config/chipseq_peak_calling_dynamic.json",
    "config/chipseq_peak_calling_policy.json",
    "config/chipseq_peak_calling_plan.tsv",
    "config/chipseq_peak_calling_plan_summary.json",
    "snapshots/chipseq/incremental_planning_validation_001/chipseq_analysis_plan.tsv",
    "workflow/rules/chipseq_peak_calling.smk",
    "workflow/scripts/build_chipseq_peak_calling_plan.py",
    "workflow/scripts/chipseq_peak_calling_dynamic.py",
    "workflow/scripts/prepare_chipseq_peak_calling_runtime.py",
    "workflow/envs/chipseq_peak_calling.yaml",
    "workflow/envs/chipseq_phantompeakqualtools.yaml",
    "workflow/slurm/chipseq_peak_calling.sbatch",
    "workflow/slurm/submit_chipseq_peak_calling.py",
    "workflow/tests/test_chipseq_peak_calling_dynamic.py",
    "workflow/tests/test_chipseq_peak_calling_dynamic_runtime.py",
    "workflow/tests/test_chipseq_peak_calling_plan.py",
    "workflow/tests/test_chipseq_phase6c_operational.py",
]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def local_validation():
    branch = subprocess.check_output(
        ["git", "branch", "--show-current"], cwd=ROOT, text=True
    ).strip()
    if branch != EXPECTED_BRANCH:
        raise ValueError(f"Expected branch {EXPECTED_BRANCH}; observed {branch}")
    for name in FILES:
        if not (ROOT / name).is_file():
            raise ValueError(f"Missing dynamic peak-calling source: {name}")
    subprocess.run(
        ["bash", "-n", str(ROOT / "workflow/slurm/chipseq_peak_calling.sbatch")],
        check=True,
    )
    subprocess.run(
        ["python3", "-m", "unittest", "workflow.tests.test_chipseq_peak_calling_plan",
         "workflow.tests.test_chipseq_peak_calling_dynamic",
         "workflow.tests.test_chipseq_peak_calling_dynamic_runtime",
         "workflow.tests.test_chipseq_phase6c_operational"],
        cwd=ROOT,
        check=True,
    )
    subprocess.run(
        ["python3", "workflow/scripts/prepare_chipseq_peak_calling_runtime.py", "--check"],
        cwd=ROOT,
        check=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--check", action="store_true")
    modes.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    local_validation()
    if args.check:
        print("[OK] Dynamic peak-calling check passed; no job submitted.")
        return

    parent = ROOT / "results/chipseq/peak_calling/slurm"
    parent.mkdir(parents=True, exist_ok=True)
    submission = Path(tempfile.mkdtemp(prefix="submission_", dir=parent))
    project = submission / "project"
    for name in FILES:
        target = project / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
    subprocess.run(
        ["python3", str(ROOT / "workflow/scripts/prepare_chipseq_peak_calling_runtime.py"),
         "--output", str(project / "config/chipseq_peak_calling_runtime.json")],
        cwd=ROOT,
        check=True,
    )
    subprocess.run(
        ["snakemake", "--snakefile", "workflow/rules/chipseq_peak_calling.smk",
         "chipseq_peak_calling_all", "--cores", "1", "--dry-run"],
        cwd=project,
        check=True,
    )
    manifest = submission / "source.sha256"
    manifest.write_text("".join(
        f"{sha256(path)}  {path.relative_to(project).as_posix()}\n"
        for path in sorted(project.rglob("*")) if path.is_file()
    ), encoding="utf-8")
    (submission / "base_commit.txt").write_bytes(
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT)
    )
    (submission / "git_status.txt").write_bytes(
        subprocess.check_output(["git", "status", "--short"], cwd=ROOT)
    )
    result = subprocess.check_output([
        "sbatch", "--parsable", f"--chdir={project}",
        f"--output={submission}/slurm_%j.out",
        str(project / "workflow/slurm/chipseq_peak_calling.sbatch"),
    ], text=True).strip()
    job = result.split(";", 1)[0]
    if not re.fullmatch(r"[0-9]+", job):
        raise ValueError(f"Unexpected sbatch response: {result}")
    (submission / "job_id.txt").write_text(job + "\n", encoding="utf-8")
    (parent / "latest_submission.txt").write_text(str(submission) + "\n", encoding="utf-8")
    print(f"JOB SUBMITTED: {job}")


if __name__ == "__main__":
    main()
