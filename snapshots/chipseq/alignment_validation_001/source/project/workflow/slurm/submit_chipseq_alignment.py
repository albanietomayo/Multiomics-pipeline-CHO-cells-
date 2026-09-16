#!/usr/bin/env python3
"""Snapshot and validate a ChIP pilot alignment submission; optionally submit."""
import argparse
import importlib.util
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FILES = [
    "environment.yml", "config/config.yaml", "config/samples.tsv", "config/chipseq_pilot.json",
    "config/chipseq_experimental_eligibility.json", "config/chipseq_alignment.json",
    "workflow/rules/chipseq_alignment.smk", "workflow/rules/chipseq_eligibility.smk",
    "workflow/scripts/validate_chipseq_eligibility.py", "workflow/scripts/fetch_reference_genome.py",
    "workflow/scripts/chipseq_alignment_support.py", "workflow/envs/reference.yaml",
    "workflow/envs/chipseq_alignment.yaml", "workflow/slurm/chipseq_alignment.sbatch",
    "workflow/slurm/submit_chipseq_alignment.py", "docs/chipseq_alignment_pilot.md",
    "snapshots/chipseq/control_validation_001/chipseq_conditions.tsv",
    "snapshots/chipseq/control_validation_001/control_candidate_summary.json",
    "snapshots/chipseq/control_validation_001/request_plan_provenance.json",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    branch = subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip()
    if branch != "chipseq-metadata-audit":
        raise ValueError(f"Unexpected branch: {branch}")
    for name in FILES:
        if not (ROOT / name).is_file():
            raise ValueError(f"Missing source: {name}")
    spec = importlib.util.spec_from_file_location("chip_support", ROOT / "workflow/scripts/chipseq_alignment_support.py")
    support = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(support)
    plan = support.prepare(ROOT)
    subprocess.run(["bash", "-n", str(ROOT / "workflow/slurm/chipseq_alignment.sbatch")], check=True)
    parent = ROOT / "results/chipseq/alignment/slurm"
    parent.mkdir(parents=True, exist_ok=True)
    submission = Path(tempfile.mkdtemp(prefix="submission_", dir=parent))
    project = submission / "project"
    for name in FILES:
        target = project / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
    support.dump(project / "config/chipseq_alignment_inputs.json", plan)
    for command, name in ((["git", "rev-parse", "HEAD"], "base_commit.txt"),
                          (["git", "status", "--short"], "git_status.txt"),
                          (["git", "diff", "--binary", "--", *FILES], "source_changes.patch")):
        (submission / name).write_bytes(subprocess.check_output(command, cwd=ROOT))
    hashes = []
    for path in sorted(project.rglob("*")):
        if path.is_file():
            hashes.append(f"{support.sha(path)}  {path.relative_to(project).as_posix()}\n")
    (submission / "source.sha256").write_text("".join(hashes))
    subprocess.run(["python3", "workflow/scripts/validate_chipseq_eligibility.py",
                    "--decisions", "config/chipseq_experimental_eligibility.json",
                    "--samples", "config/samples.tsv",
                    "--conditions", "snapshots/chipseq/control_validation_001/chipseq_conditions.tsv",
                    "--pilot", "config/chipseq_pilot.json",
                    "--report", str(submission / "eligibility_preflight.json")], cwd=project, check=True)
    print("[OK] Source snapshot, eligibility, upstream reports and FASTQ availability checked.", flush=True)
    print("The complete FASTQ hashes will be verified on the compute node before alignment.", flush=True)
    print(f"LOGS: {submission}", flush=True)
    if args.check:
        print("[OK] Check only. No job submitted; VERA execution validation remains pending.")
        return
    result = subprocess.check_output(["sbatch", "--parsable", f"--chdir={project}",
        f"--output={submission}/slurm_%j.out", str(project / "workflow/slurm/chipseq_alignment.sbatch")], text=True).strip()
    job = result.split(";")[0]
    if not re.fullmatch(r"[0-9]+", job):
        raise ValueError(f"Unexpected sbatch response: {result}")
    (submission / "job_id.txt").write_text(job + "\n")
    (parent / "latest_submission.txt").write_text(str(submission) + "\n")
    print(f"JOB SUBMITTED: {job}", flush=True)
    subprocess.run(["squeue", "-j", job, "-o", "%.18i %.26j %.12T %.10M %R"], check=False)


if __name__ == "__main__":
    main()
