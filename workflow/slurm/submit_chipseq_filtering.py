#!/usr/bin/env python3
"""Prepare an immutable SINGLE ChIP filtering submission; optionally submit it."""
import argparse
import importlib.util
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FILES = [
    "environment.yml", "config/samples.tsv", "config/chipseq_pilot.json",
    "config/chipseq_experimental_eligibility.json", "config/chipseq_filtering.json",
    "workflow/rules/chipseq_eligibility.smk", "workflow/rules/chipseq_filtering.smk",
    "workflow/scripts/validate_chipseq_eligibility.py", "workflow/scripts/chipseq_filtering_support.py",
    "workflow/envs/chipseq_duplicates.yaml", "workflow/envs/chipseq_filtering.yaml",
    "workflow/slurm/chipseq_filtering.sbatch", "workflow/slurm/submit_chipseq_filtering.py",
    "workflow/tests/test_chipseq_filtering.py", "docs/chipseq_filtering_pilot.md",
    "snapshots/chipseq/control_validation_001/chipseq_conditions.tsv",
    "snapshots/chipseq/control_validation_001/control_candidate_summary.json",
    "snapshots/chipseq/control_validation_001/request_plan_provenance.json",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--check", action="store_true")
    modes.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    if subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip() != "chipseq-metadata-audit":
        raise ValueError("Expected branch chipseq-metadata-audit")
    spec = importlib.util.spec_from_file_location("chip_filter_support", ROOT / "workflow/scripts/chipseq_filtering_support.py")
    support = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(support)
    plan = support.prepare(ROOT)
    cfg = support.configuration(ROOT)
    evidence = support.relative(cfg["alignment_evidence"])
    files = FILES + [(evidence / name).as_posix() for name in ("SHA256SUMS.txt", "job/output.sha256", "job/job_status.tsv")]
    for name in files:
        if not (ROOT / name).is_file():
            raise ValueError(f"Missing source: {name}")
    subprocess.run(["bash", "-n", str(ROOT / "workflow/slurm/chipseq_filtering.sbatch")], check=True)
    parent = ROOT / "results/chipseq/filtering/slurm"
    parent.mkdir(parents=True, exist_ok=True)
    submission = Path(tempfile.mkdtemp(prefix="submission_", dir=parent))
    project = submission / "project"
    for name in files:
        target = project / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
    support.dump(project / "config/chipseq_filtering_inputs.json", plan)
    for name, digest in plan["small_files"].items():
        relative = Path(name).relative_to("outputs")
        target = project / "inputs/upstream" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(Path(plan["source_root"]) / name, target)
        support.check(target, digest)
        if relative.parts[0] == "reference":
            reference = project / "inputs" / relative
            reference.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(target, reference)
    for command, name in ((["git", "rev-parse", "HEAD"], "base_commit.txt"),
                          (["git", "status", "--short"], "git_status.txt"),
                          (["git", "diff", "--binary", "--", *files], "source_changes.patch")):
        (submission / name).write_bytes(subprocess.check_output(command, cwd=ROOT))
    hashes = "".join(f"{support.sha(p)}  {p.relative_to(project).as_posix()}\n"
                     for p in sorted(project.rglob("*")) if p.is_file())
    (submission / "source.sha256").write_text(hashes)
    subprocess.run(["python3", "workflow/scripts/validate_chipseq_eligibility.py",
                    "--decisions", "config/chipseq_experimental_eligibility.json",
                    "--samples", "config/samples.tsv", "--conditions",
                    "snapshots/chipseq/control_validation_001/chipseq_conditions.tsv",
                    "--pilot", "config/chipseq_pilot.json", "--report", str(submission / "eligibility_preflight.json")],
                   cwd=project, check=True)
    print("[OK] Eligibility, archived alignment evidence, reports and BAM availability checked.", flush=True)
    print("Full BAM hashes will be verified on the compute node before marking duplicates.", flush=True)
    print(f"LOGS: {submission}", flush=True)
    if args.check:
        print("[OK] Check only; no job submitted. VERA workflow validation remains pending.")
        return
    result = subprocess.check_output(["sbatch", "--parsable", f"--chdir={project}",
        f"--output={submission}/slurm_%j.out", str(project / "workflow/slurm/chipseq_filtering.sbatch")], text=True).strip()
    job = result.split(";")[0]
    if not re.fullmatch(r"[0-9]+", job):
        raise ValueError(f"Unexpected sbatch response: {result}")
    (submission / "job_id.txt").write_text(job + "\n")
    (parent / "latest_submission.txt").write_text(str(submission) + "\n")
    print(f"JOB SUBMITTED: {job}", flush=True)
    subprocess.run(["squeue", "-j", job, "-o", "%.18i %.26j %.12T %.10M %R"], check=False)


if __name__ == "__main__":
    main()
