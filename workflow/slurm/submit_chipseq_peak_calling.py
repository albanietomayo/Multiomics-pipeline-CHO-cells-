#!/usr/bin/env python3
"""Prepare an immutable ChIP-seq peak-calling submission; optionally submit it."""

import argparse
import importlib.util
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

FILES = [
    "environment.yml",
    "config/samples.tsv",
    "config/chipseq_pilot.json",
    "config/chipseq_experimental_eligibility.json",
    "config/chipseq_peak_calling.json",
    "workflow/rules/chipseq_eligibility.smk",
    "workflow/rules/chipseq_peak_calling.smk",
    "workflow/scripts/validate_chipseq_eligibility.py",
    "workflow/scripts/chipseq_peak_calling_support.py",
    "workflow/envs/chipseq_peak_calling.yaml",
    "workflow/slurm/chipseq_peak_calling.sbatch",
    "workflow/slurm/submit_chipseq_peak_calling.py",
    "workflow/tests/test_chipseq_peak_calling.py",
    "docs/chipseq_peak_calling_pilot.md",
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

    branch = subprocess.check_output(
        ["git", "branch", "--show-current"],
        cwd=ROOT,
        text=True,
    ).strip()

    if branch != "chipseq-metadata-audit":
        raise ValueError(f"Unexpected branch: {branch}")

    spec = importlib.util.spec_from_file_location(
        "chip_peak_support",
        ROOT / "workflow/scripts/chipseq_peak_calling_support.py",
    )

    support = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(support)

    cfg = support.configuration(ROOT)
    evidence = ROOT / cfg["filtering_evidence"]

    evidence_files = [
        path.relative_to(ROOT).as_posix()
        for path in sorted(evidence.rglob("*"))
        if path.is_file()
    ]

    files = FILES + evidence_files

    if len(files) != len(set(files)):
        raise ValueError("Duplicate source paths in peak-calling submission")

    for name in files:
        if not (ROOT / name).is_file():
            raise ValueError(f"Missing source: {name}")

    subprocess.run(
        ["bash", "-n", str(ROOT / "workflow/slurm/chipseq_peak_calling.sbatch")],
        check=True,
    )

    subprocess.run(
        [
            "python3",
            "-m",
            "py_compile",
            "workflow/scripts/chipseq_peak_calling_support.py",
            "workflow/slurm/submit_chipseq_peak_calling.py",
            "workflow/tests/test_chipseq_peak_calling.py",
        ],
        cwd=ROOT,
        check=True,
    )

    subprocess.run(
        ["python3", "workflow/tests/test_chipseq_peak_calling.py"],
        cwd=ROOT,
        check=True,
    )

    plan = support.prepare(ROOT)

    parent = ROOT / "results/chipseq/peak_calling/slurm"
    parent.mkdir(parents=True, exist_ok=True)

    submission = Path(
        tempfile.mkdtemp(
            prefix="submission_",
            dir=parent,
        )
    )

    project = submission / "project"

    for name in files:
        target = project / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)

    support.dump(
        project / "config/chipseq_peak_calling_inputs.json",
        plan,
    )

    reference_source = Path(plan["reference"]["source_path"])
    reference_target = project / "inputs/reference/genome_plus_mt.fa.fai"

    reference_target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(reference_source, reference_target)

    support.check(
        reference_target,
        plan["reference"]["sha256"],
    )

    commands = [
        (["git", "rev-parse", "HEAD"], "base_commit.txt"),
        (["git", "status", "--short"], "git_status.txt"),
        (
            ["git", "diff", "--binary", "--", *files],
            "source_changes.patch",
        ),
    ]

    for command, name in commands:
        (submission / name).write_bytes(
            subprocess.check_output(command, cwd=ROOT)
        )

    hashes = "".join(
        f"{support.sha(path)}  {path.relative_to(project).as_posix()}\n"
        for path in sorted(project.rglob("*"))
        if path.is_file()
    )

    (submission / "source.sha256").write_text(hashes)

    subprocess.run(
        [
            "python3",
            "workflow/scripts/validate_chipseq_eligibility.py",
            "--decisions",
            "config/chipseq_experimental_eligibility.json",
            "--samples",
            "config/samples.tsv",
            "--conditions",
            "snapshots/chipseq/control_validation_001/chipseq_conditions.tsv",
            "--pilot",
            "config/chipseq_pilot.json",
            "--report",
            str(submission / "eligibility_preflight.json"),
        ],
        cwd=project,
        check=True,
    )

    print(
        "[OK] Tests, eligibility, filtering evidence, "
        "reference span and BAM availability checked.",
        flush=True,
    )

    print(
        "Complete filtered BAM and CSI hashes will be verified "
        "on the compute node before MACS3.",
        flush=True,
    )

    print(f"LOGS: {submission}", flush=True)

    if args.check:
        print(
            "[OK] Check only; no job submitted. "
            "VERA MACS3 validation remains pending."
        )
        return

    result = subprocess.check_output(
        [
            "sbatch",
            "--parsable",
            f"--chdir={project}",
            f"--output={submission}/slurm_%j.out",
            str(project / "workflow/slurm/chipseq_peak_calling.sbatch"),
        ],
        text=True,
    ).strip()

    job = result.split(";")[0]

    if not re.fullmatch(r"[0-9]+", job):
        raise ValueError(f"Unexpected sbatch response: {result}")

    (submission / "job_id.txt").write_text(job + "\n")
    (parent / "latest_submission.txt").write_text(str(submission) + "\n")

    print(f"JOB SUBMITTED: {job}", flush=True)

    subprocess.run(
        [
            "squeue",
            "-j",
            job,
            "-o",
            "%.18i %.26j %.12T %.10M %R",
        ],
        check=False,
    )


if __name__ == "__main__":
    main()
