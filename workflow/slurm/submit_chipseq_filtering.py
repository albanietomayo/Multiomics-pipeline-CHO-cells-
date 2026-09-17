#!/usr/bin/env python3

"""Prepare a dynamic ChIP-seq filtering submission.

--check validates and snapshots the dynamic alignment upstream
without calling sbatch.

--submit performs the same validation and then submits filtering.
"""

import argparse
import importlib.util
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(
    __file__
).resolve().parents[2]

EXPECTED_BRANCH = (
    "chipseq-incremental-generalization"
)

FILES = [
    "environment.yml",
    "config/chipseq_filtering.json",
    "workflow/rules/chipseq_filtering.smk",
    "workflow/scripts/chipseq_filtering_support.py",
    "workflow/envs/chipseq_duplicates.yaml",
    "workflow/envs/chipseq_filtering.yaml",
    "workflow/slurm/chipseq_filtering.sbatch",
    "workflow/slurm/submit_chipseq_filtering.py",
    "workflow/tests/test_chipseq_filtering.py",
    "workflow/tests/test_chipseq_filtering_dynamic_runtime.py",
    "workflow/tests/test_chipseq_filtering_operational.py",
    "docs/chipseq_incremental_generalization.md",
]


def load_support():
    script = (
        ROOT
        / "workflow"
        / "scripts"
        / "chipseq_filtering_support.py"
    )

    spec = importlib.util.spec_from_file_location(
        "chipseq_filtering_support",
        script,
    )

    module = importlib.util.module_from_spec(
        spec
    )

    spec.loader.exec_module(
        module
    )

    return module


def main():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    modes = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    modes.add_argument(
        "--check",
        action="store_true",
    )

    modes.add_argument(
        "--submit",
        action="store_true",
    )

    args = parser.parse_args()

    branch = subprocess.check_output(
        [
            "git",
            "branch",
            "--show-current",
        ],
        cwd=ROOT,
        text=True,
    ).strip()

    if branch != EXPECTED_BRANCH:
        raise ValueError(
            f"Expected branch {EXPECTED_BRANCH}; "
            f"observed {branch}"
        )

    support = load_support()

    plan = support.prepare(
        ROOT
    )

    for name in FILES:

        if not (
            ROOT
            / name
        ).is_file():

            raise ValueError(
                f"Missing source: {name}"
            )

    subprocess.run(
        [
            "bash",
            "-n",
            str(
                ROOT
                / "workflow"
                / "slurm"
                / "chipseq_filtering.sbatch"
            ),
        ],
        check=True,
    )

    parent = (
        ROOT
        / "results"
        / "chipseq"
        / "filtering"
        / "slurm"
    )

    parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    submission = Path(
        tempfile.mkdtemp(
            prefix="submission_",
            dir=parent,
        )
    )

    project = (
        submission
        / "project"
    )

    for name in FILES:

        source = (
            ROOT
            / name
        )

        target = (
            project
            / name
        )

        target.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        shutil.copyfile(
            source,
            target,
        )

    support.dump(
        project
        / "config"
        / "chipseq_filtering_inputs.json",
        plan,
    )

    source_root = Path(
        plan[
            "source_root"
        ]
    )

    for name, digest in sorted(
        plan[
            "small_files"
        ].items()
    ):

        relative = Path(
            name
        ).relative_to(
            "outputs"
        )

        source = (
            source_root
            / name
        )

        target = (
            project
            / "inputs"
            / "upstream"
            / relative
        )

        target.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        shutil.copyfile(
            source,
            target,
        )

        support.check(
            target,
            digest,
        )

        if (
            relative.parts
            and relative.parts[0]
            == "reference"
        ):

            reference = (
                project
                / "inputs"
                / relative
            )

            reference.parent.mkdir(
                parents=True,
                exist_ok=True,
            )

            shutil.copyfile(
                target,
                reference,
            )

            support.check(
                reference,
                digest,
            )

    upstream_record = {
        "schema_version": 1,
        "alignment_submission": plan[
            "alignment_submission"
        ],
        "alignment_job_id": plan[
            "alignment_job_id"
        ],
        "alignment_job_directory": plan[
            "source_root"
        ],
        "alignment_source_manifest_sha256": plan[
            "source_manifest_sha256"
        ],
        "alignment_job_status_sha256": plan[
            "job_status_sha256"
        ],
        "alignment_run_count": plan[
            "run_count"
        ],
        "alignment_role_counts": plan[
            "role_counts"
        ],
        "library_layout": plan[
            "library_layout"
        ],
        "instrument_platform": plan[
            "instrument_platform"
        ],
        "full_bam_sha256_verification_deferred_to_compute": True,
        "upstream_raw_bam_cleanup_authorized": False,
    }

    support.dump(
        submission
        / "upstream_alignment.json",
        upstream_record,
    )

    commands = (
        (
            [
                "git",
                "rev-parse",
                "HEAD",
            ],
            "base_commit.txt",
        ),
        (
            [
                "git",
                "status",
                "--short",
            ],
            "git_status.txt",
        ),
        (
            [
                "git",
                "diff",
                "--binary",
                "--",
                *FILES,
            ],
            "source_changes.patch",
        ),
    )

    for command, name in commands:

        (
            submission
            / name
        ).write_bytes(
            subprocess.check_output(
                command,
                cwd=ROOT,
            )
        )

    hashes = "".join(
        (
            f"{support.sha(path)}  "
            f"{path.relative_to(project).as_posix()}\n"
        )
        for path in sorted(
            project.rglob("*")
        )
        if path.is_file()
    )

    (
        submission
        / "source.sha256"
    ).write_text(
        hashes,
        encoding="utf-8",
    )

    print(
        "[OK] Dynamic ChIP filtering source snapshot "
        "and alignment upstream validated.",
        flush=True,
    )

    print(
        "FILTERING RUNS: "
        f"{plan['run_count']} "
        f"(IP={plan['role_counts']['ip']}, "
        f"Input={plan['role_counts']['input']})",
        flush=True,
    )

    print(
        "UPSTREAM ALIGNMENT: "
        + plan[
            "alignment_submission"
        ],
        flush=True,
    )

    if args.check:

        print(
            f"CHECK DIRECTORY: {submission}",
            flush=True,
        )

        print(
            "[OK] Check only. No job submitted.",
            flush=True,
        )

        return

    result = subprocess.check_output(
        [
            "sbatch",
            "--parsable",
            f"--chdir={project}",
            f"--output={submission}/slurm_%j.out",
            str(
                project
                / "workflow"
                / "slurm"
                / "chipseq_filtering.sbatch"
            ),
        ],
        text=True,
    ).strip()

    job = result.split(
        ";"
    )[0]

    if not re.fullmatch(
        r"[0-9]+",
        job,
    ):
        raise ValueError(
            f"Unexpected sbatch response: {result}"
        )

    (
        submission
        / "job_id.txt"
    ).write_text(
        job + "\n",
        encoding="utf-8",
    )

    (
        parent
        / "latest_submission.txt"
    ).write_text(
        str(
            submission
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        f"JOB SUBMITTED: {job}",
        flush=True,
    )

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
