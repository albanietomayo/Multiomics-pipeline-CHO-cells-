#!/usr/bin/env python3

"""
Create a reproducible dynamic ChIP-seq alignment submission.

The alignment input plan is generated from the latest successfully
completed dynamic ChIP-seq preprocessing submission.

--check validates and snapshots the submission without calling sbatch.
--submit performs the same validation and then submits the job.
"""

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

EXPECTED_BRANCH = (
    "chipseq-incremental-generalization"
)

FILES = [
    "environment.yml",
    "config/config.yaml",
    "config/chipseq_alignment.json",
    "workflow/rules/chipseq_alignment.smk",
    "workflow/scripts/fetch_reference_genome.py",
    "workflow/scripts/chipseq_alignment_support.py",
    "workflow/scripts/prepare_chipseq_alignment_inputs.py",
    "workflow/envs/reference.yaml",
    "workflow/envs/chipseq_alignment.yaml",
    "workflow/slurm/chipseq_alignment.sbatch",
    "workflow/slurm/submit_chipseq_alignment.py",
    "docs/chipseq_incremental_generalization.md",
]


def sha256(path):
    digest = hashlib.sha256()

    with Path(path).open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def relative_path(value):
    path = Path(value)

    if (
        path.is_absolute()
        or not path.parts
        or ".." in path.parts
    ):
        raise ValueError(
            f"Expected safe project-relative path: {value}"
        )

    return path


def read_job_status(path):
    values = {}

    with Path(path).open(
        encoding="utf-8"
    ) as handle:
        for line in handle:
            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) != 2:
                raise ValueError(
                    f"Malformed job status line: {line!r}"
                )

            key, value = fields

            if key in values:
                raise ValueError(
                    f"Duplicate job status key: {key}"
                )

            values[key] = value

    required = {
        "workflow_exit_status",
        "exit_status",
    }

    if set(values) != required:
        raise ValueError(
            "Unexpected preprocessing job-status keys: "
            f"{sorted(values)}"
        )

    try:
        workflow_exit = int(
            values["workflow_exit_status"]
        )

        final_exit = int(
            values["exit_status"]
        )

    except ValueError as exc:
        raise ValueError(
            "Non-integer preprocessing exit status"
        ) from exc

    return {
        "workflow_exit_status": workflow_exit,
        "exit_status": final_exit,
    }


def resolve_preprocessing_upstream(root):
    root = Path(root).resolve()

    config_path = (
        root
        / "config"
        / "chipseq_alignment.json"
    )

    config = json.loads(
        config_path.read_text(
            encoding="utf-8"
        )
    )

    if config.get("schema_version") != 1:
        raise ValueError(
            "Unsupported alignment configuration schema"
        )

    pointer_rel = relative_path(
        config[
            "preprocessing_latest_submission"
        ]
    )

    pointer = root / pointer_rel

    if not pointer.is_file():
        raise ValueError(
            "No submitted dynamic ChIP preprocessing is "
            f"available: {pointer}"
        )

    raw_submission = pointer.read_text(
        encoding="utf-8"
    ).strip()

    if not raw_submission:
        raise ValueError(
            "Dynamic preprocessing latest-submission "
            "pointer is empty"
        )

    submission = Path(
        raw_submission
    )

    if not submission.is_absolute():
        submission = root / submission

    submission = submission.resolve()

    allowed_parent = (
        root
        / "results"
        / "chipseq"
        / "preprocessing"
        / "slurm"
    ).resolve()

    if submission.parent != allowed_parent:
        raise ValueError(
            "Preprocessing submission is outside the "
            "expected project SLURM directory"
        )

    if not re.fullmatch(
        r"submission_[A-Za-z0-9_-]+",
        submission.name,
    ):
        raise ValueError(
            f"Unexpected preprocessing submission name: "
            f"{submission.name}"
        )

    job_id_file = (
        submission
        / "job_id.txt"
    )

    if not job_id_file.is_file():
        raise ValueError(
            "Latest preprocessing submission has no "
            "job_id.txt and is therefore not a completed "
            "submitted workflow"
        )

    job_id = job_id_file.read_text(
        encoding="utf-8"
    ).strip()

    if not re.fullmatch(
        r"[0-9]+",
        job_id,
    ):
        raise ValueError(
            f"Invalid preprocessing job id: {job_id!r}"
        )

    job_dir = (
        submission
        / f"job_{job_id}"
    )

    status_file = (
        job_dir
        / "job_status.tsv"
    )

    if not status_file.is_file():
        raise ValueError(
            "Preprocessing job status is missing: "
            f"{status_file}"
        )

    status = read_job_status(
        status_file
    )

    if (
        status["workflow_exit_status"] != 0
        or status["exit_status"] != 0
    ):
        raise ValueError(
            "Latest preprocessing submission did not "
            f"complete successfully: {status}"
        )

    project = (
        submission
        / "project"
    ).resolve()

    if not project.is_dir():
        raise ValueError(
            "Preprocessing project directory is missing: "
            f"{project}"
        )

    output_manifest = (
        job_dir
        / "output.sha256"
    )

    if not output_manifest.is_file():
        raise ValueError(
            "Successful preprocessing job is missing "
            f"output.sha256: {output_manifest}"
        )

    verification = subprocess.run(
        [
            "sha256sum",
            "--check",
            str(output_manifest),
        ],
        cwd=project,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    if verification.returncode != 0:
        raise ValueError(
            "Preprocessing output SHA-256 validation "
            "failed:\n"
            + verification.stdout
        )

    artifacts = {}

    for key in (
        "processing_runs",
        "qc_by_fastq",
        "qc_by_job",
    ):
        relative = relative_path(
            config[key]
        )

        absolute = (
            project
            / relative
        )

        if not absolute.is_file():
            raise ValueError(
                f"Required preprocessing artifact is "
                f"missing: {absolute}"
            )

        artifacts[key] = {
            "relative": str(relative),
            "absolute": str(
                absolute.resolve()
            ),
            "sha256": sha256(
                absolute
            ),
        }

    return {
        "schema_version": 1,
        "preprocessing_submission": str(
            submission
        ),
        "preprocessing_project": str(
            project
        ),
        "preprocessing_job_id": job_id,
        "preprocessing_job_directory": str(
            job_dir
        ),
        "job_status": status,
        "job_status_sha256": sha256(
            status_file
        ),
        "output_manifest": str(
            output_manifest
        ),
        "output_manifest_sha256": sha256(
            output_manifest
        ),
        "output_manifest_verified": True,
        "artifacts": artifacts,
    }


def create_alignment_plan(
    root,
    upstream,
    output,
):
    root = Path(root).resolve()
    output = Path(output)

    config = json.loads(
        (
            root
            / "config"
            / "chipseq_alignment.json"
        ).read_text(
            encoding="utf-8"
        )
    )

    artifacts = upstream[
        "artifacts"
    ]

    command = [
        sys.executable,
        str(
            root
            / "workflow"
            / "scripts"
            / "prepare_chipseq_alignment_inputs.py"
        ),
        "--project-root",
        upstream[
            "preprocessing_project"
        ],
        "--processing-runs",
        artifacts[
            "processing_runs"
        ]["absolute"],
        "--qc-by-fastq",
        artifacts[
            "qc_by_fastq"
        ]["absolute"],
        "--qc-by-job",
        artifacts[
            "qc_by_job"
        ]["absolute"],
        "--output",
        str(output),
    ]

    subprocess.run(
        command,
        cwd=root,
        check=True,
    )

    plan = json.loads(
        output.read_text(
            encoding="utf-8"
        )
    )

    if plan.get("schema_version") != 1:
        raise ValueError(
            "Unexpected generated alignment-plan schema"
        )

    if plan.get("run_count", 0) < 1:
        raise ValueError(
            "Generated alignment plan is empty"
        )

    if (
        plan.get("library_layout")
        != "SINGLE"
    ):
        raise ValueError(
            "Operational alignment currently supports "
            "SINGLE only"
        )

    if (
        plan.get("instrument_platform")
        != "ILLUMINA"
    ):
        raise ValueError(
            "Operational alignment currently supports "
            "ILLUMINA only"
        )

    expected_project = str(
        Path(
            upstream[
                "preprocessing_project"
            ]
        ).resolve()
    )

    if (
        str(
            Path(
                plan["source_root"]
            ).resolve()
        )
        != expected_project
    ):
        raise ValueError(
            "Alignment plan source_root does not match "
            "validated preprocessing project"
        )

    if config[
        "scope"
    ] != "dynamic_single_end_alignment_and_qc":
        raise ValueError(
            "Unexpected operational alignment scope"
        )

    return plan


def main():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    mode = parser.add_mutually_exclusive_group(
        required=True
    )

    mode.add_argument(
        "--check",
        action="store_true",
    )

    mode.add_argument(
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
            f"Unexpected branch: {branch}"
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
                / "chipseq_alignment.sbatch"
            ),
        ],
        check=True,
    )

    upstream = (
        resolve_preprocessing_upstream(
            ROOT
        )
    )

    parent = (
        ROOT
        / "results"
        / "chipseq"
        / "alignment"
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
        target = (
            project
            / name
        )

        target.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        shutil.copy2(
            ROOT / name,
            target,
        )

    plan_path = (
        project
        / "config"
        / "chipseq_alignment_inputs.json"
    )

    plan = create_alignment_plan(
        ROOT,
        upstream,
        plan_path,
    )

    (
        submission
        / "upstream_preprocessing.json"
    ).write_text(
        json.dumps(
            {
                **upstream,
                "alignment_run_count": (
                    plan["run_count"]
                ),
                "alignment_role_counts": (
                    plan["role_counts"]
                ),
                "alignment_plan_sha256": (
                    sha256(plan_path)
                ),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    (
        submission
        / "base_commit.txt"
    ).write_bytes(
        subprocess.check_output(
            [
                "git",
                "rev-parse",
                "HEAD",
            ],
            cwd=ROOT,
        )
    )

    (
        submission
        / "git_status.txt"
    ).write_bytes(
        subprocess.check_output(
            [
                "git",
                "status",
                "--short",
            ],
            cwd=ROOT,
        )
    )

    (
        submission
        / "source_changes.patch"
    ).write_bytes(
        subprocess.check_output(
            [
                "git",
                "diff",
                "--binary",
                "--",
                *FILES,
            ],
            cwd=ROOT,
        )
    )

    source_hashes = []

    for path in sorted(
        project.rglob("*")
    ):
        if path.is_file():
            source_hashes.append(
                f"{sha256(path)}  "
                f"{path.relative_to(project).as_posix()}\n"
            )

    (
        submission
        / "source.sha256"
    ).write_text(
        "".join(
            source_hashes
        ),
        encoding="utf-8",
    )

    source_validation = subprocess.run(
        [
            "sha256sum",
            "--check",
            str(
                submission
                / "source.sha256"
            ),
        ],
        cwd=project,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    (
        submission
        / "source_validation.log"
    ).write_text(
        source_validation.stdout,
        encoding="utf-8",
    )

    if source_validation.returncode != 0:
        raise ValueError(
            "Alignment source snapshot SHA validation "
            "failed"
        )

    subprocess.run(
        [
            sys.executable,
            "-m",
            "py_compile",
            "workflow/scripts/"
            "chipseq_alignment_support.py",
            "workflow/scripts/"
            "prepare_chipseq_alignment_inputs.py",
            "workflow/slurm/"
            "submit_chipseq_alignment.py",
        ],
        cwd=project,
        check=True,
    )

    print(
        "[OK] Dynamic ChIP alignment source snapshot "
        "and preprocessing upstream validated.",
        flush=True,
    )

    print(
        "ALIGNMENT RUNS: "
        f"{plan['run_count']} "
        f"(IP={plan['role_counts']['ip']}, "
        f"Input={plan['role_counts']['input']})",
        flush=True,
    )

    print(
        f"UPSTREAM PREPROCESSING: "
        f"{upstream['preprocessing_submission']}",
        flush=True,
    )

    print(
        f"CHECK DIRECTORY: {submission}",
        flush=True,
    )

    if args.check:
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
                / "chipseq_alignment.sbatch"
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
        str(submission) + "\n",
        encoding="utf-8",
    )

    print(
        f"JOB SUBMITTED: {job}",
        flush=True,
    )

    print(
        f"LOGS: {submission}",
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
