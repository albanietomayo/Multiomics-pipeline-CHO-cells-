#!/usr/bin/env python3

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(
    __file__
).resolve().parents[2]

SUBMITTER = (
    ROOT
    / "workflow"
    / "slurm"
    / "submit_chipseq_alignment.py"
)

SUPPORT = (
    ROOT
    / "workflow"
    / "scripts"
    / "chipseq_alignment_support.py"
)

CONFIG = (
    ROOT
    / "config"
    / "chipseq_alignment.json"
)

SPEC = importlib.util.spec_from_file_location(
    "alignment_submitter",
    SUBMITTER,
)

MODULE = importlib.util.module_from_spec(
    SPEC
)

SPEC.loader.exec_module(
    MODULE
)


def sha(path):
    return hashlib.sha256(
        Path(path).read_bytes()
    ).hexdigest()


def build_upstream(
    root,
    *,
    workflow_exit=0,
    final_exit=0,
):
    root = Path(root)

    config_dir = (
        root
        / "config"
    )

    config_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    config = json.loads(
        CONFIG.read_text(
            encoding="utf-8"
        )
    )

    (
        config_dir
        / "chipseq_alignment.json"
    ).write_text(
        json.dumps(
            config,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    slurm = (
        root
        / "results"
        / "chipseq"
        / "preprocessing"
        / "slurm"
    )

    submission = (
        slurm
        / "submission_TEST01"
    )

    project = (
        submission
        / "project"
    )

    job_id = "123456"

    job_dir = (
        submission
        / f"job_{job_id}"
    )

    project.mkdir(
        parents=True,
        exist_ok=True,
    )

    job_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    (
        submission
        / "job_id.txt"
    ).write_text(
        job_id + "\n",
        encoding="utf-8",
    )

    (
        job_dir
        / "job_status.tsv"
    ).write_text(
        (
            "workflow_exit_status\t"
            f"{workflow_exit}\n"
            "exit_status\t"
            f"{final_exit}\n"
        ),
        encoding="utf-8",
    )

    for relative in (
        config["processing_runs"],
        config["qc_by_fastq"],
        config["qc_by_job"],
    ):
        path = (
            project
            / relative
        )

        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        path.write_text(
            "fixture\n",
            encoding="utf-8",
        )

    output_manifest = (
        job_dir
        / "output.sha256"
    )

    lines = []

    for path in sorted(
        (
            project
            / "results"
        ).rglob("*")
    ):
        if path.is_file():
            lines.append(
                f"{sha(path)}  "
                f"{path.relative_to(project).as_posix()}\n"
            )

    output_manifest.write_text(
        "".join(lines),
        encoding="utf-8",
    )

    slurm.mkdir(
        parents=True,
        exist_ok=True,
    )

    (
        slurm
        / "latest_submission.txt"
    ).write_text(
        str(submission.resolve())
        + "\n",
        encoding="utf-8",
    )

    return submission


class OperationalAlignmentTests(
    unittest.TestCase
):

    def test_successful_upstream_resolves(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)

            build_upstream(
                root
            )

            result = (
                MODULE.resolve_preprocessing_upstream(
                    root
                )
            )

            self.assertEqual(
                result[
                    "job_status"
                ],
                {
                    "workflow_exit_status": 0,
                    "exit_status": 0,
                },
            )

            self.assertTrue(
                result[
                    "output_manifest_verified"
                ]
            )

            self.assertEqual(
                set(
                    result[
                        "artifacts"
                    ]
                ),
                {
                    "processing_runs",
                    "qc_by_fastq",
                    "qc_by_job",
                },
            )

    def test_failed_upstream_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)

            build_upstream(
                root,
                workflow_exit=1,
                final_exit=1,
            )

            with self.assertRaises(
                ValueError
            ):
                MODULE.resolve_preprocessing_upstream(
                    root
                )

    def test_missing_latest_pointer_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)

            (
                root
                / "config"
            ).mkdir(
                parents=True,
                exist_ok=True,
            )

            (
                root
                / "config"
                / "chipseq_alignment.json"
            ).write_bytes(
                CONFIG.read_bytes()
            )

            with self.assertRaises(
                ValueError
            ):
                MODULE.resolve_preprocessing_upstream(
                    root
                )

    def test_configuration_has_no_legacy_pilot_upstream(self):
        obj = json.loads(
            CONFIG.read_text(
                encoding="utf-8"
            )
        )

        for key in (
            "source_project",
            "source_checksums",
            "source_job_status",
        ):
            self.assertNotIn(
                key,
                obj,
            )

        self.assertEqual(
            obj["scope"],
            "dynamic_single_end_alignment_and_qc",
        )

    def test_support_has_no_legacy_prepare_or_pilot(self):
        text = SUPPORT.read_text(
            encoding="utf-8"
        )

        self.assertNotIn(
            "def prepare(",
            text,
        )

        self.assertNotIn(
            "chipseq_pilot.json",
            text,
        )

        self.assertNotIn(
            'add_parser("prepare")',
            text,
        )


if __name__ == "__main__":
    unittest.main()
