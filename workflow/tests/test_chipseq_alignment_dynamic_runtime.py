#!/usr/bin/env python3

import hashlib
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path


ROOT = Path(
    __file__
).resolve().parents[2]

SUPPORT = (
    ROOT
    / "workflow"
    / "scripts"
    / "chipseq_alignment_support.py"
)

SNAKEFILE = (
    ROOT
    / "workflow"
    / "rules"
    / "chipseq_alignment.smk"
)

SPEC = importlib.util.spec_from_file_location(
    "chipseq_alignment_support",
    SUPPORT,
)

MODULE = importlib.util.module_from_spec(
    SPEC
)

SPEC.loader.exec_module(
    MODULE
)


def digest(path):
    return hashlib.sha256(
        Path(path).read_bytes()
    ).hexdigest()


class DynamicAlignmentRuntimeTests(
    unittest.TestCase
):

    def test_snakefile_no_longer_requires_exactly_two_runs(self):
        text = SNAKEFILE.read_text(
            encoding="utf-8"
        )

        self.assertNotIn(
            "Expected exactly two distinct pilot runs",
            text,
        )

        self.assertNotIn(
            "rule chipseq_align_pilot_run:",
            text,
        )

        self.assertIn(
            "rule chipseq_align_run:",
            text,
        )

        self.assertNotIn(
            "chipseq_eligibility.smk",
            text,
        )

        self.assertNotIn(
            "CHIP_ELIGIBILITY_",
            text,
        )

        self.assertNotIn(
            "chipseq_pilot.json",
            text,
        )

        self.assertNotIn(
            "pilot_eligibility.json",
            text,
        )

        self.assertIn(
            'ALIGNMENT_PLAN = "config/chipseq_alignment_inputs.json"',
            text,
        )

    def test_stage_copies_multiple_runs_and_verifies_hashes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)

            source_root = (
                root
                / "preprocessing"
            )

            source_root.mkdir()

            rows = []

            for index, role in (
                (100, "input"),
                (101, "ip"),
                (102, "ip"),
            ):
                run = f"SRR{index}"

                relative = (
                    Path("results")
                    / "preprocessing"
                    / "fastp"
                    / run
                    / "SINGLE"
                    / f"{run}.fastq.gz"
                )

                source = (
                    source_root
                    / relative
                )

                source.parent.mkdir(
                    parents=True,
                    exist_ok=True,
                )

                source.write_bytes(
                    f"FASTQ-{run}\n".encode()
                )

                rows.append({
                    "run_accession": run,
                    "role": role,
                    "source_relative": str(relative),
                    "destination": f"inputs/{run}.fastq.gz",
                    "reads": 100,
                    "bytes": source.stat().st_size,
                    "sha256": digest(source),
                    "library_layout": "SINGLE",
                    "instrument_platform": "ILLUMINA",
                })

            (
                root
                / "config"
            ).mkdir()

            (
                root
                / "config"
                / "chipseq_alignment_inputs.json"
            ).write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "source_root": str(
                            source_root
                        ),
                        "runs": rows,
                    }
                )
                + "\n",
                encoding="utf-8",
            )

            original = Path.cwd()

            try:
                os.chdir(root)

                MODULE.stage()

            finally:
                os.chdir(original)

            verified = json.loads(
                (
                    root
                    / "inputs"
                    / "verified.json"
                ).read_text(
                    encoding="utf-8"
                )
            )

            self.assertEqual(
                verified["run_count"],
                3,
            )

            self.assertEqual(
                {
                    row["run_accession"]
                    for row in verified["runs"]
                },
                {
                    "SRR100",
                    "SRR101",
                    "SRR102",
                },
            )

            for row in rows:
                target = (
                    root
                    / row["destination"]
                )

                self.assertTrue(
                    target.is_file()
                )

                self.assertEqual(
                    digest(target),
                    row["sha256"],
                )

    def test_stage_rejects_duplicate_run(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)

            (
                root
                / "config"
            ).mkdir()

            (
                root
                / "config"
                / "chipseq_alignment_inputs.json"
            ).write_text(
                json.dumps({
                    "schema_version": 1,
                    "source_root": str(root),
                    "runs": [
                        {
                            "run_accession": "SRR100",
                            "role": "input",
                            "source_relative": "a.fastq.gz",
                            "destination": "inputs/SRR100.fastq.gz",
                            "reads": 1,
                            "bytes": 1,
                            "sha256": "0" * 64,
                        },
                        {
                            "run_accession": "SRR100",
                            "role": "ip",
                            "source_relative": "b.fastq.gz",
                            "destination": "inputs/SRR100.fastq.gz",
                            "reads": 1,
                            "bytes": 1,
                            "sha256": "0" * 64,
                        },
                    ],
                })
                + "\n",
                encoding="utf-8",
            )

            original = Path.cwd()

            try:
                os.chdir(root)

                with self.assertRaises(
                    ValueError
                ):
                    MODULE.stage()

            finally:
                os.chdir(original)

    def test_stage_rejects_paired_layout(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)

            (
                root
                / "config"
            ).mkdir()

            (
                root
                / "config"
                / "chipseq_alignment_inputs.json"
            ).write_text(
                json.dumps({
                    "schema_version": 1,
                    "source_root": str(root),
                    "runs": [
                        {
                            "run_accession": "SRR100",
                            "role": "input",
                            "source_relative": "a.fastq.gz",
                            "destination": "inputs/SRR100.fastq.gz",
                            "reads": 1,
                            "bytes": 1,
                            "sha256": "0" * 64,
                            "library_layout": "PAIRED",
                            "instrument_platform": "ILLUMINA",
                        },
                    ],
                })
                + "\n",
                encoding="utf-8",
            )

            original = Path.cwd()

            try:
                os.chdir(root)

                with self.assertRaises(
                    ValueError
                ):
                    MODULE.stage()

            finally:
                os.chdir(original)


if __name__ == "__main__":
    unittest.main()
