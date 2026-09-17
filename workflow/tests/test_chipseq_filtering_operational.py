#!/usr/bin/env python3

import hashlib
import importlib.util
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path


ROOT = Path(
    __file__
).resolve().parents[2]

SUPPORT_PATH = (
    ROOT
    / "workflow"
    / "scripts"
    / "chipseq_filtering_support.py"
)

CONFIG_PATH = (
    ROOT
    / "config"
    / "chipseq_filtering.json"
)

SUBMITTER_PATH = (
    ROOT
    / "workflow"
    / "slurm"
    / "submit_chipseq_filtering.py"
)

SBATCH_PATH = (
    ROOT
    / "workflow"
    / "slurm"
    / "chipseq_filtering.sbatch"
)


SPEC = importlib.util.spec_from_file_location(
    "chipseq_filtering_support",
    SUPPORT_PATH,
)

S = importlib.util.module_from_spec(
    SPEC
)

SPEC.loader.exec_module(
    S
)


def write_json(
    path,
    obj,
):
    path = Path(
        path
    )

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    path.write_text(
        json.dumps(
            obj,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def sha(
    path,
):
    return hashlib.sha256(
        Path(
            path
        ).read_bytes()
    ).hexdigest()


class FilteringOperationalTests(
    unittest.TestCase
):

    def setUp(self):

        self.temp = (
            tempfile.TemporaryDirectory()
        )

        self.root = Path(
            self.temp.name
        )

        (
            self.root
            / "config"
        ).mkdir(
            parents=True,
            exist_ok=True,
        )

        shutil.copyfile(
            CONFIG_PATH,
            self.root
            / S.CONFIG,
        )

    def tearDown(self):
        self.temp.cleanup()

    def make_alignment(
        self,
        *,
        stage="completed",
        exit_status="0",
        missing_bam=None,
        qc_role_override=None,
    ):

        parent = (
            self.root
            / "results"
            / "chipseq"
            / "alignment"
            / "slurm"
        )

        submission = (
            parent
            / "submission_TEST001"
        )

        job_id = "999001"

        job = (
            submission
            / f"job_{job_id}"
        )

        job.mkdir(
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
            job
            / "job_status.tsv"
        ).write_text(
            f"stage\t{stage}\n"
            f"exit_status\t{exit_status}\n",
            encoding="utf-8",
        )

        (
            job
            / "output_copy_validation.log"
        ).write_text(
            "[OK] Outputs copied and SHA-256 verified: "
            + str(
                job
                / "outputs"
            )
            + "\n",
            encoding="utf-8",
        )

        runs = [
            (
                "SRR100",
                "input",
            ),
            (
                "SRR101",
                "ip",
            ),
            (
                "SRR102",
                "ip",
            ),
        ]

        input_rows = []

        for run, role in runs:

            input_rows.append({
                "run_accession": run,
                "role": role,
                "reads": 100,
            })

        write_json(
            job
            / "outputs"
            / "input_provenance.json",
            {
                "schema_version": 1,
                "run_count": 3,
                "role_counts": {
                    "ip": 2,
                    "input": 1,
                },
                "library_layout": "SINGLE",
                "instrument_platform": "ILLUMINA",
                "runs": input_rows,
            },
        )

        write_json(
            job
            / "outputs"
            / "alignment_parameters.json",
            {
                "schema_version": 1,
                "diagnostic_mapq": 30,
                "mitochondrial_accession": "NC_007936.1",
            },
        )

        write_json(
            job
            / "outputs"
            / "reference"
            / "reference_provenance.json",
            {
                "mitochondrial_accession": "NC_007936.1",
            },
        )

        (
            job
            / "outputs"
            / "reference"
            / "genome_plus_mt.fa.fai"
        ).write_text(
            "chr1\t1000\t0\t0\t0\n"
            "NC_007936.1\t16284\t0\t0\t0\n",
            encoding="utf-8",
        )

        (
            job
            / "outputs"
            / "alignment_qc.tsv"
        ).write_text(
            "run_accession\trole\n"
            "SRR100\tinput\n"
            "SRR101\tip\n"
            "SRR102\tip\n",
            encoding="utf-8",
        )

        write_json(
            job
            / "outputs"
            / "verified_fastq_inputs.json",
            {
                "schema_version": 1,
                "run_count": 3,
            },
        )

        for run, role in runs:

            effective_role = role

            if (
                qc_role_override
                and run
                == qc_role_override[0]
            ):
                effective_role = (
                    qc_role_override[1]
                )

            write_json(
                job
                / "outputs"
                / run
                / "alignment_qc.json",
                {
                    "run_accession": run,
                    "role": effective_role,
                    "diagnostic_mapq": 30,
                    "alignment_records": 100,
                    "nonprimary_records": 0,
                    "primary_reads": 100,
                    "mapped_reads": 90,
                    "nuclear_mapq_ge_threshold": 80,
                    "duplicate_marking": "not_performed",
                    "bam_filtering": "not_performed",
                },
            )

            bam = (
                job
                / "outputs"
                / run
                / "raw.sorted.bam"
            )

            bam.parent.mkdir(
                parents=True,
                exist_ok=True,
            )

            bam.write_bytes(
                (
                    f"synthetic-bam-{run}\n"
                ).encode()
            )

            (
                job
                / "outputs"
                / run
                / "raw.sorted.bam.csi"
            ).write_bytes(
                (
                    f"synthetic-index-{run}\n"
                ).encode()
            )

        manifest_lines = []

        for path in sorted(
            (
                job
                / "outputs"
            ).rglob("*")
        ):

            if not path.is_file():
                continue

            relative = (
                Path("outputs")
                / path.relative_to(
                    job
                    / "outputs"
                )
            )

            manifest_lines.append(
                f"{sha(path)}  {relative.as_posix()}\n"
            )

        (
            job
            / "output.sha256"
        ).write_text(
            "".join(
                manifest_lines
            ),
            encoding="utf-8",
        )

        if missing_bam:

            (
                job
                / "outputs"
                / missing_bam
                / "raw.sorted.bam"
            ).unlink()

        parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        (
            parent
            / "latest_submission.txt"
        ).write_text(
            str(
                submission.resolve()
            )
            + "\n",
            encoding="utf-8",
        )

        return (
            submission,
            job,
        )

    def test_successful_alignment_builds_dynamic_filtering_plan(self):

        submission, job = (
            self.make_alignment()
        )

        plan = S.prepare(
            self.root
        )

        self.assertEqual(
            plan[
                "alignment_submission"
            ],
            str(
                submission.resolve()
            ),
        )

        self.assertEqual(
            plan[
                "source_root"
            ],
            str(
                job.resolve()
            ),
        )

        self.assertEqual(
            plan[
                "run_count"
            ],
            3,
        )

        self.assertEqual(
            plan[
                "role_counts"
            ],
            {
                "ip": 2,
                "input": 1,
            },
        )

        self.assertEqual(
            plan[
                "library_layout"
            ],
            "SINGLE",
        )

        self.assertEqual(
            plan[
                "instrument_platform"
            ],
            "ILLUMINA",
        )

        self.assertTrue(
            plan[
                "full_bam_sha256_verification_deferred_to_compute"
            ]
        )

        self.assertFalse(
            plan[
                "upstream_raw_bam_cleanup_authorized"
            ]
        )

        self.assertEqual(
            len(
                plan[
                    "runs"
                ]
            ),
            3,
        )

    def test_missing_latest_alignment_pointer_fails_closed(self):

        with self.assertRaisesRegex(
            ValueError,
            "No submitted dynamic ChIP alignment",
        ):

            S.prepare(
                self.root
            )

    def test_noncompleted_alignment_fails_closed(self):

        self.make_alignment(
            stage="alignment",
            exit_status="1",
        )

        with self.assertRaisesRegex(
            ValueError,
            "did not complete successfully",
        ):

            S.prepare(
                self.root
            )

    def test_missing_raw_bam_fails_closed(self):

        self.make_alignment(
            missing_bam="SRR101",
        )

        with self.assertRaisesRegex(
            ValueError,
            "Missing upstream alignment output",
        ):

            S.prepare(
                self.root
            )

    def test_qc_role_mismatch_fails_closed(self):

        self.make_alignment(
            qc_role_override=(
                "SRR101",
                "input",
            ),
        )

        with self.assertRaisesRegex(
            ValueError,
            "Unexpected upstream SINGLE alignment QC",
        ):

            S.prepare(
                self.root
            )

    def test_operational_layer_contains_no_pilot_dependencies(self):

        combined = "\n".join([
            CONFIG_PATH.read_text(
                encoding="utf-8"
            ),
            SUPPORT_PATH.read_text(
                encoding="utf-8"
            ),
            SUBMITTER_PATH.read_text(
                encoding="utf-8"
            ),
            SBATCH_PATH.read_text(
                encoding="utf-8"
            ),
        ])

        forbidden = [
            "chipseq_pilot.json",
            "chipseq_experimental_eligibility",
            "chipseq-metadata-audit",
            "alignment_evidence",
            "single_end_pilot",
            "Expected one input and one distinct IP",
        ]

        for token in forbidden:

            self.assertNotIn(
                token,
                combined,
            )


if __name__ == "__main__":
    unittest.main()
