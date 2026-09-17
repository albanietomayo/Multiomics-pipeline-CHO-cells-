#!/usr/bin/env python3

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(
    __file__
).resolve().parents[2]

SCRIPT = (
    ROOT
    / "workflow"
    / "scripts"
    / "prepare_chipseq_alignment_inputs.py"
)

SPEC = importlib.util.spec_from_file_location(
    "alignment_inputs",
    SCRIPT,
)

MODULE = importlib.util.module_from_spec(
    SPEC
)

SPEC.loader.exec_module(
    MODULE
)


PROCESSING_FIELDS = [
    "run_accession",
    "library_role",
    "study_accession",
    "declared_target",
    "library_layout",
    "instrument_platform",
    "ready_analysis_count",
    "analysis_ids",
    "fastq_file_count",
]

FASTQ_FIELDS = [
    "study_accession",
    "run_accession",
    "omics",
    "fastq_role",
    "processed_fastq",
    "post_total_sequences",
]

JOB_FIELDS = [
    "run_accession",
    "omics",
    "fastq_role",
    "reads_after",
]


def write_tsv(
    path,
    fields,
    rows,
):
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(
            rows
        )


class AlignmentInputTests(
    unittest.TestCase
):

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()

        self.root = Path(
            self.temp.name
        )

        self.processing = (
            self.root
            / "processing.tsv"
        )

        self.fastq_qc = (
            self.root
            / "fastq.tsv"
        )

        self.job_qc = (
            self.root
            / "job.tsv"
        )

    def tearDown(self):
        self.temp.cleanup()

    def create_fixture(
        self,
        processing=None,
        fastq=None,
        jobs=None,
    ):
        if processing is None:
            processing = [
                {
                    "run_accession": "SRR100",
                    "library_role": "input",
                    "study_accession": "PRJ1",
                    "declared_target": "",
                    "library_layout": "SINGLE",
                    "instrument_platform": "ILLUMINA",
                    "ready_analysis_count": "2",
                    "analysis_ids": (
                        "SRR101;SRR102"
                    ),
                    "fastq_file_count": "1",
                },
                {
                    "run_accession": "SRR101",
                    "library_role": "ip",
                    "study_accession": "PRJ1",
                    "declared_target": "H3K4me3",
                    "library_layout": "SINGLE",
                    "instrument_platform": "ILLUMINA",
                    "ready_analysis_count": "1",
                    "analysis_ids": "SRR101",
                    "fastq_file_count": "1",
                },
            ]

        if fastq is None:
            fastq = []

            for row in processing:
                run = row[
                    "run_accession"
                ]

                relative = (
                    Path(
                        "results/preprocessing/fastp"
                    )
                    / run
                    / "SINGLE"
                    / f"{run}.fastq.gz"
                )

                target = (
                    self.root
                    / relative
                )

                target.parent.mkdir(
                    parents=True,
                    exist_ok=True,
                )

                target.write_bytes(
                    f"FASTQ-{run}\n".encode()
                )

                fastq.append({
                    "study_accession": row[
                        "study_accession"
                    ],
                    "run_accession": run,
                    "omics": "ChIP-seq",
                    "fastq_role": "SINGLE",
                    "processed_fastq": str(
                        relative
                    ),
                    "post_total_sequences": "100",
                })

        if jobs is None:
            jobs = [
                {
                    "run_accession": row[
                        "run_accession"
                    ],
                    "omics": "ChIP-seq",
                    "fastq_role": "SINGLE",
                    "reads_after": "100",
                }
                for row in processing
            ]

        write_tsv(
            self.processing,
            PROCESSING_FIELDS,
            processing,
        )

        write_tsv(
            self.fastq_qc,
            FASTQ_FIELDS,
            fastq,
        )

        write_tsv(
            self.job_qc,
            JOB_FIELDS,
            jobs,
        )

        return (
            processing,
            fastq,
            jobs,
        )

    def build(self):
        return MODULE.build_plan(
            project_root=self.root,
            processing_runs=self.processing,
            qc_by_fastq=self.fastq_qc,
            qc_by_job=self.job_qc,
        )

    def test_valid_dynamic_plan(self):
        self.create_fixture()

        plan = self.build()

        self.assertEqual(
            plan["run_count"],
            2,
        )

        self.assertEqual(
            plan["role_counts"],
            {
                "ip": 1,
                "input": 1,
            },
        )

        self.assertEqual(
            [
                row["run_accession"]
                for row in plan["runs"]
            ],
            [
                "SRR100",
                "SRR101",
            ],
        )

    def test_pilot_accessions_are_not_special_cases(self):
        processing = [
            {
                "run_accession": "SRR20770287",
                "library_role": "input",
                "study_accession": "PRJNA865478",
                "declared_target": "",
                "library_layout": "SINGLE",
                "instrument_platform": "ILLUMINA",
                "ready_analysis_count": "1",
                "analysis_ids": "SRR20770297",
                "fastq_file_count": "1",
            },
            {
                "run_accession": "SRR20770297",
                "library_role": "ip",
                "study_accession": "PRJNA865478",
                "declared_target": "H3K4me3",
                "library_layout": "SINGLE",
                "instrument_platform": "ILLUMINA",
                "ready_analysis_count": "1",
                "analysis_ids": "SRR20770297",
                "fastq_file_count": "1",
            },
        ]

        self.create_fixture(
            processing=processing
        )

        plan = self.build()

        self.assertEqual(
            {
                row["run_accession"]
                for row in plan["runs"]
            },
            {
                "SRR20770287",
                "SRR20770297",
            },
        )

    def test_missing_fastq_qc_run_fails_closed(self):
        processing, fastq, jobs = (
            self.create_fixture()
        )

        write_tsv(
            self.fastq_qc,
            FASTQ_FIELDS,
            fastq[:-1],
        )

        with self.assertRaises(
            ValueError
        ):
            self.build()

    def test_extra_job_qc_run_fails_closed(self):
        processing, fastq, jobs = (
            self.create_fixture()
        )

        jobs.append({
            "run_accession": "SRR999",
            "omics": "ChIP-seq",
            "fastq_role": "SINGLE",
            "reads_after": "100",
        })

        write_tsv(
            self.job_qc,
            JOB_FIELDS,
            jobs,
        )

        with self.assertRaises(
            ValueError
        ):
            self.build()

    def test_paired_layout_requires_review(self):
        processing, fastq, jobs = (
            self.create_fixture()
        )

        processing[1][
            "library_layout"
        ] = "PAIRED"

        write_tsv(
            self.processing,
            PROCESSING_FIELDS,
            processing,
        )

        with self.assertRaises(
            ValueError
        ):
            self.build()

    def test_non_illumina_requires_review(self):
        processing, fastq, jobs = (
            self.create_fixture()
        )

        processing[1][
            "instrument_platform"
        ] = "OTHER"

        write_tsv(
            self.processing,
            PROCESSING_FIELDS,
            processing,
        )

        with self.assertRaises(
            ValueError
        ):
            self.build()

    def test_read_count_mismatch_fails_closed(self):
        processing, fastq, jobs = (
            self.create_fixture()
        )

        jobs[1][
            "reads_after"
        ] = "99"

        write_tsv(
            self.job_qc,
            JOB_FIELDS,
            jobs,
        )

        with self.assertRaises(
            ValueError
        ):
            self.build()

    def test_path_escape_fails_closed(self):
        processing, fastq, jobs = (
            self.create_fixture()
        )

        fastq[1][
            "processed_fastq"
        ] = "../outside.fastq.gz"

        write_tsv(
            self.fastq_qc,
            FASTQ_FIELDS,
            fastq,
        )

        with self.assertRaises(
            ValueError
        ):
            self.build()


if __name__ == "__main__":
    unittest.main()
