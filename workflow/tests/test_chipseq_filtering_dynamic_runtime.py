#!/usr/bin/env python3

import importlib.util
import json
import os
import shutil
import tempfile
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(
    __file__
).resolve().parents[2]

SUPPORT = (
    ROOT
    / "workflow"
    / "scripts"
    / "chipseq_filtering_support.py"
)

SNAKEFILE = (
    ROOT
    / "workflow"
    / "rules"
    / "chipseq_filtering.smk"
)

CONFIG = (
    ROOT
    / "config"
    / "chipseq_filtering.json"
)

SPEC = importlib.util.spec_from_file_location(
    "chipseq_filtering_support",
    SUPPORT,
)

S = importlib.util.module_from_spec(
    SPEC
)

SPEC.loader.exec_module(
    S
)


def make_runs():
    rows = []

    for index in range(22):

        run = (
            f"SRR{50000000 + index}"
        )

        role = (
            "input"
            if index < 4
            else "ip"
        )

        rows.append({
            "run_accession": run,
            "role": role,
            "source_relative": (
                f"outputs/{run}/raw.sorted.bam"
            ),
            "sha256": (
                f"{index + 1:064x}"
            ),
            "bytes": 1000 + index,
            "input_reads": 1000,
            "mapped_reads": 900,
            "nuclear_mapq_ge_threshold": 800,
        })

    return rows


class DynamicFilteringRuntimeTests(
    unittest.TestCase
):

    def setUp(self):
        self.previous = Path.cwd()

        self.temporary = (
            tempfile.TemporaryDirectory()
        )

        os.chdir(
            self.temporary.name
        )

        Path(
            "config"
        ).mkdir()

        shutil.copyfile(
            CONFIG,
            S.CONFIG,
        )

    def tearDown(self):
        os.chdir(
            self.previous
        )

        self.temporary.cleanup()

    def write_plan(
        self,
        rows=None,
    ):
        if rows is None:
            rows = make_runs()

        S.dump(
            S.PLAN,
            {
                "schema_version": 1,
                "source_root": "/synthetic/alignment/job",
                "library_layout": "SINGLE",
                "instrument_platform": "ILLUMINA",
                "configuration_sha256": S.sha(
                    S.CONFIG
                ),
                "runs": rows,
            },
        )

    def test_item_for_accepts_dynamic_22_run_cohort(self):
        rows = make_runs()

        self.write_plan(
            rows
        )

        observed = []

        for expected in rows:
            _, item = S.item_for(
                expected[
                    "run_accession"
                ]
            )

            self.assertEqual(
                item[
                    "run_accession"
                ],
                expected[
                    "run_accession"
                ],
            )

            observed.append(
                item["role"]
            )

        self.assertEqual(
            Counter(observed),
            {
                "ip": 18,
                "input": 4,
            },
        )

    def test_duplicate_run_is_rejected_by_runtime_lookup(self):
        rows = make_runs()

        rows.append(
            dict(
                rows[0]
            )
        )

        self.write_plan(
            rows
        )

        with self.assertRaisesRegex(
            ValueError,
            "Unknown or duplicate filtering run",
        ):
            S.item_for(
                rows[0][
                    "run_accession"
                ]
            )

    def test_validated_scientific_policy_is_preserved(self):
        cfg = S.configuration()

        self.assertEqual(
            cfg[
                "min_mapq"
            ],
            30,
        )

        self.assertEqual(
            cfg[
                "exclude_flags"
            ],
            3844,
        )

        self.assertIs(
            cfg[
                "exclude_mapq_255"
            ],
            True,
        )

        self.assertEqual(
            cfg[
                "mitochondrial_accession"
            ],
            "NC_007936.1",
        )

        self.assertIs(
            cfg[
                "optical_duplicate_detection"
            ],
            False,
        )

        self.assertEqual(
            cfg[
                "duplicate_scoring_strategy"
            ],
            "SUM_OF_BASE_QUALITIES",
        )

        self.assertEqual(
            cfg["threads"],
            4,
        )

        self.assertEqual(
            cfg[
                "picard_heap_mb"
            ],
            12000,
        )

    def test_filtering_snakefile_runtime_is_pilot_independent(self):
        text = SNAKEFILE.read_text(
            encoding="utf-8"
        )

        forbidden = [
            "Expected two distinct pilot runs",
            "chipseq_eligibility.smk",
            "CHIP_ELIGIBILITY_",
            "chipseq_pilot.json",
            "pilot_eligibility.json",
        ]

        for token in forbidden:
            self.assertNotIn(
                token,
                text,
            )

        required = [
            'FILTER_INPUTS.get("library_layout") != "SINGLE"',
            'FILTER_INPUTS.get("instrument_platform") != "ILLUMINA"',
            'ROLE_COUNTS["ip"] < 1',
            'ROLE_COUNTS["input"] < 1',
            'f"outputs/{run}/raw.sorted.bam"',
            "rule chipseq_mark_duplicates:",
            "rule chipseq_filter_bam:",
            "rule chipseq_filtered_qc:",
        ]

        for token in required:
            self.assertIn(
                token,
                text,
            )


if __name__ == "__main__":
    unittest.main()
