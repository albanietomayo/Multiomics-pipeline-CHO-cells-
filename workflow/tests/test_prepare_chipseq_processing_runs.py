#!/usr/bin/env python3

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

SPEC = importlib.util.spec_from_file_location(
    "processing_runs",
    ROOT
    / "workflow/scripts/prepare_chipseq_processing_runs.py",
)

M = importlib.util.module_from_spec(
    SPEC
)

SPEC.loader.exec_module(
    M
)


def sample(
    run,
    study="STUDY1",
    layout="SINGLE",
    omics="ChIP-seq",
    platform="ILLUMINA",
):
    return {
        "run_accession": run,
        "study_accession": study,
        "omics": omics,
        "library_layout": layout,
        "instrument_platform": platform,
        "fastq_ftp": (
            f"ftp.sra.ebi.ac.uk/{run}.fastq.gz"
        ),
        "fastq_bytes": "100",
        "fastq_md5": (
            "0123456789abcdef0123456789abcdef"
        ),
    }


def analysis(
    analysis_id,
    ip,
    control,
    status="ready",
    study="STUDY1",
    target="H3K4me3",
    layout="SINGLE",
):
    return {
        "analysis_id": analysis_id,
        "study_accession": study,
        "ip_run_accession": ip,
        "control_run_accession": control,
        "declared_target": target,
        "library_layout": layout,
        "instrument_platform": "ILLUMINA",
        "analysis_status": status,
    }


class ProcessingRunSetTest(
    unittest.TestCase
):

    def test_shared_input_is_deduplicated(self):
        plan = [
            analysis(
                "SRR10000001",
                "SRR10000001",
                "SRR10000003",
                target="H3K4me3",
            ),
            analysis(
                "SRR10000002",
                "SRR10000002",
                "SRR10000003",
                target="H3K27ac",
            ),
        ]

        samples = [
            sample("SRR10000001"),
            sample("SRR10000002"),
            sample("SRR10000003"),
        ]

        rows, ready = M.build_processing_rows(
            plan,
            samples,
        )

        self.assertEqual(
            ready,
            2,
        )

        self.assertEqual(
            len(rows),
            3,
        )

        by_run = {
            row["run_accession"]: row
            for row in rows
        }

        self.assertEqual(
            by_run[
                "SRR10000003"
            ][
                "library_role"
            ],
            "input",
        )

        self.assertEqual(
            by_run[
                "SRR10000003"
            ][
                "ready_analysis_count"
            ],
            2,
        )

    def test_non_ready_analyses_are_excluded(self):
        plan = [
            analysis(
                "SRR10000001",
                "SRR10000001",
                "SRR10000003",
            ),
            analysis(
                "SRR10000002",
                "SRR10000002",
                "",
                status="blocked_no_control",
            ),
        ]

        samples = [
            sample("SRR10000001"),
            sample("SRR10000002"),
            sample("SRR10000003"),
        ]

        rows, ready = M.build_processing_rows(
            plan,
            samples,
        )

        self.assertEqual(
            ready,
            1,
        )

        self.assertEqual(
            {
                row["run_accession"]
                for row in rows
            },
            {
                "SRR10000001",
                "SRR10000003",
            },
        )

    def test_missing_sample_is_rejected(self):
        plan = [
            analysis(
                "SRR10000001",
                "SRR10000001",
                "SRR10000003",
            )
        ]

        samples = [
            sample("SRR10000001"),
        ]

        with self.assertRaisesRegex(
            ValueError,
            "missing from sample catalogue",
        ):
            M.build_processing_rows(
                plan,
                samples,
            )

    def test_conflicting_role_is_rejected(self):
        plan = [
            analysis(
                "SRR10000001",
                "SRR10000001",
                "SRR10000003",
            ),
            analysis(
                "SRR10000003",
                "SRR10000003",
                "SRR10000004",
                target="H3K27ac",
            ),
        ]

        samples = [
            sample("SRR10000001"),
            sample("SRR10000003"),
            sample("SRR10000004"),
        ]

        with self.assertRaisesRegex(
            ValueError,
            "conflicting roles",
        ):
            M.build_processing_rows(
                plan,
                samples,
            )

    def test_unsupported_layout_is_rejected(self):
        plan = [
            analysis(
                "SRR10000001",
                "SRR10000001",
                "SRR10000003",
                layout="PAIRED",
            )
        ]

        samples = [
            sample(
                "SRR10000001",
                layout="PAIRED",
            ),
            sample(
                "SRR10000003",
                layout="PAIRED",
            ),
        ]

        with self.assertRaisesRegex(
            ValueError,
            "unsupported layout",
        ):
            M.build_processing_rows(
                plan,
                samples,
            )

    def test_ready_analysis_without_control_is_rejected(self):
        plan = [
            analysis(
                "SRR10000001",
                "SRR10000001",
                "",
            )
        ]

        samples = [
            sample("SRR10000001"),
        ]

        with self.assertRaisesRegex(
            ValueError,
            "lacks IP/Input",
        ):
            M.build_processing_rows(
                plan,
                samples,
            )

    def test_no_ready_analyses_is_rejected(self):
        plan = [
            analysis(
                "SRR10000001",
                "SRR10000001",
                "",
                status="review_required",
            )
        ]

        samples = [
            sample("SRR10000001"),
        ]

        with self.assertRaisesRegex(
            ValueError,
            "no ready analyses",
        ):
            M.build_processing_rows(
                plan,
                samples,
            )


if __name__ == "__main__":
    unittest.main()
