#!/usr/bin/env python3

import csv
import importlib.util
import json
import tempfile
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

SCRIPT = (
    ROOT /
    "workflow/scripts/"
    "build_chipseq_peak_calling_plan.py"
)

POLICY = (
    ROOT /
    "config/chipseq_peak_calling_policy.json"
)

REAL_ANALYSIS_PLAN = (
    ROOT / "snapshots/chipseq/incremental_planning_validation_001/"
    "chipseq_analysis_plan.tsv"
)


spec = importlib.util.spec_from_file_location(
    "peak_plan",
    SCRIPT,
)

peak_plan = importlib.util.module_from_spec(spec)
spec.loader.exec_module(peak_plan)


FIELDS = [
    "analysis_id",
    "study_accession",
    "ip_run_accession",
    "control_run_accession",
    "declared_target",
    "library_layout",
    "analysis_status",
]


def write_plan(path, rows):
    with Path(path).open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=FIELDS,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)


def row(
    analysis_id,
    study,
    target,
    status="ready",
    layout="SINGLE",
    control=None,
):
    return {
        "analysis_id": analysis_id,
        "study_accession": study,
        "ip_run_accession": analysis_id,
        "control_run_accession": (
            control or analysis_id + "_INPUT"
        ),
        "declared_target": target,
        "library_layout": layout,
        "analysis_status": status,
    }


class PeakCallingPlanTests(unittest.TestCase):

    def test_validated_real_cohort_is_inherited_exactly(self):
        result = peak_plan.build_plan(REAL_ANALYSIS_PLAN, POLICY)

        with REAL_ANALYSIS_PLAN.open(newline="", encoding="utf-8") as handle:
            upstream = [
                item for item in csv.DictReader(handle, delimiter="\t")
                if item["analysis_status"] == "ready"
            ]

        self.assertEqual(len(result), 18)
        self.assertEqual(Counter(x["study_accession"] for x in result), {
            "PRJEB9291": 12, "PRJNA865478": 6,
        })
        self.assertEqual(Counter(x["peak_mode"] for x in result), {
            "narrow": 8, "broad": 10,
        })
        self.assertEqual(Counter(x["fragment_size_policy"] for x in result), {
            "fixed": 6, "phantompeakqualtools": 12,
        })
        self.assertEqual(
            len({x[key] for x in result for key in (
                "ip_run_accession", "control_run_accession",
            )}),
            22,
        )

        inherited = {
            x["analysis_id"]: (x["ip_run_accession"], x["control_run_accession"])
            for x in result
        }
        expected = {
            x["analysis_id"]: (x["ip_run_accession"], x["control_run_accession"])
            for x in upstream
        }
        self.assertEqual(inherited, expected)

    def build(
        self,
        rows,
        policy_path=POLICY,
    ):
        with tempfile.TemporaryDirectory() as td:
            plan_path = (
                Path(td) /
                "analysis.tsv"
            )

            write_plan(
                plan_path,
                rows,
            )

            return peak_plan.build_plan(
                plan_path,
                policy_path,
            )

    def test_fixed_prjna_narrow(self):
        result = self.build([
            row(
                "SRR1",
                "PRJNA865478",
                "H3K4me3",
            )
        ])

        self.assertEqual(
            len(result),
            1,
        )

        item = result[0]

        self.assertEqual(
            item["peak_mode"],
            "narrow",
        )

        self.assertEqual(
            item["fragment_size_policy"],
            "fixed",
        )

        self.assertEqual(
            item["fragment_size_bp"],
            "147",
        )

        self.assertEqual(
            item["execution_status"],
            "ready_for_peak_calling",
        )

        self.assertEqual(
            item["broad_cutoff"],
            "",
        )

    def test_fixed_prjna_broad(self):
        result = self.build([
            row(
                "SRR2",
                "PRJNA865478",
                "H3K9me3",
            )
        ])

        item = result[0]

        self.assertEqual(
            item["peak_mode"],
            "broad",
        )

        self.assertEqual(
            item["broad_cutoff"],
            "0.1",
        )

        self.assertEqual(
            item["fragment_size_bp"],
            "147",
        )

    def test_prjeb_requires_fragment_estimation(self):
        result = self.build([
            row(
                "ERR1",
                "PRJEB9291",
                "H3K36me3",
            )
        ])

        item = result[0]

        self.assertEqual(
            item["peak_mode"],
            "broad",
        )

        self.assertEqual(
            item["fragment_size_policy"],
            "phantompeakqualtools",
        )

        self.assertEqual(
            item["fragment_size_bp"],
            "",
        )

        self.assertEqual(
            item["fragment_size_status"],
            "requires_estimation",
        )

        self.assertEqual(
            item["execution_status"],
            "requires_fragment_estimation",
        )

    def test_prjeb_narrow(self):
        result = self.build([
            row(
                "ERR2",
                "PRJEB9291",
                "H3K27ac",
            )
        ])

        self.assertEqual(
            result[0]["peak_mode"],
            "narrow",
        )

    def test_non_ready_rows_not_promoted(self):
        result = self.build([
            row(
                "A",
                "PRJNA865478",
                "H3K4me3",
                status="blocked_no_control",
            ),
            row(
                "B",
                "PRJNA865478",
                "H3K4me3",
                status="review_required",
            ),
        ])

        self.assertEqual(
            result,
            [],
        )

    def test_unknown_study_fails_closed(self):
        with self.assertRaisesRegex(
            ValueError,
            "No peak-calling study policy",
        ):
            self.build([
                row(
                    "X",
                    "PRJUNKNOWN",
                    "H3K4me3",
                )
            ])

    def test_unknown_target_fails_closed(self):
        with self.assertRaisesRegex(
            ValueError,
            "No peak-calling target policy",
        ):
            self.build([
                row(
                    "X",
                    "PRJEB9291",
                    "UNKNOWN",
                )
            ])

    def test_missing_control_fails_closed(self):
        bad = row(
            "X",
            "PRJEB9291",
            "H3K4me3",
        )

        bad["control_run_accession"] = ""

        with self.assertRaisesRegex(
            ValueError,
            "Missing control_run_accession",
        ):
            self.build([bad])

    def test_same_ip_control_fails_closed(self):
        bad = row(
            "X",
            "PRJEB9291",
            "H3K4me3",
            control="X",
        )

        with self.assertRaisesRegex(
            ValueError,
            "IP and Input are identical",
        ):
            self.build([bad])

    def test_unsupported_layout_fails_closed(self):
        with self.assertRaisesRegex(
            ValueError,
            "Unsupported library layout",
        ):
            self.build([
                row(
                    "X",
                    "PRJEB9291",
                    "H3K4me3",
                    layout="PAIRED",
                )
            ])

    def test_duplicate_analysis_fails_closed(self):
        rows = [
            row(
                "X",
                "PRJEB9291",
                "H3K4me3",
            ),
            row(
                "X",
                "PRJEB9291",
                "H3K27ac",
            ),
        ]

        with self.assertRaisesRegex(
            ValueError,
            "Duplicate ready analysis_id",
        ):
            self.build(rows)

    def test_duplicate_ready_ip_fails_closed(self):
        first = row("A", "PRJEB9291", "H3K4me3")
        second = row("B", "PRJEB9291", "H3K27ac")
        second["ip_run_accession"] = first["ip_run_accession"]
        with self.assertRaisesRegex(ValueError, "Duplicate ready IP analysis"):
            self.build([first, second])

    def test_unsafe_identifier_fails_closed(self):
        bad = row("../escape", "PRJEB9291", "H3K4me3")
        with self.assertRaisesRegex(ValueError, "Unsafe analysis_id"):
            self.build([bad])

    def test_non_boolean_policy_value_fails_closed(self):
        policy = json.loads(POLICY.read_text())
        policy["technical_defaults"]["spmr"] = 1
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "policy.json"
            path.write_text(json.dumps(policy))
            with self.assertRaisesRegex(ValueError, "must be a boolean"):
                self.build([row("A", "PRJEB9291", "H3K4me3")], path)

    def test_unsupported_genome_size_policy_fails_closed(self):
        policy = json.loads(POLICY.read_text())
        policy["technical_defaults"]["effective_genome_size_policy"] = "magic"
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "policy.json"
            path.write_text(json.dumps(policy))
            with self.assertRaisesRegex(ValueError, "genome-size policy"):
                self.build([row("A", "PRJEB9291", "H3K4me3")], path)

    def test_invalid_fixed_extsize_fails_closed(self):
        policy = json.loads(
            POLICY.read_text()
        )

        policy[
            "studies"
        ][
            "PRJNA865478"
        ][
            "fragment_size_policy"
        ][
            "extsize"
        ] = 0

        with tempfile.TemporaryDirectory() as td:
            bad_policy = (
                Path(td) /
                "policy.json"
            )

            bad_policy.write_text(
                json.dumps(policy),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError,
                "positive extsize",
            ):
                self.build(
                    [
                        row(
                            "X",
                            "PRJNA865478",
                            "H3K4me3",
                        )
                    ],
                    bad_policy,
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
