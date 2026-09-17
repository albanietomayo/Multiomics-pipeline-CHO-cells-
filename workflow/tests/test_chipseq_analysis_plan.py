#!/usr/bin/env python3

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

SPEC = importlib.util.spec_from_file_location(
    "chipseq_analysis_plan",
    ROOT / "workflow/scripts/build_chipseq_analysis_plan.py",
)

M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


class AnalysisPlanTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

        self.policy = self.dir / "policy.json"
        self.policy.write_text(
            json.dumps({
                "schema_version": 1,
                "scope": "incremental_chipseq_analysis_planning",
                "supported_library_layouts": ["SINGLE"],
                "condition_match_fields": list(M.CONDITION_FIELDS),
            }),
            encoding="utf-8",
        )

    def tearDown(self):
        self.tmp.cleanup()

    def dataset(
        self,
        candidate_status="single_candidate_needs_review",
        candidate_runs=None,
        issue_runs=None,
        ip_eligibility="retained_by_rules",
        control_eligibility="retained_by_rules",
        control_time="3",
    ):
        if candidate_runs is None:
            candidate_runs = ["CTRL1"]

        if issue_runs is None:
            issue_runs = []

        runs = [
            {
                "run_accession": "IP1",
                "study_accession": "STUDY1",
                "experiment_accession": "EXP1",
                "sample_accession": "S1",
                "omics": "ChIP-seq",
                "experiment_title": "",
                "sample_title": "",
                "library_name": "",
                "library_layout": "SINGLE",
                "instrument_platform": "ILLUMINA",
            },
            {
                "run_accession": "CTRL1",
                "study_accession": "STUDY1",
                "experiment_accession": "EXP2",
                "sample_accession": "S2",
                "omics": "ChIP-seq",
                "experiment_title": "",
                "sample_title": "",
                "library_name": "",
                "library_layout": "SINGLE",
                "instrument_platform": "ILLUMINA",
            },
        ]

        annotations = [
            {
                "run_accession": "IP1",
                "study_accession": "STUDY1",
                "experiment_accession": "EXP1",
                "sample_accession": "S1",
                "declared_target": "H3K4me3",
                "library_role": "ip",
                "label_status": "resolved_from_labels",
                "target_candidates": '["H3K4me3"]',
                "control_role_candidates": "[]",
            },
            {
                "run_accession": "CTRL1",
                "study_accession": "STUDY1",
                "experiment_accession": "EXP2",
                "sample_accession": "S2",
                "declared_target": "",
                "library_role": "input",
                "label_status": "resolved_from_labels",
                "target_candidates": "[]",
                "control_role_candidates": '["input"]',
            },
        ]

        condition_fields = [
            "run_accession",
            "study_accession",
            "sample_accession",
            "experiment_accession",
            "library_role",
            "declared_target",
            "library_layout",
            "instrument_platform",
            "condition_status",
            "condition_rule",
            "label_context",
            "time_value",
            "time_unit",
            "documented_issue",
            "evidence_json",
        ]

        conditions = [
            {
                "run_accession": "IP1",
                "study_accession": "STUDY1",
                "sample_accession": "S1",
                "experiment_accession": "EXP1",
                "library_role": "ip",
                "declared_target": "H3K4me3",
                "library_layout": "SINGLE",
                "instrument_platform": "ILLUMINA",
                "condition_status": "resolved_from_labels",
                "condition_rule": "title_time_label",
                "label_context": "cho-k1",
                "time_value": "3",
                "time_unit": "days",
                "documented_issue": "false",
                "evidence_json": "[]",
            },
            {
                "run_accession": "CTRL1",
                "study_accession": "STUDY1",
                "sample_accession": "S2",
                "experiment_accession": "EXP2",
                "library_role": "input",
                "declared_target": "",
                "library_layout": "SINGLE",
                "instrument_platform": "ILLUMINA",
                "condition_status": "resolved_from_labels",
                "condition_rule": "title_time_label",
                "label_context": "cho-k1",
                "time_value": control_time,
                "time_unit": "days",
                "documented_issue": "false",
                "evidence_json": "[]",
            },
        ]

        candidates = [{
            "ip_run_accession": "IP1",
            "study_accession": "STUDY1",
            "candidate_status": candidate_status,
            "candidate_count": str(len(candidate_runs)),
            "candidate_input_runs": json.dumps(candidate_runs),
            "candidates_with_documented_issue": json.dumps(issue_runs),
            "pairing_approved": "false",
        }]

        eligibility = {
            "schema_version": 1,
            "scope": "experimental_protocol_eligibility_only",
            "runs": [
                {
                    "run_accession": "IP1",
                    "study_accession": "STUDY1",
                    "library_role": "ip",
                    "declared_target": "H3K4me3",
                    "experimental_eligibility": ip_eligibility,
                },
                {
                    "run_accession": "CTRL1",
                    "study_accession": "STUDY1",
                    "library_role": "input",
                    "declared_target": "",
                    "experimental_eligibility": control_eligibility,
                },
            ],
        }

        paths = {
            "runs": self.dir / "runs.tsv",
            "annotations": self.dir / "annotations.tsv",
            "conditions": self.dir / "conditions.tsv",
            "candidates": self.dir / "candidates.tsv",
            "eligibility": self.dir / "eligibility.json",
        }

        write_tsv(
            paths["runs"],
            list(runs[0]),
            runs,
        )

        write_tsv(
            paths["annotations"],
            list(annotations[0]),
            annotations,
        )

        write_tsv(
            paths["conditions"],
            condition_fields,
            conditions,
        )

        write_tsv(
            paths["candidates"],
            list(candidates[0]),
            candidates,
        )

        paths["eligibility"].write_text(
            json.dumps(eligibility),
            encoding="utf-8",
        )

        return paths

    def build(self, paths):
        return M.build_plan(
            paths["runs"],
            paths["annotations"],
            paths["conditions"],
            paths["candidates"],
            paths["eligibility"],
            self.policy,
        )

    def test_unique_exact_condition_control_is_ready(self):
        plan = self.build(self.dataset())

        self.assertEqual(len(plan), 1)
        self.assertEqual(plan[0]["analysis_status"], "ready")
        self.assertEqual(plan[0]["control_run_accession"], "CTRL1")
        self.assertEqual(
            plan[0]["pairing_rule"],
            "unique_exact_condition_input",
        )

    def test_no_control_is_blocked(self):
        paths = self.dataset(
            candidate_status="no_candidate",
            candidate_runs=[],
        )

        plan = self.build(paths)

        self.assertEqual(
            plan[0]["analysis_status"],
            "blocked_no_control",
        )
        self.assertEqual(
            plan[0]["status_reason"],
            "no_compatible_input_control",
        )

    def test_multiple_controls_require_review(self):
        paths = self.dataset(
            candidate_status="multiple_candidates_needs_review",
            candidate_runs=["CTRL1", "CTRL2"],
        )

        candidates, fields = M.read_tsv(paths["candidates"])
        candidates[0]["candidate_count"] = "2"

        write_tsv(
            paths["candidates"],
            fields,
            candidates,
        )

        plan = self.build(paths)

        self.assertEqual(
            plan[0]["analysis_status"],
            "review_required",
        )

    def test_documented_input_issue_requires_review(self):
        plan = self.build(
            self.dataset(issue_runs=["CTRL1"])
        )

        self.assertEqual(
            plan[0]["analysis_status"],
            "review_required",
        )
        self.assertEqual(
            plan[0]["status_reason"],
            "unique_input_has_documented_issue",
        )

    def test_nonretained_input_requires_review(self):
        plan = self.build(
            self.dataset(
                control_eligibility="review_required",
            )
        )

        self.assertEqual(
            plan[0]["analysis_status"],
            "review_required",
        )
        self.assertEqual(
            plan[0]["status_reason"],
            "unique_input_not_retained_by_eligibility",
        )

    def test_candidate_condition_mismatch_is_rejected(self):
        paths = self.dataset(control_time="7")

        with self.assertRaisesRegex(
            ValueError,
            "exact condition key",
        ):
            self.build(paths)

    def test_new_upstream_classified_exact_pair_enters_plan_automatically(self):
        paths = self.dataset()

        runs, run_fields = M.read_tsv(paths["runs"])
        annotations, annotation_fields = M.read_tsv(
            paths["annotations"]
        )
        conditions, condition_fields = M.read_tsv(
            paths["conditions"]
        )
        candidates, candidate_fields = M.read_tsv(
            paths["candidates"]
        )
        eligibility = M.load_json(paths["eligibility"])

        runs.extend([
            {
                **runs[0],
                "run_accession": "IP_NEW",
                "experiment_accession": "EXP_NEW_IP",
                "sample_accession": "S_NEW_IP",
            },
            {
                **runs[1],
                "run_accession": "CTRL_NEW",
                "experiment_accession": "EXP_NEW_CTRL",
                "sample_accession": "S_NEW_CTRL",
            },
        ])

        annotations.extend([
            {
                **annotations[0],
                "run_accession": "IP_NEW",
                "experiment_accession": "EXP_NEW_IP",
                "sample_accession": "S_NEW_IP",
            },
            {
                **annotations[1],
                "run_accession": "CTRL_NEW",
                "experiment_accession": "EXP_NEW_CTRL",
                "sample_accession": "S_NEW_CTRL",
            },
        ])

        conditions.extend([
            {
                **conditions[0],
                "run_accession": "IP_NEW",
                "experiment_accession": "EXP_NEW_IP",
                "sample_accession": "S_NEW_IP",
                "time_value": "7",
            },
            {
                **conditions[1],
                "run_accession": "CTRL_NEW",
                "experiment_accession": "EXP_NEW_CTRL",
                "sample_accession": "S_NEW_CTRL",
                "time_value": "7",
            },
        ])

        candidates.append({
            **candidates[0],
            "ip_run_accession": "IP_NEW",
            "candidate_count": "1",
            "candidate_input_runs": '["CTRL_NEW"]',
            "candidates_with_documented_issue": "[]",
        })

        eligibility["runs"].extend([
            {
                **eligibility["runs"][0],
                "run_accession": "IP_NEW",
            },
            {
                **eligibility["runs"][1],
                "run_accession": "CTRL_NEW",
            },
        ])

        write_tsv(paths["runs"], run_fields, runs)
        write_tsv(
            paths["annotations"],
            annotation_fields,
            annotations,
        )
        write_tsv(
            paths["conditions"],
            condition_fields,
            conditions,
        )
        write_tsv(
            paths["candidates"],
            candidate_fields,
            candidates,
        )

        paths["eligibility"].write_text(
            json.dumps(eligibility),
            encoding="utf-8",
        )

        plan = self.build(paths)

        by_ip = {
            row["ip_run_accession"]: row
            for row in plan
        }

        self.assertEqual(
            by_ip["IP_NEW"]["analysis_status"],
            "ready",
        )
        self.assertEqual(
            by_ip["IP_NEW"]["control_run_accession"],
            "CTRL_NEW",
        )


if __name__ == "__main__":
    unittest.main()
