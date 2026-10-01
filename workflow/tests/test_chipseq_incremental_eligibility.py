#!/usr/bin/env python3

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

SPEC = importlib.util.spec_from_file_location(
    "incremental_eligibility",
    ROOT
    / "workflow/scripts/build_chipseq_incremental_eligibility.py",
)

M = importlib.util.module_from_spec(
    SPEC
)

SPEC.loader.exec_module(
    M
)


def write_tsv(path, fields, rows):
    with Path(path).open(
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


class IncrementalEligibilityTest(
    unittest.TestCase
):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(
            self.tmp.name
        )

        self.paths = {
            name: self.dir / name
            for name in (
                "runs.tsv",
                "annotations.tsv",
                "conditions.tsv",
                "baseline_conditions.tsv",
                "baseline.json",
                "study_evidence.json",
                "policy.json",
            )
        }

        self.run_fields = [
            "run_accession",
            "study_accession",
            "experiment_accession",
            "sample_accession",
            "omics",
            "experiment_title",
            "sample_title",
            "library_name",
            "library_layout",
            "instrument_platform",
        ]

        self.annotation_fields = [
            "run_accession",
            "study_accession",
            "experiment_accession",
            "sample_accession",
            "declared_target",
            "library_role",
            "label_status",
            "target_candidates",
            "control_role_candidates",
        ]

        self.condition_fields = [
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

        self.runs = [
            self.run_row(
                "IP1",
                "STUDY1",
            ),
            self.run_row(
                "CTRL1",
                "STUDY1",
            ),
        ]

        self.annotations = [
            self.annotation_row(
                "IP1",
                "STUDY1",
                "ip",
                "H3K4me3",
            ),
            self.annotation_row(
                "CTRL1",
                "STUDY1",
                "input",
                "",
            ),
        ]

        self.conditions = [
            self.condition_row(
                "IP1",
                "STUDY1",
                "ip",
                "H3K4me3",
            ),
            self.condition_row(
                "CTRL1",
                "STUDY1",
                "input",
                "",
            ),
        ]

        baseline = {
            "schema_version": 1,
            "scope": "experimental_protocol_eligibility_only",
            "runs": [
                self.decision_row(
                    "IP1",
                    "STUDY1",
                    "ip",
                    "H3K4me3",
                ),
                self.decision_row(
                    "CTRL1",
                    "STUDY1",
                    "input",
                    "",
                ),
            ],
        }

        evidence = {
            "schema_version": 1,
            "studies": {
                "STUDY1": {
                    "evidence_type": "publication_methods",
                }
            },
        }

        policy = {
            "schema_version": 1,
            "scope": "incremental_experimental_protocol_eligibility",
            "strategy": "reviewed_protocol_envelope",
            "require_resolved_labels": True,
            "require_resolved_conditions": True,
            "reject_documented_issues_for_new_runs": True,
            "default_new_run_status": "review_required",
            "automatic_studies": {
                "STUDY1": {
                    "automatic_new_runs": True,
                    "required_evidence_type": "publication_methods",
                }
            },
        }

        write_tsv(
            self.paths["runs.tsv"],
            self.run_fields,
            self.runs,
        )

        write_tsv(
            self.paths["annotations.tsv"],
            self.annotation_fields,
            self.annotations,
        )

        write_tsv(
            self.paths["conditions.tsv"],
            self.condition_fields,
            self.conditions,
        )

        write_tsv(
            self.paths["baseline_conditions.tsv"],
            self.condition_fields,
            self.conditions,
        )

        self.paths[
            "baseline.json"
        ].write_text(
            json.dumps(
                baseline
            ),
            encoding="utf-8",
        )

        self.paths[
            "study_evidence.json"
        ].write_text(
            json.dumps(
                evidence
            ),
            encoding="utf-8",
        )

        self.paths[
            "policy.json"
        ].write_text(
            json.dumps(
                policy
            ),
            encoding="utf-8",
        )

    def tearDown(self):
        self.tmp.cleanup()

    def run_row(
        self,
        run,
        study,
    ):
        return {
            "run_accession": run,
            "study_accession": study,
            "experiment_accession": run + "_EXP",
            "sample_accession": run + "_SAMPLE",
            "omics": "ChIP-seq",
            "experiment_title": "",
            "sample_title": "",
            "library_name": "",
            "library_layout": "SINGLE",
            "instrument_platform": "ILLUMINA",
        }

    def annotation_row(
        self,
        run,
        study,
        role,
        target,
        label_status="resolved_from_labels",
    ):
        return {
            "run_accession": run,
            "study_accession": study,
            "experiment_accession": run + "_EXP",
            "sample_accession": run + "_SAMPLE",
            "declared_target": target,
            "library_role": role,
            "label_status": label_status,
            "target_candidates": (
                json.dumps(
                    [target]
                )
                if target
                else "[]"
            ),
            "control_role_candidates": (
                '["input"]'
                if role == "input"
                else "[]"
            ),
        }

    def condition_row(
        self,
        run,
        study,
        role,
        target,
        status="resolved_from_labels",
        documented_issue="false",
        time_value="3",
    ):
        return {
            "run_accession": run,
            "study_accession": study,
            "sample_accession": run + "_SAMPLE",
            "experiment_accession": run + "_EXP",
            "library_role": role,
            "declared_target": target,
            "library_layout": "SINGLE",
            "instrument_platform": "ILLUMINA",
            "condition_status": status,
            "condition_rule": "title_time_label",
            "label_context": "cho-k1",
            "time_value": time_value,
            "time_unit": "days",
            "documented_issue": documented_issue,
            "evidence_json": "[]",
        }

    def decision_row(
        self,
        run,
        study,
        role,
        target,
    ):
        return {
            "run_accession": run,
            "study_accession": study,
            "library_role": role,
            "declared_target": target,
            "experimental_eligibility": "retained_by_rules",
            "rule_id": "baseline_rule",
            "reason": "curated_baseline",
            "evidence_scope": "study_methods_and_run_metadata",
        }

    def append_new(
        self,
        run,
        study="STUDY1",
        role="ip",
        target="H3K4me3",
        label_status="resolved_from_labels",
        condition_status="resolved_from_labels",
        documented_issue="false",
        time_value="3",
    ):
        runs, _ = M.read_tsv(
            self.paths["runs.tsv"]
        )

        annotations, _ = M.read_tsv(
            self.paths["annotations.tsv"]
        )

        conditions, _ = M.read_tsv(
            self.paths["conditions.tsv"]
        )

        runs.append(
            self.run_row(
                run,
                study,
            )
        )

        annotations.append(
            self.annotation_row(
                run,
                study,
                role,
                target,
                label_status,
            )
        )

        conditions.append(
            self.condition_row(
                run,
                study,
                role,
                target,
                condition_status,
                documented_issue,
                time_value,
            )
        )

        write_tsv(
            self.paths["runs.tsv"],
            self.run_fields,
            runs,
        )

        write_tsv(
            self.paths["annotations.tsv"],
            self.annotation_fields,
            annotations,
        )

        write_tsv(
            self.paths["conditions.tsv"],
            self.condition_fields,
            conditions,
        )

    def build(self):
        return M.build(
            self.paths["runs.tsv"],
            self.paths["annotations.tsv"],
            self.paths["conditions.tsv"],
            self.paths["baseline.json"],
            self.paths["baseline_conditions.tsv"],
            self.paths["study_evidence.json"],
            self.paths["policy.json"],
        )

    def by_run(self):
        return {
            row["run_accession"]: row
            for row in self.build()
        }

    def test_baseline_decisions_are_preserved(self):
        decisions = self.by_run()

        self.assertEqual(
            decisions["IP1"][
                "experimental_eligibility"
            ],
            "retained_by_rules",
        )

        self.assertEqual(
            decisions["IP1"][
                "decision_origin"
            ],
            "baseline_curated_review",
        )

    def test_new_matching_ip_is_retained(self):
        self.append_new(
            "IP_NEW"
        )

        decision = self.by_run()[
            "IP_NEW"
        ]

        self.assertEqual(
            decision[
                "experimental_eligibility"
            ],
            "retained_by_rules",
        )

        self.assertEqual(
            decision[
                "decision_origin"
            ],
            "incremental_reviewed_envelope",
        )

    def test_new_matching_input_is_retained(self):
        self.append_new(
            "CTRL_NEW",
            role="input",
            target="",
        )

        decision = self.by_run()[
            "CTRL_NEW"
        ]

        self.assertEqual(
            decision[
                "experimental_eligibility"
            ],
            "retained_by_rules",
        )

    def test_new_study_requires_review(self):
        self.append_new(
            "IP_NEW",
            study="STUDY_NEW",
        )

        decision = self.by_run()[
            "IP_NEW"
        ]

        self.assertEqual(
            decision[
                "experimental_eligibility"
            ],
            "review_required",
        )

        self.assertEqual(
            decision[
                "reason"
            ],
            "new_study_requires_protocol_review",
        )

    def test_new_condition_requires_review(self):
        self.append_new(
            "IP_NEW",
            time_value="7",
        )

        decision = self.by_run()[
            "IP_NEW"
        ]

        self.assertEqual(
            decision[
                "experimental_eligibility"
            ],
            "review_required",
        )

        self.assertEqual(
            decision[
                "reason"
            ],
            "run_outside_reviewed_protocol_envelope",
        )

    def test_unresolved_label_requires_review(self):
        self.append_new(
            "IP_NEW",
            label_status="unresolved",
        )

        decision = self.by_run()[
            "IP_NEW"
        ]

        self.assertEqual(
            decision[
                "reason"
            ],
            "label_annotation_unresolved",
        )

    def test_documented_issue_requires_review(self):
        self.append_new(
            "IP_NEW",
            documented_issue="true",
        )

        decision = self.by_run()[
            "IP_NEW"
        ]

        self.assertEqual(
            decision[
                "reason"
            ],
            "documented_issue_requires_review",
        )

    def test_missing_baseline_run_is_rejected(self):
        runs, _ = M.read_tsv(
            self.paths["runs.tsv"]
        )

        annotations, _ = M.read_tsv(
            self.paths["annotations.tsv"]
        )

        conditions, _ = M.read_tsv(
            self.paths["conditions.tsv"]
        )

        runs = [
            row
            for row in runs
            if row["run_accession"]
            != "IP1"
        ]

        annotations = [
            row
            for row in annotations
            if row["run_accession"]
            != "IP1"
        ]

        conditions = [
            row
            for row in conditions
            if row["run_accession"]
            != "IP1"
        ]

        write_tsv(
            self.paths["runs.tsv"],
            self.run_fields,
            runs,
        )

        write_tsv(
            self.paths["annotations.tsv"],
            self.annotation_fields,
            annotations,
        )

        write_tsv(
            self.paths["conditions.tsv"],
            self.condition_fields,
            conditions,
        )

        with self.assertRaisesRegex(
            ValueError,
            "disappeared",
        ):
            self.build()


if __name__ == "__main__":
    unittest.main()
