#!/usr/bin/env python3
"""Static and operational guards for the dynamic Phase 6B rule graph."""

import csv
import importlib.util
import unittest
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PLAN = ROOT / "config/chipseq_peak_calling_plan.tsv"
RULES = ROOT / "workflow/rules/chipseq_peak_calling.smk"
PREPARE = ROOT / "workflow/scripts/prepare_chipseq_peak_calling_runtime.py"

spec = importlib.util.spec_from_file_location("prepare_peak_runtime", PREPARE)
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


class DynamicRuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with PLAN.open(encoding="utf-8", newline="") as handle:
            cls.rows = list(csv.DictReader(handle, delimiter="\t"))
        cls.rules = RULES.read_text(encoding="utf-8")

    def test_plan_drives_exact_dynamic_counts(self):
        self.assertEqual(len(self.rows), 18)
        self.assertEqual(Counter(row["peak_mode"] for row in self.rows), {
            "narrow": 8, "broad": 10,
        })
        self.assertEqual(Counter(row["fragment_size_policy"] for row in self.rows), {
            "fixed": 6, "phantompeakqualtools": 12,
        })

    def test_prjna_has_zero_fragment_estimation_dependencies(self):
        prjna = [row for row in self.rows if row["study_accession"] == "PRJNA865478"]
        self.assertEqual(len(prjna), 6)
        self.assertTrue(all(row["fragment_size_policy"] == "fixed" for row in prjna))
        self.assertTrue(all(row["fragment_size_bp"] == "147" for row in prjna))
        self.assertTrue(all(row["execution_status"] == "ready_for_peak_calling" for row in prjna))

    def test_prjeb_all_depend_on_own_estimate(self):
        prjeb = [row for row in self.rows if row["study_accession"] == "PRJEB9291"]
        self.assertEqual(len(prjeb), 12)
        self.assertTrue(all(row["fragment_size_policy"] == "phantompeakqualtools" for row in prjeb))
        self.assertTrue(all(row["fragment_size_bp"] == "" for row in prjeb))
        self.assertTrue(all(row["execution_status"] == "requires_fragment_estimation" for row in prjeb))

    def test_rule_graph_separates_narrow_and_broad_contracts(self):
        self.assertIn("rule chipseq_call_narrow_peaks:", self.rules)
        self.assertIn("rule chipseq_call_broad_peaks:", self.rules)
        broad = self.rules.split("rule chipseq_call_broad_peaks:", 1)[1].split(
            "rule chipseq_peak_frip:", 1
        )[0]
        self.assertIn("peaks.gappedPeak", broad)
        self.assertIn("peaks.broadPeak", broad)
        self.assertNotIn("summits", broad)

    def test_missing_real_filtering_pointer_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "remains fail-closed"):
            prepare.completed_filtering_job(ROOT / "definitely_missing_peak_test_pointer.txt")


if __name__ == "__main__":
    unittest.main(verbosity=2)
