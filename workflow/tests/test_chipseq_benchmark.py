#!/usr/bin/env python3
"""Regression tests for the non-executing ChIP-seq benchmark control layer."""

import csv
import importlib.util
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
SUPPORT_PATH = ROOT / "workflow/scripts/chipseq_benchmark.py"
SUBMIT_PATH = ROOT / "workflow/slurm/submit_chipseq_benchmark.py"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SUPPORT = load("chipseq_benchmark_support_test", SUPPORT_PATH)
SUBMIT = load("chipseq_benchmark_submit_test", SUBMIT_PATH)


class ChipseqBenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.temp = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def fixture_root(self):
        root = self.temp / "repository"
        (root / "config").mkdir(parents=True)
        (root / "benchmarks").mkdir()
        shutil.copy2(ROOT / "config/samples.tsv", root / "config/samples.tsv")
        shutil.copy2(SUPPORT.SELECTION, root / "benchmarks/selection.tsv")
        return root

    def test_exact_manifest_membership_count_and_uniqueness(self):
        rows = SUPPORT.validate_manifest()
        self.assertEqual(len(rows), 3)
        self.assertEqual(
            {row["benchmark_class"]: row["run_accession"] for row in rows},
            SUPPORT.EXPECTED,
        )
        self.assertEqual(len({row["run_accession"] for row in rows}), 3)

    def test_exact_bytes_md5_layout_and_urls(self):
        expected = {
            "ERR868176": ("1273744805", "099ff78a49450e5e9fa92c6557a65b73"),
            "ERR868152": ("1678499646", "5e0becc282f6a3b34e5659842f2695af"),
            "SRR20770294": ("2193948334", "b8930b0ee55c390b54bebd21ae46ce66"),
        }
        for row in SUPPORT.validate_manifest():
            self.assertEqual((row["fastq_bytes"], row["fastq_md5"]), expected[row["run_accession"]])
            self.assertEqual(row["library_layout"], "SINGLE")
            self.assertEqual(row["instrument_platform"], "ILLUMINA")
            self.assertTrue(row["fastq_ftp"].startswith("ftp.sra.ebi.ac.uk/"))
            self.assertTrue(row["fastq_ftp"].endswith("/" + row["filename"]))

    def test_manifest_is_deterministic_from_authoritative_metadata(self):
        expected = SUPPORT.tsv_text(SUPPORT.MANIFEST_FIELDS, SUPPORT.generate_manifest_rows())
        self.assertEqual(SUPPORT.MANIFEST.read_text(encoding="utf-8"), expected)

    def test_missing_selected_run_fails_closed(self):
        root = self.fixture_root()
        source = root / "config/samples.tsv"
        lines = source.read_text(encoding="utf-8").splitlines()
        source.write_text("\n".join(line for line in lines if not line.startswith("ERR868176\t")) + "\n")
        with mock.patch.object(SUPPORT, "git_value", return_value="a" * 40):
            with self.assertRaisesRegex(ValueError, "found 0"):
                SUPPORT.generate_manifest_rows(source, root / "benchmarks/selection.tsv", root)

    def test_duplicate_selected_run_fails_closed(self):
        root = self.fixture_root()
        selection = root / "benchmarks/selection.tsv"
        selection.write_text(
            "benchmark_class\trun_accession\nP10\tERR868176\nP50\tERR868176\nP90\tSRR20770294\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "must be exactly"):
            SUPPORT.load_selection(selection)

    def test_duplicate_authoritative_run_fails_closed(self):
        root = self.fixture_root()
        source = root / "config/samples.tsv"
        lines = source.read_text(encoding="utf-8").splitlines()
        duplicate = next(line for line in lines if line.startswith("ERR868152\t"))
        source.write_text("\n".join(lines + [duplicate]) + "\n", encoding="utf-8")
        with mock.patch.object(SUPPORT, "git_value", return_value="a" * 40):
            with self.assertRaisesRegex(ValueError, "found 2"):
                SUPPORT.generate_manifest_rows(source, root / "benchmarks/selection.tsv", root)

    def test_metadata_mismatch_fails_manifest_validation(self):
        root = self.fixture_root()
        manifest = root / "benchmarks/manifest.tsv"
        with mock.patch.object(SUPPORT, "git_value", return_value="a" * 40):
            SUPPORT.write_manifest(manifest, root / "config/samples.tsv", root / "benchmarks/selection.tsv", root)
            text = (root / "config/samples.tsv").read_text(encoding="utf-8")
            (root / "config/samples.tsv").write_text(text.replace("1273744805", "1273744806", 1), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "differs"):
                SUPPORT.validate_manifest(manifest, root / "config/samples.tsv", root / "benchmarks/selection.tsv", root)

    def test_unsafe_accession_and_class_fail_closed(self):
        rows = SUPPORT.validate_manifest()
        with self.assertRaisesRegex(ValueError, "Unsafe benchmark class"):
            SUPPORT.manifest_item(rows, "P99")
        tampered = [dict(row) for row in rows]
        tampered[0]["run_accession"] = "SRR1"
        with self.assertRaisesRegex(ValueError, "mismatched"):
            SUPPORT.manifest_item(tampered, "P10")

    def test_storage_preflight_reports_known_and_unknown_components(self):
        row = SUPPORT.manifest_item(SUPPORT.validate_manifest(), "P10")
        usage = shutil._ntuple_diskusage(total=10_000_000_000, used=1, free=9_999_999_999)
        with mock.patch.object(SUPPORT.shutil, "disk_usage", return_value=usage):
            report = SUPPORT.storage_preflight(row, self.temp / "output", self.temp / "scratch")
        self.assertEqual(report["known_compressed_input_bytes"], 1273744805)
        self.assertEqual(report["downstream_expansion_bytes"], "unknown_pending_empirical_benchmark")
        self.assertTrue(report["minimum_known_storage_satisfied"])

    def test_storage_minimum_and_output_collision_fail_closed(self):
        row = SUPPORT.manifest_item(SUPPORT.validate_manifest(), "P90")
        usage = shutil._ntuple_diskusage(total=100, used=1, free=99)
        with mock.patch.object(SUPPORT.shutil, "disk_usage", return_value=usage):
            with self.assertRaisesRegex(ValueError, "known compressed input"):
                SUPPORT.storage_preflight(row, self.temp / "output", self.temp / "scratch")
        collision = self.temp / "output" / "P90_SRR20770294"
        collision.mkdir(parents=True)
        with self.assertRaisesRegex(ValueError, "collision"):
            SUPPORT.storage_preflight(row, self.temp / "output", self.temp / "scratch")

    def test_check_mode_has_no_submission_side_effect(self):
        plan = {"mode": "check"}
        argv = [str(SUBMIT_PATH), "--benchmark-class", "P10"]
        with mock.patch.object(sys, "argv", argv), \
                mock.patch.object(SUBMIT, "local_validation", return_value=plan), \
                mock.patch.object(SUBMIT, "prepare_submission") as prepare, \
                mock.patch.object(SUBMIT.subprocess, "check_output") as check_output, \
                mock.patch("builtins.print"):
            SUBMIT.main()
        prepare.assert_not_called()
        check_output.assert_not_called()

    def test_submit_requires_clean_source_and_forbids_dirty_override(self):
        argv = [str(SUBMIT_PATH), "--benchmark-class", "P10", "--submit", "--development-dirty-check"]
        with mock.patch.object(sys, "argv", argv), self.assertRaisesRegex(ValueError, "forbidden"):
            SUBMIT.main()

    def test_sacct_parser_handles_representative_rows(self):
        raw = self.temp / "sacct.psv"
        output = self.temp / "sacct.tsv"
        raw.write_text(
            "JobID|JobName|State|Elapsed|TotalCPU|AllocCPUS|MaxRSS|AveRSS|ExitCode\n"
            "123|chipseq|COMPLETED|00:01:02|00:03:00|8|123456K|100000K|0:0\n",
            encoding="utf-8",
        )
        SUPPORT.parse_sacct(raw, output)
        rows = SUPPORT.read_tsv(output)
        self.assertEqual(rows[0]["MaxRSS"], "123456K")
        self.assertEqual(rows[0]["ExitCode"], "0:0")

    def test_sacct_missing_values_are_explicit(self):
        raw = self.temp / "sacct.psv"
        output = self.temp / "sacct.tsv"
        raw.write_text(
            "JobID|JobName|State|Elapsed|TotalCPU|AllocCPUS|MaxRSS|AveRSS|ExitCode\n"
            "123.batch|batch|RUNNING|00:00:02||8|||\n",
            encoding="utf-8",
        )
        SUPPORT.parse_sacct(raw, output)
        row = SUPPORT.read_tsv(output)[0]
        self.assertEqual(row["TotalCPU"], "not_available")
        self.assertEqual(row["MaxRSS"], "not_available")
        self.assertEqual(row["ExitCode"], "not_available")

    def test_sacct_missing_column_fails_clearly(self):
        raw = self.temp / "sacct.psv"
        raw.write_text("JobID|State\n123|COMPLETED\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "lacks requested columns"):
            SUPPORT.parse_sacct(raw, self.temp / "out.tsv")

    def test_pending_rows_are_not_completed_metrics(self):
        summary = SUPPORT.read_tsv(SUPPORT.BENCHMARK_DIR / "benchmark_summary.tsv")
        plan = SUPPORT.read_tsv(SUPPORT.BENCHMARK_DIR / "benchmark_execution_plan.tsv")
        metrics = (SUPPORT.BENCHMARK_DIR / "benchmark_metrics.tsv").read_text(encoding="utf-8").splitlines()
        self.assertTrue(all(row["status"] == "pending_not_run" for row in summary + plan))
        self.assertEqual(len(metrics), 1)

    def test_completion_marker_detects_changed_output(self):
        output = self.temp / "artifact.txt"
        marker = self.temp / "complete.json"
        output.write_text("valid\n", encoding="utf-8")
        SUPPORT.write_completion("synthetic", marker, [output])
        SUPPORT.verify_completion("synthetic", marker)
        output.write_text("changed\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "changed"):
            SUPPORT.verify_completion("synthetic", marker)


if __name__ == "__main__":
    unittest.main(verbosity=2)
