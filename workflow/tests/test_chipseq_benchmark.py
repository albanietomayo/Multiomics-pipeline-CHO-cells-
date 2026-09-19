#!/usr/bin/env python3
"""Regression tests for the non-executing ChIP-seq benchmark control layer."""

import csv
import gzip
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
SBATCH_PATH = ROOT / "workflow/slurm/chipseq_benchmark_p10_p50_p90.sbatch"


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
        self.reference_fixture_count = 0

    def tearDown(self):
        self.temporary.cleanup()

    def fixture_root(self):
        root = self.temp / "repository"
        (root / "config").mkdir(parents=True)
        (root / "benchmarks").mkdir()
        shutil.copy2(ROOT / "config/samples.tsv", root / "config/samples.tsv")
        shutil.copy2(SUPPORT.SELECTION, root / "benchmarks/selection.tsv")
        return root

    def fixture_reference(self):
        self.reference_fixture_count += 1
        root = (self.temp / f"shared_reference_{self.reference_fixture_count}").resolve()
        index = root / "bowtie2_index"
        index.mkdir(parents=True)
        combined = b">chr1\nACGT\n>NC_007936.1\nACGT\n"
        with (root / "genome_plus_mt.fa.gz").open("wb") as compressed:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=compressed, mtime=0
            ) as handle:
                handle.write(combined)
        (root / "genome_plus_mt.fa.fai").write_text(
            "chr1\t4\t6\t4\t5\nNC_007936.1\t4\t24\t4\t5\n",
            encoding="utf-8",
        )
        synthetic_content_hashes = {
            "nuclear": "1" * 64,
            "mitochondrial": "2" * 64,
            "mapping": "3" * 64,
        }
        provenance = {
            "annotation_release_independently_verified": False,
            "configured_annotation_release": 104,
            "nuclear_accession": SUPPORT.REFERENCE_IDENTITY["refseq_accession"],
            "mitochondrial_accession": SUPPORT.REFERENCE_IDENTITY["mitochondrial_accession"],
            "nuclear_sha256": synthetic_content_hashes["nuclear"],
            "mitochondrial_sha256": synthetic_content_hashes["mitochondrial"],
            "mapping_sha256": synthetic_content_hashes["mapping"],
            "mitochondrial_length": 4,
        }
        (root / "reference_provenance.json").write_text(
            json.dumps(provenance) + "\n", encoding="utf-8"
        )
        for suffix in ("1", "2", "3", "4", "rev.1", "rev.2"):
            (index / f"genome_plus_mt.{suffix}.bt2l").write_text(
                f"index-{suffix}\n", encoding="utf-8"
            )
        source_hashes = {
            relative: SUPPORT.sha256(root / relative)
            for relative in SUPPORT.REFERENCE_SOURCE_SHA256
        }
        patches = [
            mock.patch.dict(SUPPORT.REFERENCE_SOURCE_SHA256, source_hashes, clear=True),
            mock.patch.dict(SUPPORT.REFERENCE_CONTENT_SHA256, synthetic_content_hashes, clear=True),
            mock.patch.object(SUPPORT, "EXPECTED_NUCLEAR_SPAN_BP", 4),
            mock.patch.object(SUPPORT, "EXPECTED_MITOCHONDRIAL_LENGTH_BP", 4),
        ]
        for patcher in patches:
            patcher.start()
            self.addCleanup(patcher.stop)
        SUPPORT.write_reference_inventory(root)
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
        argv = [str(SUBMIT_PATH), "--benchmark-class", "P10", "--shared-reference-root", "/fixture"]
        with mock.patch.object(sys, "argv", argv), \
                mock.patch.object(SUBMIT, "local_validation", return_value=plan), \
                mock.patch.object(SUBMIT, "prepare_submission") as prepare, \
                mock.patch.object(SUBMIT.subprocess, "check_output") as check_output, \
                mock.patch("builtins.print"):
            SUBMIT.main()
        prepare.assert_not_called()
        check_output.assert_not_called()

    def test_submit_requires_clean_source_and_forbids_dirty_override(self):
        argv = [
            str(SUBMIT_PATH), "--benchmark-class", "P10", "--shared-reference-root",
            "/fixture", "--submit", "--development-dirty-check",
        ]
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


    def test_shared_reference_inventory_is_complete_and_fai_derived(self):
        root = self.fixture_reference()
        inventory = SUPPORT.reference_inventory(root)
        self.assertEqual(inventory["identity"]["refseq_accession"], "GCF_003668045.3")
        self.assertEqual(inventory["fai"]["nuclear_span_bp"], 4)
        self.assertEqual(len(inventory["required_index_components"]), 6)
        self.assertEqual(
            {item["relative_path"] for item in inventory["files"]},
            set(SUPPORT.REFERENCE_FILES),
        )

    def test_annotations_and_uncompressed_fasta_are_not_operational_inputs(self):
        root = self.fixture_reference()
        inventory = SUPPORT.reference_inventory(root)
        paths = {item["relative_path"] for item in inventory["files"]}
        self.assertNotIn("annotation.gff3", paths)
        self.assertNotIn("annotation.gtf", paths)
        self.assertNotIn("genome_plus_mt.fa", paths)
        self.assertIn("genome_plus_mt.fa.gz", paths)

    def test_shared_reference_source_hash_and_inventory_fail_closed(self):
        root = self.fixture_reference()
        (root / "genome_plus_mt.fa.gz").write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "source hash mismatch"):
            SUPPORT.reference_inventory(root)

        root = self.fixture_reference()
        inventory_path = root / SUPPORT.REFERENCE_INVENTORY
        inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
        inventory["build"]["version"] = "different"
        inventory_path.write_text(json.dumps(inventory) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "immutable reference inventory"):
            SUPPORT.reference_inventory(root)

    def test_shared_reference_validation_fails_closed(self):
        root = self.fixture_reference()
        (root / "bowtie2_index/genome_plus_mt.4.bt2l").unlink()
        with self.assertRaisesRegex(ValueError, "Missing or empty"):
            SUPPORT.reference_inventory(root)
        with self.assertRaisesRegex(ValueError, "explicit absolute"):
            SUPPORT.reference_inventory(Path("relative/reference"))

    def test_shared_reference_change_after_submission_fails_closed(self):
        root = self.fixture_reference()
        expected = self.temp / "submission_plan.json"
        expected.write_text(
            json.dumps({"shared_reference": SUPPORT.reference_inventory(root)}),
            encoding="utf-8",
        )
        (root / "bowtie2_index/genome_plus_mt.1.bt2l").write_text("changed\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "immutable (reference|submission) inventory"):
            SUPPORT.validate_reference_inventory(root, expected)

    def test_stage_sizes_and_hashes_recorded_before_transient_disappears(self):
        transient = self.temp / "raw.bam"
        marker = self.temp / "stage.json"
        transient.write_bytes(b"bam-bytes")
        SUPPORT.write_artifact_record("alignment", marker, [transient], [])
        transient.unlink()
        record = json.loads(marker.read_text(encoding="utf-8"))
        self.assertEqual(record["artifacts"][0]["bytes"], 9)
        self.assertEqual(record["artifacts"][0]["retention"], "transient")
        self.assertTrue(record["artifacts"][0]["regenerable"])

    def test_persistence_copy_is_verified_and_collision_protected(self):
        source = self.temp / "filtered.bam"
        source.write_bytes(b"filtered")
        destination = self.temp / "persistent"
        outputs = SUPPORT.persist_files(destination, [f"filtered/filtered.bam={source}"])
        self.assertEqual(outputs[0].read_bytes(), b"filtered")
        with self.assertRaisesRegex(ValueError, "collision"):
            SUPPORT.persist_files(destination, [f"filtered/filtered.bam={source}"])

    def test_transient_raw_fastq_policy(self):
        script = SBATCH_PATH.read_text(encoding="utf-8")
        self.assertNotIn("data/raw/$CHIP_FILENAME=", script)
        self.assertIn("run_stage acquisition_validation", script)

    def test_transient_processed_fastq_and_no_alignment_copy_policy(self):
        script = SBATCH_PATH.read_text(encoding="utf-8")
        self.assertNotIn("preprocessing/$CHIP_RUN.fastq.gz=", script)
        self.assertNotIn("cp @PROCESSED@", script)
        self.assertIn("-U @PROCESSED@", script)

    def test_transient_raw_and_duplicate_marked_bam_policy(self):
        script = SBATCH_PATH.read_text(encoding="utf-8")
        persistence = script.split('CHIP_STAGE="successful_persistence"', 1)[1]
        self.assertNotIn("raw.sorted.bam=", persistence)
        self.assertNotIn("dupmarked.sorted.bam=", persistence)
        self.assertIn("dupmarked.sorted.bam", script)

    def test_retained_filtered_bam_qc_metrics_and_provenance_policy(self):
        script = SBATCH_PATH.read_text(encoding="utf-8")
        for required in (
            "filtered/filtered.bam=", "filtered/filtered.bam.csi=",
            "qc/fastp/", "qc/fastqc/", "picard_metrics.txt",
            "benchmark_metrics.tsv", "reference_provenance.json",
            "persistent_footprint.tsv", "pipeline_commit.txt",
            "environment_explicit.txt",
        ):
            self.assertIn(required, script)

    def test_no_persistent_per_run_reference_or_index_duplication(self):
        script = SBATCH_PATH.read_text(encoding="utf-8")
        self.assertIn("--shared-reference-root", script)
        self.assertNotIn("bowtie2-build", script)
        self.assertNotIn('CHIP_FASTA=', script)
        self.assertIn("genome_plus_mt.fa.fai.part", script)
        self.assertNotIn("persist_work", script)

    def test_completion_marker_follows_successful_persistence(self):
        script = SBATCH_PATH.read_text(encoding="utf-8")
        self.assertLess(script.index('"$CHIP_SUPPORT" persist-files'), script.index(
            '"$CHIP_SUPPORT" complete-stage --stage benchmark'
        ))
        self.assertEqual(script.count("complete-stage --stage benchmark"), 1)
        self.assertLess(
            script.index('> "$CHIP_SAVE/persistent_footprint.tsv"'),
            script.index('mv "$CHIP_COMPLETION_PENDING" "$CHIP_SAVE/benchmark.complete.json"'),
        )
        self.assertLess(
            script.index('stage\\tcompleted'),
            script.index('CHIP_PERSISTENT_BASE_BYTES='),
        )
        self.assertLess(
            script.index('mv "$CHIP_COMPLETION_PENDING" "$CHIP_SAVE/benchmark.complete.json"'),
            script.rindex("trap - EXIT"),
        )
        self.assertEqual(script.count("trap - EXIT"), 2)
        self.assertIn('--marker "$CHIP_COMPLETION_PENDING"', script)

    def test_worker_records_required_benchmark_resource_metrics(self):
        script = SBATCH_PATH.read_text(encoding="utf-8")
        for field in (
            "elapsed_seconds", "max_rss_kbytes", "user_cpu_seconds",
            "system_cpu_seconds", "baseline_scratch_bytes", "peak_scratch_bytes",
            "incremental_peak_scratch_bytes", "input_bytes", "output_bytes",
        ):
            self.assertIn(field, script)


if __name__ == "__main__":
    unittest.main(verbosity=2)
