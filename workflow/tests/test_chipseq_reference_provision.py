#!/usr/bin/env python3
"""Synthetic tests for one-time shared-reference provisioning."""

import gzip
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
PROVISION_PATH = ROOT / "workflow/scripts/chipseq_reference_provision.py"
SUBMIT_PATH = ROOT / "workflow/slurm/submit_chipseq_reference_provision.py"
SBATCH_PATH = ROOT / "workflow/slurm/chipseq_provision_shared_reference.sbatch"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PROVISION = load("chipseq_reference_provision_test", PROVISION_PATH)
SUBMIT = load("submit_chipseq_reference_provision_test", SUBMIT_PATH)


class ChipseqReferenceProvisionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.temp = Path(self.temporary.name)
        self.source = self.temp / "historical"
        self.source.mkdir()
        self.combined = b">chr1\nACGT\n>NC_007936.1\nACGT\n"
        with (self.source / "genome_plus_mt.fa.gz").open("wb") as compressed:
            with gzip.GzipFile(filename="", mode="wb", fileobj=compressed, mtime=0) as handle:
                handle.write(self.combined)
        (self.source / "genome_plus_mt.fa.fai").write_text(
            "chr1\t4\t6\t4\t5\nNC_007936.1\t4\t24\t4\t5\n",
            encoding="utf-8",
        )
        mapping_hash = PROVISION.hashlib.sha256(self.combined).hexdigest()
        provenance = {
            "annotation_release_independently_verified": False,
            "configured_annotation_release": 104,
            "nuclear_accession": PROVISION.SUPPORT.REFERENCE_IDENTITY["refseq_accession"],
            "mitochondrial_accession": "NC_007936.1",
            "nuclear_sha256": "1" * 64,
            "mitochondrial_sha256": "2" * 64,
            "mapping_sha256": mapping_hash,
            "mitochondrial_length": 4,
        }
        (self.source / "reference_provenance.json").write_text(
            json.dumps(provenance, sort_keys=True) + "\n", encoding="utf-8"
        )
        self.mapping_hash = mapping_hash
        self.patches = [
            mock.patch.object(PROVISION, "EXPECTED_MAPPING_SHA256", mapping_hash),
            mock.patch.object(PROVISION, "EXPECTED_SEQUENCE_COUNT", 2),
            mock.patch.object(PROVISION, "EXPECTED_TOTAL_BP", 8),
            mock.patch.object(PROVISION.SUPPORT, "EXPECTED_NUCLEAR_SPAN_BP", 4),
            mock.patch.object(PROVISION.SUPPORT, "EXPECTED_MITOCHONDRIAL_LENGTH_BP", 4),
            mock.patch.dict(
                PROVISION.SUPPORT.REFERENCE_CONTENT_SHA256,
                {"nuclear": "1" * 64, "mitochondrial": "2" * 64, "mapping": mapping_hash},
                clear=True,
            ),
        ]
        for patcher in self.patches:
            patcher.start()
            self.addCleanup(patcher.stop)
        self.refresh_source_contract()

    def tearDown(self):
        self.temporary.cleanup()

    def refresh_source_contract(self):
        sizes = {
            name: (self.source / name).stat().st_size
            for name in PROVISION.SUPPORT.REFERENCE_SOURCE_SHA256
        }
        hashes = {
            name: PROVISION.SUPPORT.sha256(self.source / name)
            for name in PROVISION.SUPPORT.REFERENCE_SOURCE_SHA256
        }
        if hasattr(self, "size_patch"):
            self.size_patch.stop()
            self.hash_patch.stop()
        self.size_patch = mock.patch.dict(
            PROVISION.SUPPORT.REFERENCE_SOURCE_BYTES, sizes, clear=True
        )
        self.hash_patch = mock.patch.dict(
            PROVISION.SUPPORT.REFERENCE_SOURCE_SHA256, hashes, clear=True
        )
        self.size_patch.start()
        self.hash_patch.start()
        self.addCleanup(self.size_patch.stop)
        self.addCleanup(self.hash_patch.stop)

    def reference_tree(self, name="reference"):
        root = self.temp / name
        (root / "bowtie2_index").mkdir(parents=True)
        for source in PROVISION.SUPPORT.REFERENCE_SOURCE_SHA256:
            os.link(self.source / source, root / source)
        for relative in PROVISION.SUPPORT.INDEX_COMPONENTS:
            (root / relative).write_bytes(("synthetic-" + relative).encode())
        PROVISION.SUPPORT.write_reference_inventory(root)
        return root

    def test_source_sizes_hashes_decompressed_hash_and_fai_identity(self):
        report = PROVISION.validate_source(self.source)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["decompressed_sha256"], self.mapping_hash)
        self.assertEqual(report["fai"]["sequence_count"], 2)
        self.assertEqual(report["fai"]["total_bp"], 8)
        self.assertEqual(len(report["files"]), 3)

    def test_incorrect_source_size_rejected(self):
        expected = dict(PROVISION.SUPPORT.REFERENCE_SOURCE_BYTES)
        expected["genome_plus_mt.fa.fai"] += 1
        with mock.patch.dict(PROVISION.SUPPORT.REFERENCE_SOURCE_BYTES, expected, clear=True):
            with self.assertRaisesRegex(ValueError, "size mismatch"):
                PROVISION.validate_source(self.source, verify_decompressed=False)

    def test_incorrect_source_hash_rejected(self):
        path = self.source / "genome_plus_mt.fa.fai"
        data = bytearray(path.read_bytes())
        data[-2] = ord("6")
        path.write_bytes(data)
        with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
            PROVISION.validate_source(self.source, verify_decompressed=False)

    def test_incorrect_decompressed_hash_rejected(self):
        with (self.source / "genome_plus_mt.fa.gz").open("wb") as compressed:
            with gzip.GzipFile(filename="", mode="wb", fileobj=compressed, mtime=0) as handle:
                handle.write(b"different")
        self.refresh_source_contract()
        with self.assertRaisesRegex(ValueError, "Decompressed combined FASTA"):
            PROVISION.validate_source(self.source)

    def test_fai_reference_identity_rejected(self):
        fai = self.source / "genome_plus_mt.fa.fai"
        fai.write_text("chr1\t4\t6\t4\t5\nNC_007936.1\t5\t24\t5\t6\n", encoding="utf-8")
        self.refresh_source_contract()
        with self.assertRaisesRegex(ValueError, "FAI identity"):
            PROVISION.validate_source(self.source, verify_decompressed=False)

    def test_all_six_index_components_are_required(self):
        root = self.reference_tree()
        observed = PROVISION.SUPPORT.reference_inventory(root)
        self.assertEqual(len(observed["required_index_components"]), 6)

    def test_incomplete_index_rejected(self):
        root = self.reference_tree()
        (root / PROVISION.SUPPORT.INDEX_COMPONENTS[-1]).unlink()
        with self.assertRaisesRegex(ValueError, "Missing or empty"):
            PROVISION.SUPPORT.reference_inventory(root)

    def test_complete_destination_collision_rejected(self):
        root = self.reference_tree()
        with self.assertRaisesRegex(ValueError, "Complete destination"):
            PROVISION.validate_destination(self.source, root)

    def test_partial_destination_collision_rejected(self):
        root = self.temp / "partial"
        root.mkdir()
        with self.assertRaisesRegex(ValueError, "Partial or invalid"):
            PROVISION.validate_destination(self.source, root)

    def test_stale_hidden_staging_collision_rejected(self):
        destination = self.temp / "parent/reference"
        stale = destination.parent / ".chipseq.provisioning.123.staging"
        stale.mkdir(parents=True)
        with self.assertRaisesRegex(ValueError, "Partial provisioning staging"):
            PROVISION.validate_destination(self.source, destination)

    def test_inventory_is_deterministic(self):
        first = self.reference_tree("first")
        second = self.reference_tree("second")
        self.assertEqual(
            (first / "reference_inventory.json").read_bytes(),
            (second / "reference_inventory.json").read_bytes(),
        )

    def test_no_persistent_uncompressed_fasta(self):
        script = SBATCH_PATH.read_text(encoding="utf-8")
        self.assertIn('CHIP_FASTA="$CHIP_WORK/genome_plus_mt.fa"', script)
        self.assertNotIn('"$CHIP_STAGE/genome_plus_mt.fa"', script)
        self.assertIn('[[ ! -e "$CHIP_DEST/genome_plus_mt.fa" ]]', script)

    def test_quota_calculation_and_p10_reserve(self):
        report = PROVISION.storage_report(10_000, 1_000, 2_000, 100, 1_150)
        self.assertEqual(report["current_free_bytes"], 9_000)
        self.assertEqual(report["additional_persistent_bytes_required"], 2_100)
        self.assertEqual(report["expected_free_bytes_after_publication"], 6_900)
        self.assertEqual(report["expected_free_bytes_after_publication_and_p10_reserve"], 5_750)

    def test_insufficient_quota_blocks_publication(self):
        with self.assertRaisesRegex(ValueError, "Insufficient quota"):
            PROVISION.storage_report(10_000, 8_000, 1_000, 1, 1_000)

    def test_p10_reserve_is_exact_committed_value(self):
        self.assertEqual(PROVISION.P10_RESERVE_BYTES, 1_234_803_098)

    def test_check_only_cannot_call_sbatch(self):
        report = {
            "repository": {"head": "a" * 40},
            "source": {},
            "destination": {},
            "quota": {},
        }
        with mock.patch.object(sys, "argv", [str(SUBMIT_PATH)]), \
                mock.patch.object(SUBMIT.PROVISION, "preflight", return_value=report), \
                mock.patch.object(SUBMIT.subprocess, "check_output") as check_output, \
                mock.patch("builtins.print"):
            SUBMIT.main()
        check_output.assert_not_called()

    def test_check_only_reports_external_scheduler_log_paths(self):
        report = {
            "repository": {"head": "a" * 40},
            "source": {},
            "destination": {},
            "quota": {},
        }
        log_root = self.temp / "check-only-logs"
        stdout_path = log_root / "job-%j.out"
        stderr_path = log_root / "job-%j.err"
        with mock.patch.object(sys, "argv", [str(SUBMIT_PATH)]), \
                mock.patch.object(SUBMIT.PROVISION, "preflight", return_value=report), \
                mock.patch.object(SUBMIT, "LOG_ROOT", log_root), \
                mock.patch.object(SUBMIT, "STDOUT_PATH", stdout_path), \
                mock.patch.object(SUBMIT, "STDERR_PATH", stderr_path), \
                mock.patch.object(SUBMIT.subprocess, "check_output") as check_output, \
                mock.patch("builtins.print") as output:
            SUBMIT.main()
        rendered = json.loads(output.call_args_list[0].args[0])
        logs = rendered["scheduler_logs"]
        self.assertEqual(logs["root"], str(log_root))
        self.assertEqual(logs["stdout"], str(stdout_path))
        self.assertEqual(logs["stderr"], str(stderr_path))
        self.assertTrue(logs["outside_repository"])
        self.assertFalse(log_root.exists())
        self.assertFalse(stdout_path.exists())
        self.assertFalse(stderr_path.exists())
        check_output.assert_not_called()

    def test_real_submit_uses_explicit_external_stdout_and_stderr(self):
        report = {
            "repository": {"head": "a" * 40},
            "source": {},
            "destination": {},
            "quota": {},
        }
        log_root = self.temp / "external-logs"
        self.assertFalse(log_root.exists())
        with mock.patch.object(sys, "argv", [str(SUBMIT_PATH), "--submit"]), \
                mock.patch.object(SUBMIT.PROVISION, "preflight", return_value=report), \
                mock.patch.object(SUBMIT, "LOG_ROOT", log_root), \
                mock.patch.object(SUBMIT, "STDOUT_PATH", log_root / "job-%j.out"), \
                mock.patch.object(SUBMIT, "STDERR_PATH", log_root / "job-%j.err"), \
                mock.patch.object(
                    SUBMIT.subprocess, "check_output", return_value="12345\n"
                ) as check_output, \
                mock.patch("builtins.print"):
            SUBMIT.main()
        self.assertTrue(log_root.is_dir())
        command = check_output.call_args.args[0]
        self.assertIn(f"--output={log_root / 'job-%j.out'}", command)
        self.assertIn(f"--error={log_root / 'job-%j.err'}", command)
        self.assertNotIn("slurm-%j.out", " ".join(command))
        self.assertEqual(check_output.call_count, 1)

    def test_scheduler_log_root_must_be_outside_repository(self):
        with mock.patch.object(SUBMIT, "LOG_ROOT", ROOT / "scheduler-logs"):
            with self.assertRaisesRegex(ValueError, "outside the Git repository"):
                SUBMIT.scheduler_log_report()

    def test_submitter_has_one_explicit_sbatch_call_site(self):
        source = SUBMIT_PATH.read_text(encoding="utf-8")
        self.assertEqual(source.count('"sbatch"'), 1)
        self.assertIn('f"--output={STDOUT_PATH}"', source)
        self.assertIn('f"--error={STDERR_PATH}"', source)

    def test_clean_worktree_guard_remains_fail_closed(self):
        values = {
            ("branch", "--show-current"): PROVISION.EXPECTED_BRANCH,
            ("rev-parse", "HEAD"): "a" * 40,
            ("status", "--short"): "?? scheduler-created.out",
        }
        with mock.patch.object(
            PROVISION, "git", side_effect=lambda *args: values[args]
        ), mock.patch.object(PROVISION.subprocess, "run"):
            with self.assertRaisesRegex(ValueError, "worktree is dirty"):
                PROVISION.repository_preflight()

    def test_failed_job_10319325_destination_state_is_safe_to_retry(self):
        destination = self.temp / "canonical" / "chipseq"
        report = PROVISION.validate_destination(self.source, destination)
        self.assertFalse(report["collision"])
        self.assertFalse(destination.exists())
        self.assertEqual(
            list(destination.parent.glob(".chipseq.provisioning.*.staging")), []
        )

    def test_atomic_publication_and_completion_semantics(self):
        script = SBATCH_PATH.read_text(encoding="utf-8")
        validate_stage = 'check-reference --shared-reference-root "$CHIP_STAGE"'
        publish = 'mv -T -- "$CHIP_STAGE" "$CHIP_DEST"'
        complete = 'echo "REFERENCE_PROVISIONING=COMPLETE"'
        self.assertLess(script.index(validate_stage), script.index(publish))
        self.assertLess(script.index(publish), script.index(complete))
        self.assertEqual(script.count(complete), 1)
        self.assertIn("rm -rf -- \"$CHIP_STAGE\"", script)

    def test_hardlink_reuse_and_generated_component_hashing(self):
        script = SBATCH_PATH.read_text(encoding="utf-8")
        self.assertIn('ln "$CHIP_SOURCE/$name" "$CHIP_STAGE/$name"', script)
        self.assertNotIn('cp -- "$CHIP_SOURCE/', script)
        root = self.reference_tree()
        inventory = PROVISION.SUPPORT.reference_inventory(root)
        self.assertTrue(all(item["bytes"] > 0 and len(item["sha256"]) == 64 for item in inventory["files"]))
        with mock.patch.object(PROVISION, "current_used_bytes", return_value=1_000):
            report = PROVISION.component_report(root, self.temp)
        self.assertGreater(report["inventory_file"]["bytes"], 0)
        self.assertEqual(len(report["inventory_file"]["sha256"]), 64)


if __name__ == "__main__":
    unittest.main(verbosity=2)
