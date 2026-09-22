#!/usr/bin/env python3
import importlib.util
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
SUPPORT_PATH = ROOT / "workflow/scripts/chipseq_production.py"
SBATCH = ROOT / "workflow/slurm/chipseq_production.sbatch"
SUBMIT = ROOT / "workflow/slurm/submit_chipseq_production.py"

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module

SUPPORT = load("chipseq_production_test_support", SUPPORT_PATH)

class ChipseqProductionTests(unittest.TestCase):
    def test_first_pair_is_exact_authoritative_fixed_policy(self):
        selected = SUPPORT.select_analysis("SRR20770297", "SRR20770287")
        row = selected["analysis"]
        self.assertEqual((row["study_accession"], row["declared_target"], row["peak_mode"]),
                         ("PRJNA865478", "H3K4me3", "narrow"))
        self.assertEqual((row["fragment_size_policy"], row["fragment_size_bp"]), ("fixed", "147"))
        self.assertEqual(selected["samples"]["ip"]["fastq_bytes"], 1909331888)
        self.assertEqual(selected["samples"]["input"]["fastq_bytes"], 1992978505)

    def test_wrong_control_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "authoritative pairing"):
            SUPPORT.select_analysis("SRR20770297", "SRR20770295")

    def test_prjeb9291_never_receives_fixed_fallback(self):
        selected = SUPPORT.select_analysis("ERR868151", "ERR868150")
        row = selected["analysis"]
        self.assertEqual(row["fragment_size_policy"], "phantompeakqualtools")
        self.assertEqual(row["fragment_size_bp"], "")
        with self.assertRaisesRegex(ValueError, "no fallback"):
            SUPPORT.DYNAMIC.resolve_fragment_size(row)

    def test_preflight_reserve_fits_fourteen_gib_for_first_pair(self):
        with tempfile.TemporaryDirectory() as temporary:
            scratch = Path(temporary)
            usage = shutil._ntuple_diskusage(total=20 * 1024**3, used=6 * 1024**3, free=14 * 1024**3)
            with mock.patch.object(SUPPORT.shutil, "disk_usage", return_value=usage), \
                 mock.patch.object(SUPPORT.BENCHMARK, "reference_inventory", return_value={"root": "/reference"}):
                report = SUPPORT.preflight("SRR20770297", "SRR20770287", Path("/reference"),
                                           scratch, scratch / "production")
            self.assertLess(report["scratch_structural_reserve_bytes"], 14 * 1024**3)
            self.assertGreaterEqual(report["scratch_structural_reserve_bytes"], 8 * 1024**3)

    def test_worker_encodes_validated_delete_order_and_no_reference_copy(self):
        text = SBATCH.read_text(encoding="utf-8")
        self.assertIn('samtools quickcheck -v "$align/raw.sorted.bam"', text)
        self.assertLess(text.index('samtools quickcheck -v "$root/dupmarked.sorted.bam"'),
                        text.index('rm -- "$align/raw.sorted.bam"'))
        self.assertLess(text.index('chipseq_filtering_support.py verify'),
                        text.index('rm -- "$root/dupmarked.sorted.bam"'))
        self.assertIn('INDEX="$REF/bowtie2_index/genome_plus_mt"', text)
        self.assertNotIn("cp -a \"$REF/bowtie2_index", text)
        self.assertIn("flock -n 9", text)
        self.assertLess(text.index("flock -n 9"), text.index('mkdir -p "$WORK"'))
        self.assertIn("flock -s 8", text)
        self.assertIn("flock -x 8", text)
        self.assertIn("gzip -n -c", text)
        self.assertIn('mv "$ART" "$ART_FINAL"', text)
        submit = SUBMIT.read_text(encoding="utf-8")
        self.assertLess(submit.index("reserve_production_slot"), submit.index("output = prepare"))

    def test_two_slots_reserve_and_third_fails_without_analysis_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = SUPPORT.reserve_production_slot(root, "ERR868155", "1" * 32)
            second = SUPPORT.reserve_production_slot(root, "ERR868156", "2" * 32)
            self.assertEqual({first, second}, {1, 2})
            with self.assertRaisesRegex(ValueError, "All 2"):
                SUPPORT.reserve_production_slot(root, "ERR868154", "3" * 32)
            self.assertFalse((root / "analyses" / "ERR868154").exists())
            SUPPORT.release_production_slot(root, first, "ERR868155", "1" * 32)
            self.assertEqual(SUPPORT.reserve_production_slot(root, "ERR868154", "3" * 32), first)

    def test_locked_ledger_append_is_unique(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); ledger = root / "production_ledger.tsv"
            row = root / "row.tsv"
            ledger.write_text("analysis_id\tvalidation_status\nERR868151\tPASS\n", encoding="utf-8")
            row.write_text("analysis_id\tvalidation_status\nERR868152\tPASS\n", encoding="utf-8")
            SUPPORT.append_ledger_row(ledger, row)
            self.assertEqual([item["analysis_id"] for item in SUPPORT.table(ledger)],
                             ["ERR868151", "ERR868152"])
            with self.assertRaisesRegex(ValueError, "already exists"):
                SUPPORT.append_ledger_row(ledger, row)

    def test_control_readers_share_and_exclusive_writer_waits(self):
        with tempfile.TemporaryDirectory() as temporary, \
             SUPPORT.control_lock(Path(temporary), "ERR868150", exclusive=False) as lock_path:
            shared = subprocess.run(["flock", "-sn", str(lock_path), "true"], check=False)
            exclusive = subprocess.run(["flock", "-xn", str(lock_path), "true"], check=False)
            self.assertEqual(shared.returncode, 0)
            self.assertNotEqual(exclusive.returncode, 0)

    def test_synthetic_workers_share_control_reject_third_and_release_locks(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); control_locks = root / ".control_locks"
            control_locks.mkdir(); control_path = control_locks / "ERR868150.lock"
            control_path.touch()

            def launch(slot, analysis, token, ready, delay, status):
                lock_path, _ = SUPPORT.slot_paths(root, slot)
                script = f'''set -eu
exec 9>"{lock_path}"
flock -n 9
"{sys.executable}" "{SUPPORT_PATH}" validate-slot --production-root "{root}" --slot "{slot}" --analysis "{analysis}" --token "{token}"
cleanup() {{ "{sys.executable}" "{SUPPORT_PATH}" release-slot --production-root "{root}" --slot "{slot}" --analysis "{analysis}" --token "{token}"; }}
trap cleanup EXIT
trap 'exit 143' TERM INT
exec 8>"{control_path}"
flock -s 8
touch "{ready}"
sleep "{delay}"
exit "{status}"
'''
                return subprocess.Popen(["bash", "-c", script], start_new_session=True)

            first_token, second_token = "1" * 32, "2" * 32
            first = SUPPORT.reserve_production_slot(root, "ERR868155", first_token)
            second = SUPPORT.reserve_production_slot(root, "ERR868156", second_token)
            ready_one, ready_two = root / "ready-one", root / "ready-two"
            one = launch(first, "ERR868155", first_token, ready_one, 1, 0)
            two = launch(second, "ERR868156", second_token, ready_two, 1, 0)
            for _ in range(100):
                if ready_one.exists() and ready_two.exists(): break
                time.sleep(0.02)
            self.assertTrue(ready_one.exists() and ready_two.exists())
            with self.assertRaisesRegex(ValueError, "All 2"):
                SUPPORT.reserve_production_slot(root, "ERR868154", "3" * 32)
            self.assertFalse((root / "analyses" / "ERR868154").exists())
            self.assertNotEqual(subprocess.run(
                ["flock", "-xn", str(control_path), "true"], check=False).returncode, 0)
            self.assertEqual(one.wait(timeout=5), 0)
            self.assertEqual(two.wait(timeout=5), 0)

            failure_token = "4" * 32
            failure_slot = SUPPORT.reserve_production_slot(root, "ERR868154", failure_token)
            failed = launch(failure_slot, "ERR868154", failure_token,
                            root / "ready-failure", 0, 7)
            self.assertEqual(failed.wait(timeout=5), 7)
            self.assertFalse(SUPPORT.slot_paths(root, failure_slot)[1].exists())

            signal_token = "5" * 32
            signal_slot = SUPPORT.reserve_production_slot(root, "ERR868154", signal_token)
            ready_signal = root / "ready-signal"
            signalled = launch(signal_slot, "ERR868154", signal_token, ready_signal, 30, 0)
            for _ in range(100):
                if ready_signal.exists(): break
                time.sleep(0.02)
            self.assertTrue(ready_signal.exists())
            os.killpg(signalled.pid, signal.SIGTERM)
            self.assertEqual(signalled.wait(timeout=5), 143)
            self.assertFalse(SUPPORT.slot_paths(root, signal_slot)[1].exists())

    def test_control_cleanup_requires_all_expected_analyses(self):
        with tempfile.TemporaryDirectory() as temporary, \
             mock.patch.object(SUPPORT, "verify_analysis", side_effect=ValueError("missing")) as verify:
            with self.assertRaisesRegex(ValueError, "missing"):
                SUPPORT.control_cleanup(Path(temporary), "SRR20770287", False)
            self.assertEqual(verify.call_args.args[1], "SRR20770294")

    def test_finalize_records_observed_fragment_size(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary); artifacts = output / "artifacts"; peaks = artifacts / "peaks"
            peaks.mkdir(parents=True)
            selected = SUPPORT.select_analysis("SRR20770297", "SRR20770287")
            (peaks / "SRR20770297_peaks.narrowPeak").write_text("", encoding="utf-8")
            for name in ("peak_qc.json", "provenance.json"):
                (artifacts / name).write_text("{}\n", encoding="utf-8")
            (artifacts / "parameters.json").write_text(json.dumps({"fragment_size_bp": 147}) + "\n", encoding="utf-8")
            (artifacts / "macs3.log").write_text("", encoding="utf-8")
            SUPPORT.finalize(selected, artifacts, output, "2026-09-21T00:00:00Z", "123")
            rows = SUPPORT.table(output / "production_run_manifest.tsv")
            self.assertTrue(rows)
            self.assertEqual({row["fragment_size_bp"] for row in rows}, {"147"})

if __name__ == "__main__": unittest.main()
