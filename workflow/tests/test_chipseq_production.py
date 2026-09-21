#!/usr/bin/env python3
import importlib.util
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
SUPPORT_PATH = ROOT / "workflow/scripts/chipseq_production.py"
SBATCH = ROOT / "workflow/slurm/chipseq_production.sbatch"

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
        self.assertIn("gzip -n -c", text)

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
