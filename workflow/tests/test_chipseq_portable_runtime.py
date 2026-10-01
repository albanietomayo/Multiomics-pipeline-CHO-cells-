"""Portable adapter regression tests. No sequencing production is executed."""
import gzip
import hashlib
import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("portable_test", ROOT / "workflow/scripts/chipseq_portable_runtime.py")
P = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(P)


class PortableRuntimeTests(unittest.TestCase):
    def test_exact_alignment_command_and_read_groups(self):
        self.assertEqual(P.bowtie2_command("SRR1", "reads.gz", "index"), [
            "bowtie2", "--local", "--very-sensitive-local", "--seed", "0", "-p", "6",
            "-x", "index", "-U", "reads.gz", "--rg-id", "SRR1", "--rg", "SM:SRR1",
            "--rg", "LB:SRR1", "--rg", "PL:ILLUMINA"])

    def test_duplicate_contract(self):
        cmd = P.picard_command("raw.bam", "marked.bam", "metrics", "scratch")
        for flag in ("REMOVE_DUPLICATES=false", "REMOVE_SEQUENCING_DUPLICATES=false", "READ_NAME_REGEX=null",
                     "DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES", "VALIDATION_STRINGENCY=STRICT",
                     "CREATE_INDEX=false", "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=256"):
            self.assertIn(flag, cmd)
        self.assertEqual(cmd[:4], ["picard", "-Xmx12000m", "-XX:ActiveProcessorCount=4", "MarkDuplicates"])

    def test_fastq_manifest_has_unique_twenty_two_verified_runs(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "manifest.tsv"; P.fastq_manifest(path)
            rows = P.rows(path)
            self.assertEqual(len(rows), 22)
            self.assertEqual(len({r["run_accession"] for r in rows}), 22)
            self.assertTrue(all(int(r["fastq_bytes"]) > 0 and len(r["fastq_md5"]) == 32 for r in rows))

    def test_config_blocks_frozen_output_overwrite(self):
        import yaml
        with tempfile.TemporaryDirectory() as temp:
            cfg = yaml.safe_load((ROOT / "config/chipseq.yaml").read_text())
            cfg["chipseq"]["output_root"] = "benchmarks/chipseq"
            f = Path(temp) / "bad.yaml"; f.write_text(yaml.safe_dump(cfg))
            with self.assertRaisesRegex(ValueError, "overwrite"):
                P.configuration(f)

    def test_config_cannot_redirect_frozen_pairing(self):
        import yaml
        with tempfile.TemporaryDirectory() as temp:
            cfg = yaml.safe_load((ROOT / "config/chipseq.yaml").read_text())
            cfg["chipseq"]["analysis_plan"] = "other.tsv"
            f = Path(temp) / "bad.yaml"; f.write_text(yaml.safe_dump(cfg))
            with self.assertRaisesRegex(ValueError, "cannot be redirected"):
                P.configuration(f)

    def test_reconstructed_pairing_difference_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            f = Path(temp) / "analysis.tsv"
            text = P.PRODUCTION.ANALYSIS_PLAN.read_text().replace("SRR20770287", "SRR20770295")
            f.write_text(text)
            with self.assertRaisesRegex(ValueError, "differs"):
                P.compare_plans(f, P.PRODUCTION.PROCESSING_PLAN)

    def test_same_plans_validate(self):
        self.assertEqual(P.compare_plans(P.PRODUCTION.ANALYSIS_PLAN, P.PRODUCTION.PROCESSING_PLAN)["physical_runs"], 22)

    def test_query_identity_required_before_materialization(self):
        with tempfile.TemporaryDirectory() as temp:
            f = Path(temp) / "plan.tsv"; f.write_bytes(b"wrong")
            with self.assertRaises(ValueError):
                P.metadata_record(f, temp, "ERX1")
            self.assertFalse((Path(temp) / "records").exists())

    def test_metadata_replay_uses_sanitized_snapshot_without_network(self):
        replay = mock.Mock()
        with mock.patch.object(P.subprocess, "run") as run:
            with mock.patch.object(P, "module", return_value=replay):
                P.metadata_record("plan", "output", "ERX1")
            replay.materialize.assert_called_once_with("plan", "output", "ERX1")
            run.assert_not_called()

    def test_run_propagates_tool_failure(self):
        with mock.patch.object(P.subprocess, "run", side_effect=subprocess.CalledProcessError(2, "tool")):
            with self.assertRaises(subprocess.CalledProcessError):
                P.run(["tool"])

    def test_deterministic_gzip_is_lossless(self):
        with tempfile.TemporaryDirectory() as temp:
            raw = Path(temp) / "track.bdg"; raw.write_bytes(b"chr1\t0\t10\t2\n")
            outputs = [Path(temp) / (n + ".gz") for n in ("one", "two")]
            for f in outputs:
                P.run(["gzip", "-n", "-c", raw], output_binary=f)
            self.assertEqual(outputs[0].read_bytes(), outputs[1].read_bytes())
            self.assertEqual(gzip.decompress(outputs[0].read_bytes()), raw.read_bytes())
            P.PRODUCTION.compression_record([raw], [outputs[0]], Path(temp) / "compression.json")

    def test_empty_peak_set_overlap_is_zero(self):
        with tempfile.TemporaryDirectory() as temp:
            f = Path(temp) / "empty.bed"; f.touch()
            with mock.patch.object(P.subprocess, "Popen") as proc:
                self.assertEqual(P.overlap_count("missing.bam", f), 0)
                proc.assert_not_called()

    def test_overlap_failure_is_not_a_frip_result(self):
        import io
        with tempfile.TemporaryDirectory() as temp:
            f = Path(temp) / "peaks"; f.write_text("chr1\t0\t1\n")
            process = mock.Mock(stdout=io.StringIO("read\n")); process.wait.return_value = 2
            with mock.patch.object(P.subprocess, "Popen", return_value=process):
                with self.assertRaises(RuntimeError):
                    P.overlap_count("bam", f)

    def test_slots_release_after_exception(self):
        import fcntl
        with tempfile.TemporaryDirectory() as temp:
            cfg = {"output_root": temp, "max_parallel_analyses": 1}
            with self.assertRaisesRegex(RuntimeError, "synthetic"):
                with P.production_slot(cfg):
                    with (Path(temp) / ".portable.slot.0.lock").open("a+") as probe:
                        with self.assertRaises(BlockingIOError):
                            fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    raise RuntimeError("synthetic")
            with (Path(temp) / ".portable.slot.0.lock").open("a+") as probe:
                fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def test_external_policy_never_falls_back(self):
        row = P.PRODUCTION.select_analysis("ERR868151")["analysis"]
        with self.assertRaisesRegex(ValueError, "no fallback"):
            P.DYNAMIC.resolve_fragment_size(row)

    def test_wrong_tool_version_stops_execution(self):
        with mock.patch.object(P.subprocess, "check_output", return_value="bowtie2 version 2.5.4"):
            with self.assertRaisesRegex(ValueError, "Unexpected tool version"):
                P.assert_versions()

    def test_corrupt_fastq_stops_before_any_processing(self):
        with tempfile.TemporaryDirectory() as temp:
            raw = Path(temp) / "raw.gz"; raw.write_bytes(b"corrupt")
            plan = P.PRODUCTION.select_analysis("SRR20770291")
            with mock.patch.object(P, "run") as tool:
                with self.assertRaisesRegex(ValueError, "FASTQ no longer matches"):
                    P.process_run(plan, "ip", raw, Path(temp) / "work")
                tool.assert_not_called()
            self.assertFalse((Path(temp) / "work").exists())

    def test_reference_inventory_requires_frozen_identity(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "reference_provenance.json").write_text("{}")
            (root / "reference_inventory.json").write_text(json.dumps({"identity": {}, "files": []}))
            with self.assertRaisesRegex(ValueError, "identity mismatch"):
                P.reference_validate(root)

    def test_snakemake_empty_output_scaffold_can_be_claimed(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / "analysis"
            (out / "artifacts/peaks").mkdir(parents=True)
            P.remove_output_scaffold(out, ("", "artifacts", "artifacts/peaks"))
            self.assertFalse(out.exists())

    def test_output_scaffold_never_deletes_partial_data(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / "analysis"; out.mkdir()
            f = out / "partial.bam"; f.write_bytes(b"partial")
            with self.assertRaisesRegex(ValueError, "collision"):
                P.remove_output_scaffold(out)
            self.assertEqual(f.read_bytes(), b"partial")

    def test_output_scaffold_never_follows_symlinks(self):
        with tempfile.TemporaryDirectory() as temp:
            real = Path(temp) / "reference"; real.mkdir()
            alias = Path(temp) / "alias"; alias.symlink_to(real, target_is_directory=True)
            with self.assertRaisesRegex(ValueError, "collision"):
                P.remove_output_scaffold(alias)
            self.assertTrue(real.exists())

    def test_output_scaffold_rejects_unexpected_empty_directories(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / "analysis"; (out / "unknown").mkdir(parents=True)
            with self.assertRaisesRegex(ValueError, "collision"):
                P.remove_output_scaffold(out)
            self.assertTrue((out / "unknown").exists())

    def test_runtime_per_analysis_not_shared_config(self):
        text = (ROOT / "workflow/scripts/chipseq_portable_runtime.py").read_text()
        self.assertNotIn('ROOT / "config/chipseq_alignment_inputs.json"', text)
        self.assertNotIn('ROOT / "config/chipseq_filtering_inputs.json"', text)


if __name__ == "__main__": unittest.main()
