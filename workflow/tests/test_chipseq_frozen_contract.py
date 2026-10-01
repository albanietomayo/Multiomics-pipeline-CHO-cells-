"""Regression protection for frozen evidence and unchanged public science."""
import collections
import csv
import hashlib
import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
E = ROOT / "provenance/chipseq/production_2026-09-23"
spec = importlib.util.spec_from_file_location("frozen_runtime", ROOT / "workflow/scripts/chipseq_portable_runtime.py")
P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)


class FrozenContractTests(unittest.TestCase):
    def test_frozen_input_hashes(self): P.check_manifest(E / "frozen_inputs.sha256")
    def test_existing_frozen_output_hashes(self): P.check_manifest(E / "frozen_outputs.sha256")

    def test_counts_and_marks(self):
        self.assertEqual(P.validate_contract(), {"schema_version": 1, "scope": "frozen_plan_validation", "analyses": 18,
                                                "physical_runs": 22, "shared_inputs": 4, "narrow": 8, "broad": 10})
        rows, _ = P.PRODUCTION.validated_plans()
        self.assertEqual({r["declared_target"] for r in rows}, {"H3K4me3", "H3K27ac", "H3K9me3", "H3K4me1", "H3K27me3", "H3K36me3"})

    def test_source_mapping_and_environment_mapping(self):
        records = P.rows(E / "analysis_sources.tsv")
        self.assertEqual(len(records), 18)
        expected = {
            "48d80c961acfd188033db9e1cf2972222fbc48a5": {"SRR20770294", "SRR20770296", "SRR20770297"},
            "fb62508d51efa41ab73067ba7ca57966516cd5f9": {"ERR868151", "ERR868152", "ERR868153", "ERR868154", "SRR20770291", "SRR20770292", "SRR20770293"},
            "290c35dcec4b76c258b2788a4adbc6940654d04f": {"ERR868155", "ERR868156", "ERR868171", "ERR868172", "ERR868173", "ERR868174", "ERR868175", "ERR868176"}}
        for commit, ids in expected.items():
            self.assertEqual({r["analysis_id"] for r in records if r["production_commit"] == commit}, ids)
        envs = {hashlib.sha256(f.read_bytes()).hexdigest() for f in (E / "environments").glob("*.txt")}
        self.assertEqual(len(envs), 4)
        self.assertEqual({r["environment_sha256"] for r in records}, envs)

    def test_table_a9_and_combined_pins(self):
        import yaml
        expected = {"alignment": {"bowtie2": "2.5.5", "samtools": "1.24"}, "duplicates": {"picard": "3.5.0"},
                    "filtering": {"pysam": "0.24.0"}, "phantompeakqualtools": {"phantompeakqualtools": "1.2.2"},
                    "peak_calling": {"macs3": "3.0.4", "bedtools": "2.31.1"}}
        combined = yaml.safe_load((ROOT / "workflow/envs/chipseq_production.yaml").read_text())["dependencies"]
        for stage, packages in expected.items():
            deps = yaml.safe_load((ROOT / f"workflow/envs/chipseq_{stage}.yaml").read_text())["dependencies"]
            for name, version in packages.items():
                self.assertIn(name + "=" + version, deps)
                self.assertIn(name + "=" + version, combined)
                for f in (E / "environments").glob("*.txt"):
                    self.assertIn("/" + name + "-" + version + "-", f.read_text())

    def test_reference_filtering_and_limitations(self):
        c = P.read_json(E / "frozen_contract.json")
        self.assertEqual(c["reference"], {"assembly_name": "CriGri-PICRH-1.0", "accession": "GCF_003668045.3", "mitochondrial": "NC_007936.1", "nuclear_span_bp": 2366634374})
        policy = P.read_json(ROOT / "config/chipseq_filtering.json")
        self.assertEqual((policy["min_mapq"], policy["exclude_mapq_255"], policy["exclude_flags"]), (30, True, 3844))
        self.assertFalse(policy["optical_duplicate_detection"])
        self.assertEqual(c["blacklist"], "not_performed")
        self.assertEqual(c["replicate_concordance_idr"], "not_performed")
        self.assertEqual(c["macs3"], {"qvalue": 0.01, "broad_cutoff": 0.1, "keep_dup": "all", "scale_to": "small", "spmr": True})

    def test_historical_fragment_sizes(self):
        values = {r["analysis_id"]: int(r["fragment_size_bp"]) for r in P.rows(E / "analysis_sources.tsv")}
        expected = dict(zip(["ERR868151", "ERR868152", "ERR868153", "ERR868154", "ERR868155", "ERR868156", "ERR868171", "ERR868172", "ERR868173", "ERR868174", "ERR868175", "ERR868176"], [160,155,145,145,215,185,155,150,190,160,160,155]))
        self.assertEqual({k:v for k,v in values.items() if k.startswith("ERR")}, expected)
        self.assertTrue(all(v == 147 for k,v in values.items() if k.startswith("SRR")))


if __name__ == "__main__": unittest.main()
