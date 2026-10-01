#!/usr/bin/env python3
"""Synthetic tests for dynamic Phase 6B contracts; no production tools are run."""

import importlib.util
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "dynamic_peak",
    ROOT / "workflow/scripts/chipseq_peak_calling_dynamic.py",
)
D = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(D)


def plan(study="PRJNA865478", mode="narrow", run="SRR1"):
    external = study == "PRJEB9291"
    return {
        "analysis_id": run,
        "study_accession": study,
        "ip_run_accession": run,
        "control_run_accession": "SRR9",
        "declared_target": "H3K4me3",
        "peak_mode": mode,
        "format": "BAM",
        "qvalue": "0.01",
        "broad_cutoff": "0.1" if mode == "broad" else "",
        "keep_dup": "all",
        "scale_to": "small",
        "spmr": "true",
        "store_bdg": "true",
        "call_summits": "false",
        "cutoff_analysis": "false",
        "effective_genome_size_policy": D.GENOME_SIZE_POLICY,
        "fragment_size_policy": "phantompeakqualtools" if external else "fixed",
        "fragment_size_bp": "" if external else "147",
    }


def estimate(run="ERR1", fragment=151, analysis=None, bam_path="/tmp/ERR1.filtered.bam"):
    value = {
        "schema_version": 2,
        "scope": "phantompeakqualtools_fragment_estimate",
        "analysis_id": analysis or run,
        "run_accession": run,
        "ip_bam_path": bam_path,
        "ip_bam_sha256": "0" * 64,
        "source_filename": run + ".filtered.bam",
        "num_reads": 100,
        "num_reads_semantics": (
            "PhantomPeakQualTools-reported eligible reads; preserved separately from "
            "the filtering retained-read count because formal semantics do not justify equality"
        ),
        "selected_fragment_candidate_bp": fragment,
        "fragment_size_bp": fragment,
        "fragment_candidates_bp": [fragment],
        "fragment_candidate_correlations": [0.2],
        "fragment_correlation": 0.2,
        "phantom_peak": 36.0,
        "phantom_correlation": 0.1,
        "minimum_cross_correlation_shift": -100.0,
        "minimum_cross_correlation": 0.05,
        "nsc": 1.5,
        "rsc": 1.2,
        "quality_tag": 2,
        "quality_tag_interpretation": "descriptive_tool_output_not_a_universal_pass_fail_threshold",
        "source_result_path": "/tmp/phantompeakqualtools.tsv",
        "source_result_sha256": "0" * 64,
    }
    value["record_sha256"] = D.canonical_sha256(value)
    return value


class DynamicPeakCallingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.fixture_number = 0

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name, text):
        path = self.root / name
        path.write_text(text, encoding="utf-8")
        return path

    def reference_metadata(self):
        return {
            "nuclear_accession": D.ASSEMBLY_ACCESSION,
            "nuclear_sha256": "1" * 64,
            "mitochondrial_accession": D.MTDNA,
            "mitochondrial_sha256": "2" * 64,
            "mapping_sha256": "3" * 64,
            "nuclear_sequences": 1,
            "purpose": "ChIP-seq nuclear and mitochondrial mapping",
        }

    def test_fai_derives_nuclear_span_and_reference_provenance(self):
        fai = self.write("ref.fai", "chr1\t1000\t0\t10\t11\nNC_007936.1\t50\t0\t10\t11\n")
        result = D.validate_reference(fai, self.reference_metadata(), D.sha256(fai))
        self.assertEqual(result["nuclear_span"], 1000)
        self.assertIn("not a mappability-derived", result["genome_size_description"])

    def test_fai_malformed_duplicate_missing_mt_and_invalid_lengths(self):
        cases = {
            "malformed": ("chr1\n", "Malformed"),
            "truncated": ("chr1\t10\t0\t10\n", "five"),
            "duplicate": ("chr1\t10\t0\t10\t11\nchr1\t20\t20\t10\t11\nNC_007936.1\t5\t50\t5\t6\n", "Duplicate"),
            "missing_mt": ("chr1\t10\t0\t10\t11\n", "absent"),
            "zero": ("chr1\t0\t0\t10\t11\nNC_007936.1\t5\t20\t5\t6\n", "Invalid FAI length"),
            "negative": ("chr1\t10\t-1\t10\t11\nNC_007936.1\t5\t20\t5\t6\n", "Invalid FAI offset"),
            "geometry": ("chr1\t10\t0\t10\t9\nNC_007936.1\t5\t20\t5\t6\n", "geometry"),
        }
        for name, (content, message) in cases.items():
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, message):
                D.parse_fai(self.write(name + ".fai", content))

    def test_reference_hash_tampering_fails(self):
        fai = self.write("ref.fai", "chr1\t10\t0\t10\t11\nNC_007936.1\t5\t20\t5\t6\n")
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            D.validate_reference(fai, self.reference_metadata(), "0" * 64)

    def test_reference_identity_fields_are_mandatory(self):
        fai = self.write("required.fai", "chr1\t10\t0\t10\t11\nNC_007936.1\t5\t20\t5\t6\n")
        metadata = self.reference_metadata()
        del metadata["mapping_sha256"]
        with self.assertRaisesRegex(ValueError, "mandatory"):
            D.validate_reference(fai, metadata, D.sha256(fai))

    def phantom(self, fragment="150", correlation="0.25", filename="ERR1.filtered.bam", extra=""):
        line = "\t".join([
            filename, "1000", fragment, correlation, "36", "0.10", "-100",
            "0.05", "1.5", "1.2", "2",
        ])
        return self.write("phantom.tsv", line + "\n" + extra)

    def test_formal_phantom_result_is_preserved(self):
        result = D.parse_phantompeak_table(self.phantom(), "ERR1")
        self.assertEqual(result["fragment_size_bp"], 150)
        self.assertEqual(result["fragment_candidates_bp"], [150])
        self.assertEqual(result["nsc"], 1.5)
        self.assertEqual(result["rsc"], 1.2)
        self.assertEqual(result["quality_tag"], 2)
        self.assertEqual(result["minimum_cross_correlation_shift"], -100.0)

    def test_phantom_fragment_fail_closed_cases(self):
        for value in ("NA", "0", "-1", "1.5", "bad", "", "147,,150"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                D.parse_phantompeak_table(self.phantom(fragment=value), "ERR1")

    def test_phantom_multi_candidate_formal_output(self):
        result = D.parse_phantompeak_table(
            self.phantom(fragment="147,150", correlation="0.25,0.24"), "ERR1"
        )
        self.assertEqual(result["fragment_candidates_bp"], [147, 150])
        self.assertEqual(result["fragment_candidate_correlations"], [0.25, 0.24])
        self.assertEqual(result["selected_fragment_candidate_bp"], 147)

    def test_phantom_candidate_pairing_fail_closed(self):
        cases = (("147,150", "0.25"), ("147", ""), ("147", "bad"),
                 ("147", "nan"), ("147", "inf"))
        for fragments, correlations in cases:
            with self.subTest(fragments=fragments, correlations=correlations), self.assertRaises(ValueError):
                D.parse_phantompeak_table(
                    self.phantom(fragment=fragments, correlation=correlations), "ERR1"
                )

    def test_phantom_ambiguous_and_mismatched_fail(self):
        extra = self.phantom().read_text()
        with self.assertRaisesRegex(ValueError, "exactly one"):
            D.parse_phantompeak_table(self.phantom(extra=extra), "ERR1")
        with self.assertRaisesRegex(ValueError, "does not match"):
            D.parse_phantompeak_table(self.phantom(filename="ERR2.filtered.bam"), "ERR1")

    def test_fixed_and_external_fragment_policies(self):
        fixed = D.macs3_contract(plan(), 1000)
        self.assertEqual(fixed["fragment_size_bp"], 147)
        external_row = plan("PRJEB9291", run="ERR1")
        external = D.macs3_contract(external_row, 1000, estimate())
        self.assertEqual(external["fragment_size_bp"], 151)
        self.assertEqual(external["fragment_size_source"], "phantompeakqualtools_external_estimate")

    def test_no_147_fallback_and_estimate_is_analysis_specific(self):
        row = plan("PRJEB9291", run="ERR1")
        with self.assertRaisesRegex(ValueError, "no fallback"):
            D.macs3_contract(row, 1000)
        with self.assertRaisesRegex(ValueError, "does not match"):
            D.macs3_contract(row, 1000, estimate("ERR2", 147))

    def test_tampered_fragment_result_fails_closed(self):
        result = estimate()
        result["fragment_size_bp"] = 147
        with self.assertRaisesRegex(ValueError, "tampered"):
            D.macs3_contract(plan("PRJEB9291", run="ERR1"), 1000, result)

    def test_phantom_result_is_bound_to_analysis_bam_and_raw_table(self):
        bam = self.write("ERR1.filtered.bam", "bam\n")
        table = self.phantom(fragment="147,150", correlation="0.25,0.24")
        row = plan("PRJEB9291", run="ERR1")
        result = D.parse_phantompeak_table(
            table, "ERR1", bam, row["analysis_id"], D.sha256(bam)
        )
        runtime_run = {
            "run_accession": "ERR1", "bam_path": str(bam.resolve()),
            "bam_sha256": D.sha256(bam), "retained_reads": 1000,
        }
        self.assertEqual(
            D.validate_phantom_result(row, result, runtime_run, None, table), 147
        )
        self.assertEqual(result["fragment_candidates_bp"], [147, 150])

        table.write_text(table.read_text() + "# tamper\n")
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            D.validate_phantom_result(row, result, runtime_run, None, table)

    def test_phantom_identity_and_candidate_tampering_fails(self):
        bam = self.write("ERR1.filtered.bam", "bam\n")
        table = self.phantom(fragment="147,150", correlation="0.25,0.24")
        row = plan("PRJEB9291", run="ERR1")
        runtime_run = {
            "run_accession": "ERR1", "bam_path": str(bam.resolve()),
            "bam_sha256": D.sha256(bam), "retained_reads": 1000,
        }
        original = D.parse_phantompeak_table(
            table, "ERR1", bam, row["analysis_id"], D.sha256(bam)
        )
        mutations = {
            "analysis_id": "wrong", "run_accession": "ERR2",
            "ip_bam_sha256": "f" * 64, "selected_fragment_candidate_bp": 150,
            "fragment_candidates_bp": [150, 147],
            "fragment_candidate_correlations": [0.24, 0.25],
        }
        for key, value in mutations.items():
            changed = dict(original)
            changed[key] = value
            changed["record_sha256"] = D.canonical_sha256({
                name: item for name, item in changed.items() if name != "record_sha256"
            })
            with self.subTest(key=key), self.assertRaises(ValueError):
                D.validate_phantom_result(row, changed, runtime_run, None, table)

    def test_prjna_rejects_phantom_dependency(self):
        with self.assertRaisesRegex(ValueError, "must not consume"):
            D.macs3_contract(plan(), 1000, {"run_accession": "SRR1", "fragment_size_bp": 200})

    def test_narrow_and_broad_output_contracts(self):
        narrow = D.macs3_contract(plan(mode="narrow"), 1000)
        broad = D.macs3_contract(plan(mode="broad"), 1000)
        self.assertEqual(narrow["required_peak_outputs"], ["narrowPeak", "summits"])
        self.assertEqual(broad["required_peak_outputs"], ["broadPeak", "gappedPeak"])
        self.assertNotIn("summits", broad["required_peak_outputs"])
        self.assertIn("--broad", broad["flags"])
        self.assertIn("--keep-dup", broad["flags"])

    def test_peak_validation_modes_and_zero_peaks(self):
        contigs = {"chr1": 1000, D.MTDNA: 50}
        narrow = self.write("x.narrowPeak", "chr1\t10\t20\tp\t1\t.\t2\t3\t4\t5\n")
        broad = self.write("x.broadPeak", "chr1\t10\t30\tp\t1\t.\t2\t3\t4\n")
        gapped = self.write("x.gappedPeak", "chr1\t10\t30\tp\t1\t.\t2\t3\t4\t2\t5,5\t0,15\t2\t3\t4\n")
        empty = self.write("empty.broadPeak", "")
        self.assertEqual(D.parse_peak_file(narrow, "narrow", contigs)["peak_count"], 1)
        self.assertEqual(D.parse_peak_file(broad, "broad", contigs)["peak_count"], 1)
        self.assertEqual(D.parse_peak_file(gapped, "gapped", contigs)["peak_count"], 1)
        self.assertEqual(D.parse_peak_file(empty, "broad", contigs)["peak_count"], 0)

    def test_malformed_unknown_out_of_range_and_mt_peaks(self):
        contigs = {"chr1": 100, D.MTDNA: 50}
        cases = {
            "malformed": ("chr1\t1\t2\n", "Malformed"),
            "unknown": ("chrX\t1\t2\tp\t1\t.\t2\t3\t4\n", "Unknown"),
            "range": ("chr1\t90\t101\tp\t1\t.\t2\t3\t4\n", "out of range"),
            "mt": ("NC_007936.1\t1\t2\tp\t1\t.\t2\t3\t4\n", "Mitochondrial"),
        }
        for name, (content, message) in cases.items():
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, message):
                D.parse_peak_file(self.write(name, content), "broad", contigs)

    def test_frip_edge_cases_and_denominators(self):
        result = D.descriptive_frip(100, 80, 40, 8)
        self.assertEqual(result["ip_frip"], 0.4)
        self.assertEqual(result["ip_to_input_peak_overlap_ratio"], 4.0)
        zero = D.descriptive_frip(0, 0, 0, 0)
        self.assertIsNone(zero["ip_frip"])
        with self.assertRaisesRegex(ValueError, "exceeds"):
            D.descriptive_frip(10, 10, 11, 0)
        with self.assertRaisesRegex(ValueError, "non-negative"):
            D.descriptive_frip(10, 10, -1, 0)

    def test_unsafe_overwrite_rejected(self):
        existing = self.write("result.json", "{}")
        with self.assertRaisesRegex(ValueError, "Refusing unsafe overwrite"):
            D.validate_output_destination(existing)

    def test_snakemake_owned_json_is_atomically_replaceable(self):
        output = self.root / "parameters.json"
        D.atomic_json(output, {"generation": 1})
        D.atomic_json(output, {"generation": 2})
        self.assertEqual(json.loads(output.read_text()), {"generation": 2})
        self.assertEqual(list(self.root.glob("*.part")), [])

    def runtime_fixture(self):
        self.fixture_number += 1
        root = self.root / f"runtime_{self.fixture_number}"
        root.mkdir()
        plan_fields = list(plan().keys()) + ["execution_status"]
        rows = []
        for run in ("SRR1", "SRR2"):
            item = plan(run=run)
            item["execution_status"] = "ready_for_peak_calling"
            rows.append(item)
        plan_path = root / "peak.tsv"
        with plan_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=plan_fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        policy = root / "policy.json"
        policy.write_text("{}\n")
        analysis = root / "analysis.tsv"
        with analysis.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "analysis_id", "ip_run_accession", "control_run_accession",
                    "analysis_status",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({
                    "analysis_id": row["analysis_id"],
                    "ip_run_accession": row["ip_run_accession"],
                    "control_run_accession": row["control_run_accession"],
                    "analysis_status": "ready",
                })
        summary = root / "summary.json"
        summary.write_text(json.dumps({
            "input_sha256": {
                "analysis_plan": D.sha256(analysis),
                "policy": D.sha256(policy),
            },
            "peak_calling_plan_sha256": D.sha256(plan_path),
        }))
        job = root / "job_1"
        job.mkdir()
        (job / "job_status.tsv").write_text("stage\tcompleted\nexit_status\t0\n")
        outputs = job / "outputs"
        (outputs / "reference").mkdir(parents=True)
        fai = outputs / "reference/genome_plus_mt.fa.fai"
        fai.write_text("chr1\t1000\t6\t50\t51\nNC_007936.1\t50\t1030\t50\t51\n")
        provenance = outputs / "reference/reference_provenance.json"
        provenance.write_text(json.dumps(self.reference_metadata()))
        (outputs / "input_provenance.json").write_text(json.dumps({
            "schema_version": 1,
            "source_manifest_sha256": "4" * 64,
            "small_files": {
                "outputs/reference/genome_plus_mt.fa.fai": D.sha256(fai),
                "outputs/reference/reference_provenance.json": D.sha256(provenance),
            },
            "reference_provenance": self.reference_metadata(),
        }))
        (outputs / "filtering_parameters.json").write_text("{}\n")
        roles = {"SRR1": "ip", "SRR2": "ip", "SRR9": "input"}
        for run, role in roles.items():
            folder = outputs / run
            folder.mkdir()
            (folder / "filtered.bam").write_bytes((run + " bam").encode())
            (folder / "filtered.bam.csi").write_bytes((run + " csi").encode())
            validation = {
                "run_accession": run,
                "full_read_verified": True,
                "retained_reads": 10,
                "bam_sha256": D.sha256(folder / "filtered.bam"),
                "index_sha256": D.sha256(folder / "filtered.bam.csi"),
            }
            (folder / "filtered_validation.json").write_text(json.dumps(validation))
            (folder / "filtering_qc.json").write_text(json.dumps({
                "run_accession": run, "role": role, "retained_reads": 10,
            }))
        files = sorted(path for path in outputs.rglob("*") if path.is_file())
        (job / "output.sha256").write_text("".join(
            f"{D.sha256(path)}  {path.relative_to(job).as_posix()}\n" for path in files
        ))
        return plan_path, summary, policy, analysis, job

    def refresh_output_manifest(self, job):
        outputs = job / "outputs"
        files = sorted(path for path in outputs.rglob("*") if path.is_file())
        (job / "output.sha256").write_text("".join(
            f"{D.sha256(path)}  {path.relative_to(job).as_posix()}\n" for path in files
        ))

    def test_runtime_manifest_supports_shared_inputs_and_hash_provenance(self):
        args = self.runtime_fixture()
        runtime = D.build_runtime_manifest(*args)
        self.assertEqual(runtime["analysis_count"], 2)
        self.assertEqual(runtime["run_count"], 3)
        self.assertEqual(runtime["analyses"][0]["control_run_accession"], "SRR9")
        self.assertEqual(runtime["analyses"][1]["control_run_accession"], "SRR9")
        self.assertTrue(runtime["full_bam_and_csi_hashes_verified"])
        self.assertEqual(runtime["reference"]["nuclear_span"], 1000)
        self.assertIn("analysis_plan", runtime["input_sha256"])
        self.assertIn("peak_policy", runtime["input_sha256"])
        self.assertIn("reference_fai", runtime["input_sha256"])

    def test_runtime_tampered_plan_policy_and_filtered_file_fail(self):
        for target_index, replacement in (
            (0, "tampered\n"),
            (2, "tampered\n"),
        ):
            args = list(self.runtime_fixture())
            args[target_index].write_text(replacement)
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                D.build_runtime_manifest(*args)

        args = list(self.runtime_fixture())
        (args[4] / "outputs/SRR1/filtered.bam").write_bytes(b"tampered")
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            D.build_runtime_manifest(*args)

    def test_runtime_incomplete_job_and_missing_manifest_fail(self):
        args = list(self.runtime_fixture())
        (args[4] / "job_status.tsv").write_text(
            "stage\tfiltering\nexit_status\t1\n"
        )
        with self.assertRaisesRegex(ValueError, "not successfully completed"):
            D.build_runtime_manifest(*args)

        args = list(self.runtime_fixture())
        (args[4] / "output.sha256").unlink()
        with self.assertRaisesRegex(ValueError, "no output checksum manifest"):
            D.build_runtime_manifest(*args)

    def test_runtime_missing_bam_and_csi_fail(self):
        for relative, message in (
            ("outputs/SRR1/filtered.bam", "Missing filtered bam"),
            ("outputs/SRR1/filtered.bam.csi", "Missing filtered csi"),
        ):
            args = list(self.runtime_fixture())
            (args[4] / relative).unlink()
            with self.subTest(relative=relative), self.assertRaisesRegex(
                ValueError, message
            ):
                D.build_runtime_manifest(*args)

    def test_runtime_wrong_roles_fail_after_verified_manifest(self):
        args = list(self.runtime_fixture())
        qc_path = args[4] / "outputs/SRR1/filtering_qc.json"
        qc = json.loads(qc_path.read_text())
        qc["role"] = "input"
        qc_path.write_text(json.dumps(qc))
        self.refresh_output_manifest(args[4])
        with self.assertRaisesRegex(ValueError, "non-IP filtering role"):
            D.build_runtime_manifest(*args)

    def test_runtime_reference_identity_mismatch_fails(self):
        args = list(self.runtime_fixture())
        provenance_path = args[4] / "outputs/input_provenance.json"
        provenance = json.loads(provenance_path.read_text())
        provenance["reference_provenance"]["mapping_sha256"] = "9" * 64
        provenance_path.write_text(json.dumps(provenance))
        self.refresh_output_manifest(args[4])
        with self.assertRaisesRegex(ValueError, "provenance identities disagree"):
            D.build_runtime_manifest(*args)

    def test_per_analysis_qc_and_provenance_records(self):
        runtime = D.build_runtime_manifest(*self.runtime_fixture())
        runtime_path = self.root / "runtime.json"
        D.atomic_json(runtime_path, runtime)
        parameters = D.analysis_parameters(runtime_path, "SRR1")
        parameters_path = self.root / "parameters.json"
        D.atomic_json(parameters_path, parameters)
        peaks = self.write("qc.narrowPeak", "chr1\t10\t20\tp\t1\t.\t2\t3\t4\t5\n")
        summits = self.write("summits.bed", "chr1\t15\t16\tp\t1\n")
        paths = [self.write(name, "synthetic\n") for name in (
            "peaks.xls", "treat.bdg", "control.bdg", "macs.log", "version.txt",
        )]
        count_paths = [self.write(name, value) for name, value in (
            ("ip_total", "10\n"), ("input_total", "10\n"),
            ("ip_overlap", "4\n"), ("input_overlap", "1\n"),
        )]
        qc_path, provenance_path = self.root / "qc.json", self.root / "provenance.json"
        qc = D.summarize_analysis(
            runtime_path, parameters_path, "SRR1", peaks, summits, *paths,
            *count_paths, qc_path, provenance_path,
        )
        provenance = json.loads(provenance_path.read_text())
        self.assertEqual(qc["ip_frip"], 0.4)
        self.assertEqual(qc["interpretation"], "descriptive_only_no_universal_pass_fail_threshold")
        self.assertEqual(provenance["analysis_id"], "SRR1")
        self.assertIn("peak_policy", provenance["input_sha256"])
        self.assertIn("qc.narrowPeak", provenance["output_sha256"])

    def test_macs3_parameter_tampering_fails_before_subprocess(self):
        runtime = D.build_runtime_manifest(*self.runtime_fixture())
        runtime_path = self.root / "runtime_for_macs.json"
        D.atomic_json(runtime_path, runtime)
        baseline = D.analysis_parameters(runtime_path, "SRR1")
        mutations = {
            "fragment_size_bp": 999,
            "macs3_genome_size": 999,
            "qvalue": 0.5,
            "scale_to": "large",
            "keep_dup": "1",
            "peak_mode": "broad",
            "flags": ["--broad", "--broad-cutoff", "0.9"],
            "ip_run_accession": "SRR2",
            "control_run_accession": "SRR2",
        }
        for index, (key, value) in enumerate(mutations.items()):
            parameters = dict(baseline)
            parameters[key] = value
            parameter_path = self.root / f"tampered_{index}.json"
            D.atomic_json(parameter_path, parameters)
            with self.subTest(key=key), mock.patch.object(D.subprocess, "run") as run:
                with self.assertRaisesRegex(ValueError, "reconstructed verified contract"):
                    D.run_macs3(
                        runtime_path, parameter_path, "SRR1", self.root / f"peaks_{index}",
                        self.root / f"log_{index}", self.root / f"version_{index}",
                    )
                run.assert_not_called()


if __name__ == "__main__":
    unittest.main(verbosity=2)
