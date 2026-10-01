#!/usr/bin/env python3
"""Synthetic regression tests; run with the chipseq_filtering Conda environment."""
import importlib.util
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("support", ROOT / "workflow/scripts/chipseq_filtering_support.py")
S = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(S)
RUN = "SRR20770297"
MT = "NC_007936.1"


class FilteringTests(unittest.TestCase):
    def setUp(self):
        self.previous = Path.cwd()
        self.temporary = tempfile.TemporaryDirectory()
        os.chdir(self.temporary.name)
        Path("config").mkdir()
        shutil.copyfile(ROOT / S.CONFIG, S.CONFIG)
        S.dump(S.PLAN, {"configuration_sha256": S.sha(S.CONFIG), "runs": [{
            "run_accession": RUN, "role": "ip", "input_reads": 10, "mapped_reads": 9,
            "nuclear_mapq_ge_threshold": 5}]})
        self.header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                       "SQ": [{"SN": "chr1", "LN": 10000}, {"SN": MT, "LN": 1000}],
                       "RG": [{"ID": RUN, "SM": RUN, "LB": RUN, "PL": "ILLUMINA"}]}
        # Ordered records: overlapping reasons must be counted only once in the flow.
        records = [("keep30", 0, 0, 30), ("keep_reverse", 16, 0, 42),
                   ("duplicate", 1024, 0, 42), ("low_duplicate", 1024, 0, 10),
                   ("qc_duplicate", 1536, 0, 42), ("unknown", 0, 0, 255),
                   ("low29", 0, 0, 29), ("mt_duplicate", 1024, 1, 42),
                   ("mt", 0, 1, 42), ("unmapped", 4, -1, 0)]
        with pysam.AlignmentFile("marked.bam", "wb", header=self.header) as handle:
            for index, (name, flag, tid, mapq) in enumerate(records):
                read = pysam.AlignedSegment(handle.header)
                read.query_name = name
                read.query_sequence = "A" * 25
                read.query_qualities = pysam.qualitystring_to_array("I" * 25)
                read.flag, read.reference_id, read.mapping_quality = flag, tid, mapq
                read.reference_start = index * 30 if tid >= 0 else -1
                read.cigarstring = "25M" if tid >= 0 else None
                read.set_tag("RG", RUN)
                handle.write(read)
        Path("picard.txt").write_text(
            "## METRICS CLASS\tpicard.sam.DuplicationMetrics\n"
            "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tPERCENT_DUPLICATION\n"
            f"{RUN}\t9\t0\t4\t0\t0.444444\n\n")

    def tearDown(self):
        os.chdir(self.previous)
        self.temporary.cleanup()

    def filter(self):
        S.filter_bam(RUN, "marked.bam", "picard.txt", "filtered.bam", "qc.json", "flow.tsv")

    def test_overlap_accounting_mapq_boundary_and_output(self):
        self.filter()
        qc = S.load("qc.json")
        self.assertEqual(qc["retained_reads"], 2)
        self.assertEqual([qc["removed_" + k] for k in S.REASONS], [0, 1, 1, 2, 1, 2, 1])
        self.assertEqual(qc["independent_duplicate"], 4)
        with pysam.AlignmentFile("filtered.bam", "rb") as bam:
            self.assertEqual([r.query_name for r in bam], ["keep30", "keep_reverse"])
        pysam.index("-c", "filtered.bam")
        S.verify(RUN, "filtered.bam", "qc.json", "verification.json")
        self.assertTrue(S.load("verification.json")["full_read_verified"])

    def test_picard_count_disagreement_is_rejected(self):
        Path("picard.txt").write_text(Path("picard.txt").read_text().replace("\t4\t0\t", "\t3\t0\t"))
        with self.assertRaisesRegex(ValueError, "duplicate counts disagree"):
            self.filter()
        self.assertFalse(Path("filtered.bam").exists())

    def test_upstream_count_disagreement_is_rejected(self):
        plan = S.load(S.PLAN)
        plan["runs"][0]["input_reads"] = 11
        S.dump(S.PLAN, plan)
        with self.assertRaisesRegex(ValueError, "upstream alignment counts"):
            self.filter()

    def test_paired_rejected_nonprimary_priority(self):
        with self.assertRaisesRegex(ValueError, "Paired reads"):
            S.failures(1, "chr1", 42, MT, 30)
        failed = S.failures(256 | 2048 | 1024, MT, 10, MT, 30)
        self.assertEqual(next(k for k in S.REASONS if failed[k]), "nonprimary")

    def test_reference_dictionary_and_read_group(self):
        Path("ref.fai").write_text("chr1\t10000\t0\t0\t0\n" + MT + "\t1000\t0\t0\t0\n")
        S.validate_header("marked.bam", RUN, "ref.fai", MT)
        Path("ref.fai").write_text("chr1\t9999\t0\t0\t0\n" + MT + "\t1000\t0\t0\t0\n")
        with self.assertRaisesRegex(ValueError, "sequence dictionary"):
            S.validate_header("marked.bam", RUN, "ref.fai", MT)

    def test_changed_bam_hash_rejected(self):
        digest = S.sha("marked.bam")
        with Path("marked.bam").open("ab") as handle:
            handle.write(b"unexpected")
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            S.check("marked.bam", digest)

    def test_final_verification_rejects_contamination(self):
        self.filter()
        shutil.copyfile("marked.bam", "filtered.bam")
        pysam.index("-c", "filtered.bam")
        with self.assertRaisesRegex(ValueError, "excluded record"):
            S.verify(RUN, "filtered.bam", "qc.json", "verification.json")

    def test_prepare_anchors_reports_to_dynamic_alignment_manifest(self):
        cfg = S.load(S.CONFIG)

        pointer = Path(
            cfg["alignment_latest_submission"]
        )

        parent = pointer.parent

        parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        submission = (
            parent
            / "submission_TEST001"
        )

        job_id = "999001"

        job = (
            submission
            / f"job_{job_id}"
        )

        job.mkdir(
            parents=True,
            exist_ok=True,
        )

        (
            submission
            / "job_id.txt"
        ).write_text(
            job_id + "\n"
        )

        (
            job
            / "job_status.tsv"
        ).write_text(
            "stage\tcompleted\n"
            "exit_status\t0\n"
        )

        (
            job
            / "output_copy_validation.log"
        ).write_text(
            "[OK] Outputs copied and SHA-256 verified: "
            + str(job / "outputs")
            + "\n"
        )

        roles = {
            "SRR20770287": "input",
            RUN: "ip",
        }

        S.dump(
            job
            / "outputs"
            / "input_provenance.json",
            {
                "schema_version": 1,
                "run_count": 2,
                "role_counts": {
                    "ip": 1,
                    "input": 1,
                },
                "library_layout": "SINGLE",
                "instrument_platform": "ILLUMINA",
                "runs": [
                    {
                        "run_accession": run,
                        "role": role,
                        "reads": 10,
                    }
                    for run, role
                    in roles.items()
                ],
            },
        )

        S.dump(
            job
            / "outputs"
            / "reference"
            / "reference_provenance.json",
            {
                "mitochondrial_accession": MT,
            },
        )

        reference = (
            job
            / "outputs"
            / "reference"
            / "genome_plus_mt.fa.fai"
        )

        reference.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        reference.write_text(
            "chr1\t10000\t0\t0\t0\n"
            + MT
            + "\t1000\t0\t0\t0\n"
        )

        S.dump(
            job
            / "outputs"
            / "alignment_parameters.json",
            {
                "schema_version": 1,
                "diagnostic_mapq": 30,
                "mitochondrial_accession": MT,
            },
        )

        (
            job
            / "outputs"
            / "alignment_qc.tsv"
        ).write_text(
            "run_accession\trole\n"
            "SRR20770287\tinput\n"
            f"{RUN}\tip\n"
        )

        S.dump(
            job
            / "outputs"
            / "verified_fastq_inputs.json",
            {
                "schema_version": 1,
                "run_count": 2,
            },
        )

        for run, role in roles.items():

            S.dump(
                job
                / "outputs"
                / run
                / "alignment_qc.json",
                {
                    "run_accession": run,
                    "role": role,
                    "duplicate_marking": "not_performed",
                    "bam_filtering": "not_performed",
                    "nonprimary_records": 0,
                    "primary_reads": 10,
                    "alignment_records": 10,
                    "diagnostic_mapq": 30,
                    "mapped_reads": 9,
                    "nuclear_mapq_ge_threshold": 5,
                },
            )

            bam = (
                job
                / "outputs"
                / run
                / "raw.sorted.bam"
            )

            bam.parent.mkdir(
                parents=True,
                exist_ok=True,
            )

            shutil.copyfile(
                "marked.bam",
                bam,
            )

            (
                job
                / "outputs"
                / run
                / "raw.sorted.bam.csi"
            ).write_bytes(
                b"synthetic-csi\n"
            )

        (
            job
            / "output.sha256"
        ).write_text(
            "".join(
                f"{S.sha(p)}  "
                f"{p.relative_to(job).as_posix()}\n"
                for p in sorted(
                    (job / "outputs").rglob("*")
                )
                if p.is_file()
            )
        )

        pointer.write_text(
            str(
                submission.resolve()
            )
            + "\n"
        )

        plan = S.prepare(
            Path.cwd()
        )

        self.assertEqual(
            len(
                plan["runs"]
            ),
            2,
        )

        self.assertEqual(
            plan["role_counts"],
            {
                "ip": 1,
                "input": 1,
            },
        )

        self.assertEqual(
            plan["library_layout"],
            "SINGLE",
        )

        self.assertEqual(
            plan["instrument_platform"],
            "ILLUMINA",
        )

        self.assertFalse(
            plan[
                "upstream_raw_bam_cleanup_authorized"
            ]
        )

        # The scientific intent of the historical regression test is
        # preserved: reports are anchored to output.sha256 and an
        # unexpected post-publication change must be rejected.
        (
            job
            / "outputs"
            / "alignment_qc.tsv"
        ).write_text(
            "changed report\n"
        )

        with self.assertRaisesRegex(
            ValueError,
            "SHA-256 mismatch",
        ):

            S.prepare(
                Path.cwd()
            )

    def test_publish_copies_verified_filtered_outputs(self):
        self.filter()
        directory = S.OUT / RUN
        directory.mkdir(parents=True)
        shutil.copyfile("filtered.bam", directory / "filtered.bam")
        pysam.index("-c", str(directory / "filtered.bam"))
        S.verify(RUN, str(directory / "filtered.bam"), "qc.json", directory / "filtered_validation.json")
        Path("inputs/reference").mkdir(parents=True)
        Path("inputs/reference/reference_provenance.json").write_text("{}\n")
        Path("saved").mkdir()
        S.publish("saved")
        hashes = S.checksums("saved/output.sha256")
        self.assertIn(f"outputs/{RUN}/filtered.bam", hashes)
        self.assertFalse(any("dupmarked" in name or "raw.sorted" in name for name in hashes))
        for name, digest in hashes.items():
            S.check(Path("saved") / name, digest)


if __name__ == "__main__":
    unittest.main()
