#!/usr/bin/env python3

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import subprocess
from collections import Counter


EXPECTED_TARGETS = {
    "H3K27ac": 4,
    "H3K4me3": 4,
    "H3K9me3": 4,
    "H3K27me3": 2,
    "H3K36me3": 2,
    "H3K4me1": 2,
}


def sha256_file(path):
    h = hashlib.sha256()

    with open(path, "rb") as fh:
        for chunk in iter(
            lambda: fh.read(1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


def read_tsv(path):
    with open(
        path,
        "r",
        encoding="utf-8",
        newline="",
    ) as fh:
        reader = csv.DictReader(
            fh,
            delimiter="\t",
        )
        return list(reader), list(reader.fieldnames or [])


def unique_index(rows, key, label):
    out = {}

    for row in rows:
        value = row[key]

        if value in out:
            raise RuntimeError(
                f"{label}: duplicate {key}={value}"
            )

        out[value] = row

    return out


def validate_peak_file(path, peak_file_type):
    expected_columns = {
        "narrowPeak": 10,
        "broadPeak": 9,
    }

    expected = expected_columns[peak_file_type]

    count = 0
    first_data = None

    with open(
        path,
        "r",
        encoding="utf-8",
    ) as fh:

        for line in fh:

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != expected:
                raise RuntimeError(
                    f"{path}: expected {expected} columns, "
                    f"found {len(fields)}"
                )

            start = int(fields[1])
            end = int(fields[2])

            if start < 0 or end <= start:
                raise RuntimeError(
                    f"{path}: invalid peak interval "
                    f"{fields[0]}:{start}-{end}"
                )

            if first_data is None:
                first_data = fields

            count += 1

    if count == 0:
        raise RuntimeError(
            f"{path}: peak file contains no peaks"
        )

    return count


def validate_signal_file(path):
    n_rows = 0
    first_row = None

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
    ) as fh:

        previous = None

        for line in fh:

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 4:
                raise RuntimeError(
                    f"{path}: expected 4-column bedGraph, "
                    f"found {len(fields)}"
                )

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            value = float(fields[3])

            if start < 0 or end <= start:
                raise RuntimeError(
                    f"{path}: invalid bedGraph interval "
                    f"{chrom}:{start}-{end}"
                )

            if not math.isfinite(value) or value < 0:
                raise RuntimeError(
                    f"{path}: invalid SPMR value {value}"
                )

            current = (chrom, start, end)

            if previous is not None:

                pchrom, pstart, pend = previous

                if chrom == pchrom:

                    if start < pstart:
                        raise RuntimeError(
                            f"{path}: unsorted bedGraph"
                        )

                    if start < pend:
                        raise RuntimeError(
                            f"{path}: overlapping bedGraph intervals"
                        )

            previous = current

            if first_row is None:
                first_row = fields

            n_rows += 1

    if n_rows == 0:
        raise RuntimeError(
            f"{path}: signal file contains no rows"
        )

    return n_rows


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--root",
        required=True,
    )

    parser.add_argument(
        "--manifest",
        required=True,
    )

    parser.add_argument(
        "--summary",
        required=True,
    )

    parser.add_argument(
        "--conditions",
        required=True,
    )

    parser.add_argument(
        "--design",
        required=True,
    )

    parser.add_argument(
        "--artifacts",
        required=True,
    )

    parser.add_argument(
        "--provenance",
        required=True,
    )

    args = parser.parse_args()

    manifest_rows, _ = read_tsv(args.manifest)
    summary_rows, _ = read_tsv(args.summary)
    condition_rows, _ = read_tsv(args.conditions)

    if len(manifest_rows) != 18:
        raise RuntimeError(
            f"Expected 18 manifest analyses; "
            f"found {len(manifest_rows)}"
        )

    if len(summary_rows) != 18:
        raise RuntimeError(
            f"Expected 18 summary analyses; "
            f"found {len(summary_rows)}"
        )

    manifest = unique_index(
        manifest_rows,
        "analysis_id",
        "manifest",
    )

    summary = unique_index(
        summary_rows,
        "analysis_id",
        "summary",
    )

    conditions = unique_index(
        condition_rows,
        "run_accession",
        "conditions",
    )

    if set(manifest) != set(summary):
        raise RuntimeError(
            "Manifest and production-summary "
            "analysis sets differ"
        )

    target_counts = Counter(
        row["target"]
        for row in manifest_rows
    )

    if dict(target_counts) != EXPECTED_TARGETS:
        raise RuntimeError(
            f"Unexpected target counts: "
            f"{dict(target_counts)}"
        )

    peak_modes = Counter(
        row["peak_mode"]
        for row in manifest_rows
    )

    if dict(peak_modes) != {
        "narrow": 8,
        "broad": 10,
    }:
        raise RuntimeError(
            f"Unexpected peak modes: {dict(peak_modes)}"
        )

    design_rows = []
    artifact_rows = []

    for m in manifest_rows:

        analysis_id = m["analysis_id"]
        s = summary[analysis_id]

        if m["final_status"] != "PASS":
            raise RuntimeError(
                f"{analysis_id}: final_status != PASS"
            )

        if m["reference_accession"] != "GCF_003668045.3":
            raise RuntimeError(
                f"{analysis_id}: unexpected reference"
            )

        if m["ip_run"] not in conditions:
            raise RuntimeError(
                f"{analysis_id}: IP condition metadata missing"
            )

        if m["control_run"] not in conditions:
            raise RuntimeError(
                f"{analysis_id}: control condition metadata missing"
            )

        ip = conditions[m["ip_run"]]
        ctrl = conditions[m["control_run"]]

        if ip["library_role"] != "ip":
            raise RuntimeError(
                f"{analysis_id}: treatment is not labelled IP"
            )

        if ctrl["library_role"] != "input":
            raise RuntimeError(
                f"{analysis_id}: control is not labelled input"
            )

        if ip["study_accession"] != m["study"]:
            raise RuntimeError(
                f"{analysis_id}: IP study mismatch"
            )

        if ctrl["study_accession"] != m["study"]:
            raise RuntimeError(
                f"{analysis_id}: control study mismatch"
            )

        if (
            ip["time_value"] != ctrl["time_value"]
            or ip["time_unit"] != ctrl["time_unit"]
        ):
            raise RuntimeError(
                f"{analysis_id}: IP/Input condition mismatch"
            )

        if (
            ip["declared_target"]
            and ip["declared_target"] != m["target"]
        ):
            raise RuntimeError(
                f"{analysis_id}: target mismatch"
            )

        fields_to_match = [
            "study",
            "ip_run",
            "control_run",
            "target",
            "peak_mode",
            "fragment_size_bp",
            "frip",
            "input_overlap_fraction",
            "manifest_sha256",
            "production_commit",
        ]

        for field in fields_to_match:

            if m[field] != s[field]:
                raise RuntimeError(
                    f"{analysis_id}: manifest/summary "
                    f"mismatch for {field}"
                )

        peak_path = os.path.join(
            args.root,
            m["peak_file_relative_path"],
        )

        signal_path = os.path.join(
            args.root,
            m["treat_signal_relative_path"],
        )

        for path in (peak_path, signal_path):

            if not os.path.isfile(path):
                raise RuntimeError(
                    f"{analysis_id}: missing artifact {path}"
                )

            if os.path.getsize(path) == 0:
                raise RuntimeError(
                    f"{analysis_id}: empty artifact {path}"
                )

        peak_count = validate_peak_file(
            peak_path,
            m["peak_file_type"],
        )

        expected_peak_count = int(
            s["peak_count"]
        )

        if peak_count != expected_peak_count:
            raise RuntimeError(
                f"{analysis_id}: peak count "
                f"{peak_count} != {expected_peak_count}"
            )

        signal_rows = validate_signal_file(
            signal_path,
        )

        condition_key = (
            f"{m['study']}:"
            f"{ip['time_value']}_"
            f"{ip['time_unit']}"
        )

        peak_size = os.path.getsize(
            peak_path
        )

        signal_size = os.path.getsize(
            signal_path
        )

        peak_sha = sha256_file(
            peak_path
        )

        signal_sha = sha256_file(
            signal_path
        )

        artifact_rows.extend([
            {
                "analysis_id":
                    analysis_id,
                "artifact_type":
                    "peak",
                "file_type":
                    m["peak_file_type"],
                "relative_path":
                    m["peak_file_relative_path"],
                "size_bytes":
                    str(peak_size),
                "sha256":
                    peak_sha,
                "record_count":
                    str(peak_count),
            },
            {
                "analysis_id":
                    analysis_id,
                "artifact_type":
                    "treat_pileup_spmr",
                "file_type":
                    "bedGraph.gz",
                "relative_path":
                    m["treat_signal_relative_path"],
                "size_bytes":
                    str(signal_size),
                "sha256":
                    signal_sha,
                "record_count":
                    str(signal_rows),
            },
        ])

        design_rows.append({
            "analysis_id":
                analysis_id,
            "study":
                m["study"],
            "target":
                m["target"],
            "peak_mode":
                m["peak_mode"],
            "ip_run":
                m["ip_run"],
            "control_run":
                m["control_run"],
            "condition_key":
                condition_key,
            "time_value":
                ip["time_value"],
            "time_unit":
                ip["time_unit"],
            "fragment_size_policy":
                s["fragment_size_policy"],
            "fragment_size_bp":
                m["fragment_size_bp"],
            "reference_assembly":
                m["reference_assembly"],
            "reference_accession":
                m["reference_accession"],
            "peak_file_type":
                m["peak_file_type"],
            "peak_file_relative_path":
                m["peak_file_relative_path"],
            "treat_signal_relative_path":
                m["treat_signal_relative_path"],
            "peak_count":
                s["peak_count"],
            "frip":
                m["frip"],
            "input_overlap_fraction":
                m["input_overlap_fraction"],
            "ip_input_overlap_ratio":
                s["ip_input_overlap_ratio"],
            "median_fold_enrichment":
                s["median_fold_enrichment"],
            "median_minus_log10_q":
                s["median_minus_log10_q"],
            "median_peak_width_bp":
                s["median_peak_width_bp"],
            "mean_peak_width_bp":
                s["mean_peak_width_bp"],
            "union_peak_span_bp":
                s["union_peak_span_bp"],
            "nuclear_coverage_fraction":
                s["nuclear_coverage_fraction"],
            "peak_qvalue":
                m["peak_qvalue"],
            "broad_cutoff":
                m["broad_cutoff"],
            "keep_dup":
                m["keep_dup"],
            "blacklist_status":
                m["blacklist_status"],
            "replicate_concordance_status":
                m["replicate_concordance_status"],
            "manifest_sha256":
                m["manifest_sha256"],
            "production_commit":
                m["production_commit"],
            "final_status":
                m["final_status"],
        })

    # --------------------------------------------------------
    # Validate biological structure
    # --------------------------------------------------------

    condition_targets = Counter(
        (
            row["study"],
            row["condition_key"],
            row["target"],
        )
        for row in design_rows
    )

    if any(
        n != 1
        for n in condition_targets.values()
    ):
        raise RuntimeError(
            "Study-condition-target is not unique"
        )

    if len(artifact_rows) != 36:
        raise RuntimeError(
            f"Expected 36 integration artifacts; "
            f"found {len(artifact_rows)}"
        )

    # --------------------------------------------------------
    # Write design
    # --------------------------------------------------------

    design_fields = [
        "analysis_id",
        "study",
        "target",
        "peak_mode",
        "ip_run",
        "control_run",
        "condition_key",
        "time_value",
        "time_unit",
        "fragment_size_policy",
        "fragment_size_bp",
        "reference_assembly",
        "reference_accession",
        "peak_file_type",
        "peak_file_relative_path",
        "treat_signal_relative_path",
        "peak_count",
        "frip",
        "input_overlap_fraction",
        "ip_input_overlap_ratio",
        "median_fold_enrichment",
        "median_minus_log10_q",
        "median_peak_width_bp",
        "mean_peak_width_bp",
        "union_peak_span_bp",
        "nuclear_coverage_fraction",
        "peak_qvalue",
        "broad_cutoff",
        "keep_dup",
        "blacklist_status",
        "replicate_concordance_status",
        "manifest_sha256",
        "production_commit",
        "final_status",
    ]

    with open(
        args.design,
        "w",
        encoding="utf-8",
        newline="",
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=design_fields,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(design_rows)

    # --------------------------------------------------------
    # Write artifact inventory
    # --------------------------------------------------------

    artifact_fields = [
        "analysis_id",
        "artifact_type",
        "file_type",
        "relative_path",
        "size_bytes",
        "sha256",
        "record_count",
    ]

    with open(
        args.artifacts,
        "w",
        encoding="utf-8",
        newline="",
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=artifact_fields,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(artifact_rows)

    # --------------------------------------------------------
    # Provenance
    # --------------------------------------------------------

    try:
        pipeline_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
        ).strip()

    except Exception:
        pipeline_commit = "unavailable"

    production_commits = sorted({
        row["production_commit"]
        for row in design_rows
    })

    provenance = {
        "product":
            "chipseq_integration_input_v1",

        "description":
            (
                "Validated condition-specific ChIP-seq "
                "integration input for 18 IP-Input analyses."
            ),

        "pipeline_commit_at_build":
            pipeline_commit,

        "source_root":
            args.root,

        "reference_accession":
            "GCF_003668045.3",

        "analysis_count":
            len(design_rows),

        "artifact_count":
            len(artifact_rows),

        "targets":
            EXPECTED_TARGETS,

        "production_commits":
            production_commits,

        "integration_policy": {
            "unit":
                (
                    "one condition-specific IP-Input "
                    "analysis"
                ),

            "peak_evidence":
                (
                    "MACS3 narrowPeak or broadPeak "
                    "according to validated mark policy"
                ),

            "continuous_signal":
                (
                    "MACS3 treatment pileup generated "
                    "with SPMR normalization"
                ),

            "aggregation":
                (
                    "No cross-condition or cross-study "
                    "aggregation at this stage"
                ),

            "control_lambda":
                (
                    "Not required for integration; "
                    "removed after completed QC and "
                    "production verification"
                ),
        },

        "sources": {
            "integration_manifest": {
                "path":
                    args.manifest,
                "sha256":
                    sha256_file(args.manifest),
            },

            "production_summary": {
                "path":
                    args.summary,
                "sha256":
                    sha256_file(args.summary),
            },

            "condition_metadata": {
                "path":
                    args.conditions,
                "sha256":
                    sha256_file(args.conditions),
            },

            "script": {
                "path":
                    os.path.abspath(__file__),
                "sha256":
                    sha256_file(__file__),
            },
        },
    }

    with open(
        args.provenance,
        "w",
        encoding="utf-8",
    ) as fh:

        json.dump(
            provenance,
            fh,
            indent=2,
            sort_keys=True,
        )

        fh.write("\n")

    print(f"ANALYSES={len(design_rows)}")
    print(f"ARTIFACTS={len(artifact_rows)}")
    print(
        f"PEAK_ARTIFACTS="
        f"{sum(r['artifact_type'] == 'peak' for r in artifact_rows)}"
    )
    print(
        f"SPMR_ARTIFACTS="
        f"{sum(r['artifact_type'] == 'treat_pileup_spmr' for r in artifact_rows)}"
    )

    print()
    print(
        "PASS: ChIP integration input frozen "
        "at analysis and artifact level"
    )


if __name__ == "__main__":
    main()
