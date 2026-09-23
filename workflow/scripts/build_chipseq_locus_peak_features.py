#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import os
import subprocess
from collections import Counter, defaultdict


WINDOWS = (1000, 10000, 50000)

EXPECTED_SEMANTIC = {
    "point_coordinate": 22,
    "interbase_cut_boundary": 4,
    "paired_nick_interval": 3,
    "mapped_interval": 8,
}

EXPECTED_ROLES = {
    "positive": 27,
    "negative": 4,
    "support_only": 6,
}

EXPECTED_TARGET_ANALYSES = {
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
        return (
            list(reader),
            list(reader.fieldnames or []),
        )


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


def read_peaks(path, expected_type):
    expected_columns = {
        "narrowPeak": 10,
        "broadPeak": 9,
    }[expected_type]

    peaks = defaultdict(list)
    count = 0

    with open(
        path,
        "r",
        encoding="utf-8",
    ) as fh:

        for line in fh:

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != expected_columns:
                raise RuntimeError(
                    f"{path}: expected {expected_columns} "
                    f"columns; found {len(fields)}"
                )

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            if start < 0 or end <= start:
                raise RuntimeError(
                    f"{path}: invalid interval "
                    f"{chrom}:{start}-{end}"
                )

            peaks[chrom].append(
                (start, end)
            )

            count += 1

    for chrom in peaks:
        peaks[chrom].sort()

    return peaks, count


def overlap(a_start, a_end, b_start, b_end):
    return (
        a_start < b_end
        and a_end > b_start
    )


def overlapping_peaks(intervals, qstart, qend):
    return [
        (start, end)
        for start, end in intervals
        if overlap(
            start,
            end,
            qstart,
            qend,
        )
    ]


def clipped_union_bp(intervals, qstart, qend):

    clipped = []

    for start, end in intervals:

        left = max(start, qstart)
        right = min(end, qend)

        if left < right:
            clipped.append(
                (left, right)
            )

    if not clipped:
        return 0

    clipped.sort()

    total = 0
    current_start, current_end = clipped[0]

    for start, end in clipped[1:]:

        if start <= current_end:
            current_end = max(
                current_end,
                end,
            )

        else:
            total += (
                current_end
                - current_start
            )

            current_start = start
            current_end = end

    total += (
        current_end
        - current_start
    )

    return total


def distance_anchor_to_peak(anchor0, start, end):

    if start <= anchor0 < end:
        return 0

    if anchor0 < start:
        return start - anchor0

    return anchor0 - (end - 1)


def nearest_anchor_distance(anchor0, intervals):

    if not intervals:
        return None

    return min(
        distance_anchor_to_peak(
            anchor0,
            start,
            end,
        )
        for start, end in intervals
    )


def distance_interval_to_peak(
    qstart,
    qend,
    pstart,
    pend,
):

    if overlap(
        qstart,
        qend,
        pstart,
        pend,
    ):
        return 0

    if qend <= pstart:
        return (
            pstart
            - (qend - 1)
        )

    return (
        qstart
        - (pend - 1)
    )


def nearest_interval_distance(
    qstart,
    qend,
    intervals,
):

    if not intervals:
        return None

    return min(
        distance_interval_to_peak(
            qstart,
            qend,
            start,
            end,
        )
        for start, end in intervals
    )


def fmt_distance(value):
    if value is None:
        return "NA"

    return str(value)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--root",
        required=True,
    )

    parser.add_argument(
        "--benchmark",
        required=True,
    )

    parser.add_argument(
        "--locus-audit",
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
        "--output",
        required=True,
    )

    parser.add_argument(
        "--qc",
        required=True,
    )

    parser.add_argument(
        "--provenance",
        required=True,
    )

    args = parser.parse_args()

    benchmark_rows, _ = read_tsv(
        args.benchmark
    )

    audit_rows, _ = read_tsv(
        args.locus_audit
    )

    design_rows, _ = read_tsv(
        args.design
    )

    artifact_rows, _ = read_tsv(
        args.artifacts
    )

    if len(benchmark_rows) != 37:
        raise RuntimeError(
            f"Expected 37 benchmark loci; "
            f"found {len(benchmark_rows)}"
        )

    if len(audit_rows) != 37:
        raise RuntimeError(
            f"Expected 37 locus-audit rows; "
            f"found {len(audit_rows)}"
        )

    if len(design_rows) != 18:
        raise RuntimeError(
            f"Expected 18 ChIP analyses; "
            f"found {len(design_rows)}"
        )

    benchmark = unique_index(
        benchmark_rows,
        "canonical_locus_id",
        "benchmark",
    )

    audit = unique_index(
        audit_rows,
        "canonical_locus_id",
        "locus_audit",
    )

    design = unique_index(
        design_rows,
        "analysis_id",
        "design",
    )

    if set(benchmark) != set(audit):
        raise RuntimeError(
            "Benchmark and locus-audit sets differ"
        )

    roles = Counter(
        row["benchmark_role"]
        for row in benchmark_rows
    )

    if dict(roles) != EXPECTED_ROLES:
        raise RuntimeError(
            f"Unexpected benchmark roles: "
            f"{dict(roles)}"
        )

    semantic = Counter(
        row["semantic_class"]
        for row in audit_rows
    )

    if dict(semantic) != EXPECTED_SEMANTIC:
        raise RuntimeError(
            f"Unexpected semantic classes: "
            f"{dict(semantic)}"
        )

    targets = Counter(
        row["target"]
        for row in design_rows
    )

    if dict(targets) != EXPECTED_TARGET_ANALYSES:
        raise RuntimeError(
            f"Unexpected target distribution: "
            f"{dict(targets)}"
        )

    peak_artifacts = {
        row["analysis_id"]: row
        for row in artifact_rows
        if row["artifact_type"] == "peak"
    }

    if len(peak_artifacts) != 18:
        raise RuntimeError(
            f"Expected 18 peak artifacts; "
            f"found {len(peak_artifacts)}"
        )

    if set(peak_artifacts) != set(design):
        raise RuntimeError(
            "Design and peak-artifact analysis sets differ"
        )

    output_rows = []

    # --------------------------------------------------------
    # One analysis at a time
    # --------------------------------------------------------

    for analysis in design_rows:

        analysis_id = analysis[
            "analysis_id"
        ]

        artifact = peak_artifacts[
            analysis_id
        ]

        if (
            artifact["relative_path"]
            != analysis[
                "peak_file_relative_path"
            ]
        ):
            raise RuntimeError(
                f"{analysis_id}: peak path mismatch"
            )

        if (
            artifact["file_type"]
            != analysis["peak_file_type"]
        ):
            raise RuntimeError(
                f"{analysis_id}: peak type mismatch"
            )

        peak_path = os.path.join(
            args.root,
            artifact["relative_path"],
        )

        if not os.path.isfile(peak_path):
            raise RuntimeError(
                f"{analysis_id}: missing peak file "
                f"{peak_path}"
            )

        observed_sha = sha256_file(
            peak_path
        )

        if observed_sha != artifact["sha256"]:
            raise RuntimeError(
                f"{analysis_id}: peak SHA256 changed"
            )

        peaks, peak_count = read_peaks(
            peak_path,
            artifact["file_type"],
        )

        if (
            peak_count
            != int(artifact["record_count"])
        ):
            raise RuntimeError(
                f"{analysis_id}: peak artifact "
                f"record count mismatch"
            )

        if (
            peak_count
            != int(analysis["peak_count"])
        ):
            raise RuntimeError(
                f"{analysis_id}: peak design "
                f"record count mismatch"
            )

        # ----------------------------------------------------
        # 37 canonical loci
        # ----------------------------------------------------

        for b in benchmark_rows:

            locus_id = b[
                "canonical_locus_id"
            ]

            a = audit[locus_id]

            if (
                a["locus_name"]
                != b["locus_name"]
            ):
                raise RuntimeError(
                    f"{locus_id}: locus name mismatch"
                )

            if (
                a["seqname"]
                != b["target_seqname"]
            ):
                raise RuntimeError(
                    f"{locus_id}: seqname mismatch"
                )

            if (
                a["start_1based"]
                != b["target_start"]
            ):
                raise RuntimeError(
                    f"{locus_id}: start mismatch"
                )

            if (
                a["end_1based"]
                != b["target_end"]
            ):
                raise RuntimeError(
                    f"{locus_id}: end mismatch"
                )

            start1 = int(
                b["target_start"]
            )

            end1 = int(
                b["target_end"]
            )

            width = (
                end1
                - start1
                + 1
            )

            if width != int(a["width_bp"]):
                raise RuntimeError(
                    f"{locus_id}: width mismatch"
                )

            # Canonical benchmark:
            # 1-based closed -> 0-based half-open.
            start0 = start1 - 1
            end0 = end1

            # Representative anchor:
            # floor((start + end) / 2), 1-based.
            anchor1 = (
                start1
                + end1
            ) // 2

            anchor0 = (
                anchor1
                - 1
            )

            intervals = peaks.get(
                b["target_seqname"],
                [],
            )

            exact_hits = overlapping_peaks(
                intervals,
                start0,
                end0,
            )

            exact_bp = clipped_union_bp(
                exact_hits,
                start0,
                end0,
            )

            exact_n = len(
                exact_hits
            )

            exact_overlap = int(
                exact_n > 0
            )

            exact_fraction = (
                exact_bp
                / width
            )

            anchor_in_peak = int(
                any(
                    start
                    <= anchor0
                    < end
                    for start, end in intervals
                )
            )

            nearest_anchor = (
                nearest_anchor_distance(
                    anchor0,
                    intervals,
                )
            )

            nearest_locus = (
                nearest_interval_distance(
                    start0,
                    end0,
                    intervals,
                )
            )

            row = {
                "canonical_locus_id":
                    locus_id,

                "locus_name":
                    b["locus_name"],

                "benchmark_role":
                    b["benchmark_role"],

                "best_evidence_tier":
                    b["best_evidence_tier"],

                "n_observations":
                    b["n_observations"],

                "target_assembly":
                    b["target_assembly"],

                "target_seqname":
                    b["target_seqname"],

                "target_start_1based":
                    str(start1),

                "target_end_1based":
                    str(end1),

                "locus_width_bp":
                    str(width),

                "locus_semantic_class":
                    a["semantic_class"],

                "mapping_confidence":
                    b["mapping_confidence"],

                "supporting_studies":
                    b["supporting_studies"],

                "chip_anchor_1based":
                    str(anchor1),

                "analysis_id":
                    analysis_id,

                "study":
                    analysis["study"],

                "condition_key":
                    analysis["condition_key"],

                "target":
                    analysis["target"],

                "peak_mode":
                    analysis["peak_mode"],

                "ip_run":
                    analysis["ip_run"],

                "control_run":
                    analysis["control_run"],

                "peak_file_type":
                    analysis["peak_file_type"],

                "analysis_peak_count":
                    analysis["peak_count"],

                "analysis_frip":
                    analysis["frip"],

                "analysis_input_overlap_fraction":
                    analysis[
                        "input_overlap_fraction"
                    ],

                "analysis_ip_input_overlap_ratio":
                    analysis[
                        "ip_input_overlap_ratio"
                    ],

                "analysis_median_fold_enrichment":
                    analysis[
                        "median_fold_enrichment"
                    ],

                "chip_exact_overlap":
                    str(exact_overlap),

                "chip_exact_n_peaks":
                    str(exact_n),

                "chip_exact_overlap_bp":
                    str(exact_bp),

                "chip_exact_overlap_fraction_of_locus":
                    repr(exact_fraction),

                "chip_anchor_in_peak":
                    str(anchor_in_peak),

                "chip_nearest_peak_distance_to_anchor_bp":
                    fmt_distance(
                        nearest_anchor
                    ),

                "chip_nearest_peak_distance_to_locus_bp":
                    fmt_distance(
                        nearest_locus
                    ),
            }

            for half_window in WINDOWS:

                wstart = max(
                    0,
                    anchor0
                    - half_window,
                )

                # Includes the anchor base and
                # half_window bases on each side.
                wend = (
                    anchor0
                    + half_window
                    + 1
                )

                window_bp = (
                    wend
                    - wstart
                )

                window_hits = (
                    overlapping_peaks(
                        intervals,
                        wstart,
                        wend,
                    )
                )

                enriched_bp = (
                    clipped_union_bp(
                        window_hits,
                        wstart,
                        wend,
                    )
                )

                n_peaks = len(
                    window_hits
                )

                row[
                    f"chip_n_peaks_pm{half_window}bp"
                ] = str(
                    n_peaks
                )

                row[
                    f"chip_any_peak_pm{half_window}bp"
                ] = str(
                    int(n_peaks > 0)
                )

                row[
                    f"chip_enriched_bp_pm{half_window}bp"
                ] = str(
                    enriched_bp
                )

                row[
                    f"chip_enriched_fraction_pm{half_window}bp"
                ] = repr(
                    enriched_bp
                    / window_bp
                )

            output_rows.append(
                row
            )

    # --------------------------------------------------------
    # Final structural checks
    # --------------------------------------------------------

    if len(output_rows) != 37 * 18:
        raise RuntimeError(
            f"Expected 666 rows; "
            f"found {len(output_rows)}"
        )

    pair_keys = {
        (
            row["canonical_locus_id"],
            row["analysis_id"],
        )
        for row in output_rows
    }

    if len(pair_keys) != 666:
        raise RuntimeError(
            "Duplicate locus-analysis pair"
        )

    locus_counts = Counter(
        row["canonical_locus_id"]
        for row in output_rows
    )

    if any(
        n != 18
        for n in locus_counts.values()
    ):
        raise RuntimeError(
            "Not every locus has 18 analyses"
        )

    analysis_counts = Counter(
        row["analysis_id"]
        for row in output_rows
    )

    if any(
        n != 37
        for n in analysis_counts.values()
    ):
        raise RuntimeError(
            "Not every analysis has 37 loci"
        )

    # --------------------------------------------------------
    # Write output
    # --------------------------------------------------------

    fieldnames = list(
        output_rows[0].keys()
    )

    with open(
        args.output,
        "w",
        encoding="utf-8",
        newline="",
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(
            output_rows
        )

    # --------------------------------------------------------
    # QC
    # --------------------------------------------------------

    exact_pairs = sum(
        int(row["chip_exact_overlap"])
        for row in output_rows
    )

    anchor_pairs = sum(
        int(row["chip_anchor_in_peak"])
        for row in output_rows
    )

    qc_rows = [
        (
            "locus_analysis_rows",
            len(output_rows),
        ),
        (
            "benchmark_loci",
            len(locus_counts),
        ),
        (
            "chip_analyses",
            len(analysis_counts),
        ),
        (
            "exact_overlap_pairs",
            exact_pairs,
        ),
        (
            "anchor_in_peak_pairs",
            anchor_pairs,
        ),
    ]

    for role, n in sorted(
        Counter(
            row["benchmark_role"]
            for row in output_rows
        ).items()
    ):
        qc_rows.append(
            (
                f"rows_role_{role}",
                n,
            )
        )

    for semantic_class, n in sorted(
        Counter(
            row["locus_semantic_class"]
            for row in output_rows
        ).items()
    ):
        qc_rows.append(
            (
                f"rows_semantic_{semantic_class}",
                n,
            )
        )

    for target, n in sorted(
        Counter(
            row["target"]
            for row in output_rows
        ).items()
    ):
        qc_rows.append(
            (
                f"rows_target_{target}",
                n,
            )
        )

        loci_any = len({
            row["canonical_locus_id"]
            for row in output_rows
            if (
                row["target"] == target
                and row["chip_exact_overlap"] == "1"
            )
        })

        qc_rows.append(
            (
                f"loci_with_exact_overlap_{target}",
                loci_any,
            )
        )

    for half_window in WINDOWS:

        pairs_any = sum(
            int(
                row[
                    f"chip_any_peak_pm{half_window}bp"
                ]
            )
            for row in output_rows
        )

        qc_rows.append(
            (
                f"pairs_any_peak_pm{half_window}bp",
                pairs_any,
            )
        )

    with open(
        args.qc,
        "w",
        encoding="utf-8",
        newline="",
    ) as fh:

        writer = csv.writer(
            fh,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow([
            "metric",
            "value",
        ])

        writer.writerows(
            qc_rows
        )

    # --------------------------------------------------------
    # Provenance
    # --------------------------------------------------------

    try:
        head = subprocess.check_output(
            [
                "git",
                "rev-parse",
                "HEAD",
            ],
            text=True,
        ).strip()

    except Exception:
        head = "unavailable"

    provenance = {
        "product":
            "chipseq_locus_analysis_peak_features_v1",

        "description":
            (
                "Condition-specific ChIP-seq peak "
                "features for all 37 canonical CHO "
                "benchmark loci across 18 validated "
                "IP-Input analyses."
            ),

        "repository_HEAD_at_generation":
            head,

        "dimensions": {
            "benchmark_loci":
                37,
            "chip_analyses":
                18,
            "rows":
                666,
        },

        "coordinate_policy": {
            "benchmark":
                (
                    "Canonical loci are 1-based closed "
                    "intervals."
                ),

            "peak_files":
                (
                    "MACS3 narrowPeak/broadPeak are "
                    "0-based half-open intervals."
                ),

            "conversion":
                (
                    "benchmark_start0 = target_start_1based - 1; "
                    "benchmark_end0 = target_end_1based"
                ),

            "anchor":
                (
                    "floor((start_1based + end_1based)/2); "
                    "representative context point only, "
                    "not assumed to be an exact insertion "
                    "coordinate for interval-class loci."
                ),

            "full_locus_priority":
                (
                    "Exact peak overlap is evaluated on "
                    "the full canonical interval and is "
                    "the primary locus-local feature."
                ),

            "windows":
                (
                    "Representative-anchor windows include "
                    "the anchor base plus W bases on each "
                    "side; nominal widths are 2001, 20001 "
                    "and 100001 bp for W=1000,10000,50000."
                ),
        },

        "peak_policy": {
            "narrow_marks":
                [
                    "H3K27ac",
                    "H3K4me3",
                ],

            "broad_marks":
                [
                    "H3K4me1",
                    "H3K27me3",
                    "H3K36me3",
                    "H3K9me3",
                ],

            "aggregation":
                (
                    "No aggregation across conditions, "
                    "studies or marks at this stage."
                ),

            "interpretation":
                (
                    "Descriptive regulatory-context "
                    "features; no causal or predictive "
                    "claim."
                ),
        },

        "sources": {
            "benchmark": {
                "path":
                    args.benchmark,
                "sha256":
                    sha256_file(
                        args.benchmark
                    ),
            },

            "locus_audit": {
                "path":
                    args.locus_audit,
                "sha256":
                    sha256_file(
                        args.locus_audit
                    ),
            },

            "analysis_design": {
                "path":
                    args.design,
                "sha256":
                    sha256_file(
                        args.design
                    ),
            },

            "artifact_inventory": {
                "path":
                    args.artifacts,
                "sha256":
                    sha256_file(
                        args.artifacts
                    ),
            },

            "script": {
                "path":
                    os.path.abspath(
                        __file__
                    ),
                "sha256":
                    sha256_file(
                        __file__
                    ),
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

    print(
        f"LOCUS_ANALYSIS_ROWS="
        f"{len(output_rows)}"
    )

    print(
        f"OUTPUT_COLUMNS="
        f"{len(fieldnames)}"
    )

    print(
        f"EXACT_OVERLAP_PAIRS="
        f"{exact_pairs}"
    )

    print(
        f"ANCHOR_IN_PEAK_PAIRS="
        f"{anchor_pairs}"
    )

    print()
    print(
        "PASS: condition-specific ChIP peak "
        "features built for 37 x 18 design"
    )


if __name__ == "__main__":
    main()
