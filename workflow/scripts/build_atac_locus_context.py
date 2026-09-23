#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from collections import Counter


WINDOWS = (1000, 10000, 50000)


def sha256_file(path):
    h = hashlib.sha256()

    with open(path, "rb") as fh:
        for chunk in iter(
            lambda: fh.read(1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


def load_unique(path, key):

    result = {}

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

        for row in reader:

            value = row[key]

            if value in result:
                raise RuntimeError(
                    f"{path}: duplicate {key}={value}"
                )

            result[value] = row

    return result


def as_int(row, field):
    return int(row[field])


def as_float(row, field):
    value = float(row[field])

    if not math.isfinite(value):
        raise RuntimeError(
            f"Non-finite value in {field}: {row[field]}"
        )

    return value


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--benchmark",
        required=True,
    )

    parser.add_argument(
        "--locus-audit",
        required=True,
    )

    parser.add_argument(
        "--source",
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

    benchmark = load_unique(
        args.benchmark,
        "canonical_locus_id",
    )

    locus_audit = load_unique(
        args.locus_audit,
        "canonical_locus_id",
    )

    source = load_unique(
        args.source,
        "canonical_locus_id",
    )

    # --------------------------------------------------------
    # 1. Basic structural validation
    # --------------------------------------------------------

    if len(benchmark) != 37:
        raise RuntimeError(
            f"Expected 37 benchmark loci; "
            f"found {len(benchmark)}"
        )

    if len(locus_audit) != 37:
        raise RuntimeError(
            f"Expected 37 locus-audit rows; "
            f"found {len(locus_audit)}"
        )

    if len(source) != 37:
        raise RuntimeError(
            f"Expected 37 ATAC loci; "
            f"found {len(source)}"
        )

    if set(benchmark) != set(locus_audit):
        raise RuntimeError(
            "Benchmark and locus-audit locus sets differ"
        )

    if set(benchmark) != set(source):
        raise RuntimeError(
            "Benchmark and ATAC locus sets differ"
        )

    roles = Counter(
        row["benchmark_role"]
        for row in benchmark.values()
    )

    expected_roles = {
        "positive": 27,
        "negative": 4,
        "support_only": 6,
    }

    if dict(roles) != expected_roles:
        raise RuntimeError(
            f"Unexpected benchmark roles: {dict(roles)}"
        )

    # --------------------------------------------------------
    # 2. Output schema
    # --------------------------------------------------------

    columns = [
        # Benchmark
        "canonical_locus_id",
        "locus_name",
        "benchmark_role",
        "best_evidence_tier",
        "n_observations",
        "target_assembly",
        "target_seqname",
        "target_start_1based",
        "target_end_1based",
        "locus_width_bp",
        "locus_semantic_class",
        "mapping_confidence",
        "supporting_studies",

        # Representative anchor
        "atac_anchor_1based",

        # Primary IDR <= 0.10 — complete locus interval
        "atac_idr10_exact_overlap",
        "atac_idr10_exact_n_unique_peaks",
        "atac_idr10_exact_overlap_bp",
        "atac_idr10_exact_overlap_fraction_of_locus",
        "atac_idr10_anchor_in_peak",
        "atac_idr10_nearest_peak_distance_bp",

        # IDR10 local context + continuous signal
        "atac_idr10_n_unique_peaks_pm1000bp",
        "atac_idr10_any_peak_pm1000bp",
        "atac_idr10_accessible_bp_pm1000bp",
        "atac_spmr_consensus_pm1000bp_mean",

        "atac_idr10_n_unique_peaks_pm10000bp",
        "atac_idr10_any_peak_pm10000bp",
        "atac_idr10_accessible_bp_pm10000bp",
        "atac_spmr_consensus_pm10000bp_mean",

        "atac_idr10_n_unique_peaks_pm50000bp",
        "atac_idr10_any_peak_pm50000bp",
        "atac_idr10_accessible_bp_pm50000bp",
        "atac_spmr_consensus_pm50000bp_mean",

        # Sensitivity IDR <= 0.05 — complete locus interval
        "atac_idr05_exact_overlap",
        "atac_idr05_exact_n_unique_peaks",
        "atac_idr05_exact_overlap_bp",
        "atac_idr05_exact_overlap_fraction_of_locus",
        "atac_idr05_anchor_in_peak",
        "atac_idr05_nearest_peak_distance_bp",

        # IDR05 local context
        "atac_idr05_n_unique_peaks_pm1000bp",
        "atac_idr05_any_peak_pm1000bp",
        "atac_idr05_accessible_bp_pm1000bp",

        "atac_idr05_n_unique_peaks_pm10000bp",
        "atac_idr05_any_peak_pm10000bp",
        "atac_idr05_accessible_bp_pm10000bp",

        "atac_idr05_n_unique_peaks_pm50000bp",
        "atac_idr05_any_peak_pm50000bp",
        "atac_idr05_accessible_bp_pm50000bp",
    ]

    final_rows = []

    # --------------------------------------------------------
    # 3. Preserve canonical benchmark order
    # --------------------------------------------------------

    with open(
        args.benchmark,
        "r",
        encoding="utf-8",
        newline="",
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t",
        )

        for b in reader:

            locus_id = b["canonical_locus_id"]
            l = locus_audit[locus_id]
            s = source[locus_id]

            # ------------------------------------------------
            # Cross-source identity
            # ------------------------------------------------

            exact_checks = [
                (
                    "locus_name",
                    s["locus_name"],
                    b["locus_name"],
                ),
                (
                    "benchmark_role",
                    s["benchmark_role"],
                    b["benchmark_role"],
                ),
                (
                    "target_assembly",
                    s["target_assembly"],
                    b["target_assembly"],
                ),
                (
                    "target_seqname",
                    s["target_seqname"],
                    b["target_seqname"],
                ),
                (
                    "target_start",
                    s["target_start_1based"],
                    b["target_start"],
                ),
                (
                    "target_end",
                    s["target_end_1based"],
                    b["target_end"],
                ),
            ]

            for name, observed, expected in exact_checks:

                if observed != expected:
                    raise RuntimeError(
                        f"{locus_id}: {name} mismatch: "
                        f"{observed} != {expected}"
                    )

            audit_checks = [
                (
                    "audit_locus_name",
                    l["locus_name"],
                    b["locus_name"],
                ),
                (
                    "audit_seqname",
                    l["seqname"],
                    b["target_seqname"],
                ),
                (
                    "audit_start",
                    l["start_1based"],
                    b["target_start"],
                ),
                (
                    "audit_end",
                    l["end_1based"],
                    b["target_end"],
                ),
            ]

            for name, observed, expected in audit_checks:

                if observed != expected:
                    raise RuntimeError(
                        f"{locus_id}: {name} mismatch: "
                        f"{observed} != {expected}"
                    )

            semantic_class = l["semantic_class"]

            allowed_semantic_classes = {
                "point_coordinate",
                "interbase_cut_boundary",
                "paired_nick_interval",
                "mapped_interval",
            }

            if semantic_class not in allowed_semantic_classes:
                raise RuntimeError(
                    f"{locus_id}: invalid semantic_class="
                    f"{semantic_class}"
                )

            start = int(b["target_start"])
            end = int(b["target_end"])
            width = end - start + 1

            if as_int(
                s,
                "locus_length_bp",
            ) != width:
                raise RuntimeError(
                    f"{locus_id}: locus width mismatch"
                )

            # ------------------------------------------------
            # IDR subset / monotonicity validation
            # ------------------------------------------------

            exact10 = as_int(
                s,
                "atac_idr10_exact_overlap",
            )

            exact05 = as_int(
                s,
                "atac_idr05_exact_overlap",
            )

            n_exact10 = as_int(
                s,
                "atac_idr10_exact_n_unique_peaks",
            )

            n_exact05 = as_int(
                s,
                "atac_idr05_exact_n_unique_peaks",
            )

            bp_exact10 = as_int(
                s,
                "atac_idr10_exact_overlap_bp",
            )

            bp_exact05 = as_int(
                s,
                "atac_idr05_exact_overlap_bp",
            )

            if exact05 > exact10:
                raise RuntimeError(
                    f"{locus_id}: IDR05 exact overlap "
                    "exceeds IDR10"
                )

            if n_exact05 > n_exact10:
                raise RuntimeError(
                    f"{locus_id}: IDR05 exact peak count "
                    "exceeds IDR10"
                )

            if bp_exact05 > bp_exact10:
                raise RuntimeError(
                    f"{locus_id}: IDR05 exact accessible bp "
                    "exceeds IDR10"
                )

            if bp_exact10 > width:
                raise RuntimeError(
                    f"{locus_id}: IDR10 exact accessible bp "
                    "exceeds locus width"
                )

            if bp_exact05 > width:
                raise RuntimeError(
                    f"{locus_id}: IDR05 exact accessible bp "
                    "exceeds locus width"
                )

            dist10 = as_int(
                s,
                "atac_idr10_nearest_peak_distance_bp",
            )

            dist05 = as_int(
                s,
                "atac_idr05_nearest_peak_distance_bp",
            )

            if dist05 < dist10:
                raise RuntimeError(
                    f"{locus_id}: nearest IDR05 peak is "
                    "closer than nearest IDR10 peak"
                )

            for window in WINDOWS:

                n10 = as_int(
                    s,
                    f"atac_idr10_n_unique_peaks_pm{window}bp",
                )

                n05 = as_int(
                    s,
                    f"atac_idr05_n_unique_peaks_pm{window}bp",
                )

                bp10 = as_int(
                    s,
                    f"atac_idr10_accessible_bp_pm{window}bp",
                )

                bp05 = as_int(
                    s,
                    f"atac_idr05_accessible_bp_pm{window}bp",
                )

                if n05 > n10:
                    raise RuntimeError(
                        f"{locus_id}: IDR05 peak count exceeds "
                        f"IDR10 in ±{window} bp"
                    )

                if bp05 > bp10:
                    raise RuntimeError(
                        f"{locus_id}: IDR05 accessible bp "
                        f"exceeds IDR10 in ±{window} bp"
                    )

                signal = as_float(
                    s,
                    f"atac_spmr_consensus_pm{window}bp_mean",
                )

                if signal < 0:
                    raise RuntimeError(
                        f"{locus_id}: negative SPMR signal"
                    )

            # ------------------------------------------------
            # Build final row
            # ------------------------------------------------

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
                    b["target_start"],
                "target_end_1based":
                    b["target_end"],
                "locus_width_bp":
                    str(width),
                "locus_semantic_class":
                    semantic_class,
                "mapping_confidence":
                    b["mapping_confidence"],
                "supporting_studies":
                    b["supporting_studies"],

                "atac_anchor_1based":
                    s["anchor_1based"],

                "atac_idr10_exact_overlap":
                    s["atac_idr10_exact_overlap"],
                "atac_idr10_exact_n_unique_peaks":
                    s[
                        "atac_idr10_exact_n_unique_peaks"
                    ],
                "atac_idr10_exact_overlap_bp":
                    s["atac_idr10_exact_overlap_bp"],
                "atac_idr10_exact_overlap_fraction_of_locus":
                    format(
                        bp_exact10 / width,
                        ".17g",
                    ),
                "atac_idr10_anchor_in_peak":
                    s["atac_idr10_anchor_in_peak"],
                "atac_idr10_nearest_peak_distance_bp":
                    s[
                        "atac_idr10_nearest_peak_distance_bp"
                    ],

                "atac_idr05_exact_overlap":
                    s["atac_idr05_exact_overlap"],
                "atac_idr05_exact_n_unique_peaks":
                    s[
                        "atac_idr05_exact_n_unique_peaks"
                    ],
                "atac_idr05_exact_overlap_bp":
                    s["atac_idr05_exact_overlap_bp"],
                "atac_idr05_exact_overlap_fraction_of_locus":
                    format(
                        bp_exact05 / width,
                        ".17g",
                    ),
                "atac_idr05_anchor_in_peak":
                    s["atac_idr05_anchor_in_peak"],
                "atac_idr05_nearest_peak_distance_bp":
                    s[
                        "atac_idr05_nearest_peak_distance_bp"
                    ],
            }

            for window in WINDOWS:

                n10 = as_int(
                    s,
                    f"atac_idr10_n_unique_peaks_pm{window}bp",
                )

                n05 = as_int(
                    s,
                    f"atac_idr05_n_unique_peaks_pm{window}bp",
                )

                row[
                    f"atac_idr10_n_unique_peaks_pm{window}bp"
                ] = str(n10)

                row[
                    f"atac_idr10_any_peak_pm{window}bp"
                ] = "1" if n10 > 0 else "0"

                row[
                    f"atac_idr10_accessible_bp_pm{window}bp"
                ] = s[
                    f"atac_idr10_accessible_bp_pm{window}bp"
                ]

                row[
                    f"atac_spmr_consensus_pm{window}bp_mean"
                ] = s[
                    f"atac_spmr_consensus_pm{window}bp_mean"
                ]

                row[
                    f"atac_idr05_n_unique_peaks_pm{window}bp"
                ] = str(n05)

                row[
                    f"atac_idr05_any_peak_pm{window}bp"
                ] = "1" if n05 > 0 else "0"

                row[
                    f"atac_idr05_accessible_bp_pm{window}bp"
                ] = s[
                    f"atac_idr05_accessible_bp_pm{window}bp"
                ]

            final_rows.append(row)

    # --------------------------------------------------------
    # 4. Final validation
    # --------------------------------------------------------

    if len(final_rows) != 37:
        raise RuntimeError(
            f"Expected 37 final rows; "
            f"found {len(final_rows)}"
        )

    if len({
        row["canonical_locus_id"]
        for row in final_rows
    }) != 37:
        raise RuntimeError(
            "Duplicate canonical locus IDs"
        )

    # --------------------------------------------------------
    # 5. Write output
    # --------------------------------------------------------

    with open(
        args.output,
        "w",
        encoding="utf-8",
        newline="",
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(final_rows)

    # --------------------------------------------------------
    # 6. QC
    # --------------------------------------------------------

    exact10_n = sum(
        int(row["atac_idr10_exact_overlap"])
        for row in final_rows
    )

    exact05_n = sum(
        int(row["atac_idr05_exact_overlap"])
        for row in final_rows
    )

    anchor10_n = sum(
        int(row["atac_idr10_anchor_in_peak"])
        for row in final_rows
    )

    anchor05_n = sum(
        int(row["atac_idr05_anchor_in_peak"])
        for row in final_rows
    )

    with open(
        args.qc,
        "w",
        encoding="utf-8",
    ) as fh:

        fh.write("metric\tvalue\n")
        fh.write("benchmark_loci\t37\n")

        for role in sorted(roles):
            fh.write(
                f"benchmark_role_{role}\t"
                f"{roles[role]}\n"
            )

        semantic_counts = Counter(
            row["locus_semantic_class"]
            for row in final_rows
        )

        for semantic_class in sorted(semantic_counts):
            fh.write(
                f"semantic_{semantic_class}\t"
                f"{semantic_counts[semantic_class]}\n"
            )

        fh.write(
            f"idr10_exact_overlap_loci\t{exact10_n}\n"
        )

        fh.write(
            f"idr05_exact_overlap_loci\t{exact05_n}\n"
        )

        fh.write(
            f"idr10_anchor_in_peak_loci\t{anchor10_n}\n"
        )

        fh.write(
            f"idr05_anchor_in_peak_loci\t{anchor05_n}\n"
        )

        for window in WINDOWS:

            idr10_any = sum(
                int(
                    row[
                        f"atac_idr10_any_peak_pm{window}bp"
                    ]
                )
                for row in final_rows
            )

            idr05_any = sum(
                int(
                    row[
                        f"atac_idr05_any_peak_pm{window}bp"
                    ]
                )
                for row in final_rows
            )

            fh.write(
                f"idr10_any_peak_pm{window}bp_loci\t"
                f"{idr10_any}\n"
            )

            fh.write(
                f"idr05_any_peak_pm{window}bp_loci\t"
                f"{idr05_any}\n"
            )

    # --------------------------------------------------------
    # 7. Provenance
    # --------------------------------------------------------

    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
        ).strip()
    except Exception:
        commit = "unavailable"

    provenance = {
        "product":
            "benchmark_37_atac_context_v1",

        "description":
            (
                "Reproducible ATAC-seq regulatory context "
                "for all 37 canonical CHO benchmark loci."
            ),

        "pipeline_commit":
            commit,

        "method": {
            "biological_dataset":
                "PRJNA667472 CHO-K1 bulk ATAC-seq",

            "replicates": [
                "SRR12774931",
                "SRR12774932",
            ],

            "primary_reproducible_accessibility":
                "IDR <= 0.10",

            "sensitivity_reproducible_accessibility":
                "IDR <= 0.05",

            "complete_locus_features":
                (
                    "Peak overlap evaluated against the full "
                    "canonical 1-based closed locus interval."
                ),

            "anchor_definition":
                (
                    "Integer centre of the canonical 1-based "
                    "closed locus interval; representative "
                    "context point only, not assumed to be the "
                    "exact integration coordinate for mapped "
                    "interval loci."
                ),

            "anchor_half_windows_bp":
                list(WINDOWS),

            "continuous_signal":
                (
                    "Arithmetic mean of replicate MACS3 "
                    "treatment-pileup SPMR window means."
                ),

            "interpretation":
                (
                    "Exploratory regulatory-context "
                    "characterization; no predictive or "
                    "causal claim."
                ),
        },

        "sources": {
            "benchmark": {
                "path":
                    args.benchmark,
                "sha256":
                    sha256_file(args.benchmark),
            },

            "locus_audit": {
                "path":
                    args.locus_audit,
                "sha256":
                    sha256_file(args.locus_audit),
            },

            "all37_atac_features_with_signal": {
                "path":
                    args.source,
                "sha256":
                    sha256_file(args.source),
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

    # --------------------------------------------------------
    # 8. Console
    # --------------------------------------------------------

    print("BENCHMARK_LOCI=37")
    print(
        f"OUTPUT_COLUMNS={len(columns)}"
    )

    print(
        f"IDR10_EXACT_OVERLAP_LOCI={exact10_n}"
    )

    print(
        f"IDR05_EXACT_OVERLAP_LOCI={exact05_n}"
    )

    print(
        f"IDR10_ANCHOR_IN_PEAK_LOCI={anchor10_n}"
    )

    print(
        f"IDR05_ANCHOR_IN_PEAK_LOCI={anchor05_n}"
    )

    print()
    print(
        "PASS: final ATAC context table built "
        "for all 37 loci"
    )


if __name__ == "__main__":
    main()
