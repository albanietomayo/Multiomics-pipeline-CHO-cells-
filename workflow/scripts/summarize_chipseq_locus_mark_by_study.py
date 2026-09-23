#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import math
import os
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED_INPUT_ROWS = 666
EXPECTED_INPUT_COLUMNS = 62
EXPECTED_OUTPUT_ROWS = 333
EXPECTED_LOCI = 37
EXPECTED_MARKS = 6


IDENTITY_FIELDS = [
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
    "chip_anchor_1based",
]


BINARY_FEATURES = [
    "chip_exact_overlap",
    "chip_anchor_in_peak",
    "chip_any_peak_pm1000bp",
    "chip_any_peak_pm10000bp",
    "chip_any_peak_pm50000bp",
]


MEDIAN_FEATURES = [
    "chip_exact_n_peaks",
    "chip_exact_overlap_bp",
    "chip_exact_overlap_fraction_of_locus",
    "chip_nearest_peak_distance_to_anchor_bp",
    "chip_nearest_peak_distance_to_locus_bp",

    "chip_n_peaks_pm1000bp",
    "chip_enriched_bp_pm1000bp",
    "chip_enriched_fraction_pm1000bp",

    "chip_n_peaks_pm10000bp",
    "chip_enriched_bp_pm10000bp",
    "chip_enriched_fraction_pm10000bp",

    "chip_n_peaks_pm50000bp",
    "chip_enriched_bp_pm50000bp",
    "chip_enriched_fraction_pm50000bp",
]


SPMR_FEATURES = [
    "chip_spmr_locus_mean",
    "chip_spmr_pm1000bp_mean",
    "chip_spmr_pm10000bp_mean",
    "chip_spmr_pm50000bp_mean",
]


ANALYSIS_QC_FEATURES = [
    "analysis_peak_count",
    "analysis_frip",
    "analysis_input_overlap_fraction",
    "analysis_ip_input_overlap_ratio",
    "analysis_median_fold_enrichment",
]


def fail(message):
    raise RuntimeError(message)


def sha256_file(path):
    h = hashlib.sha256()

    with Path(path).open("rb") as fh:
        for chunk in iter(
            lambda: fh.read(8 * 1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


def fmt(value):
    if not math.isfinite(value):
        fail(f"Non-finite output value: {value}")

    return format(value, ".12g")


def parse_numeric(value, field, context):
    try:
        x = float(value)
    except ValueError:
        fail(
            f"{context}: non-numeric {field}={value!r}"
        )

    if not math.isfinite(x):
        fail(
            f"{context}: non-finite {field}={value!r}"
        )

    return x


def median_numeric(rows, field, context):
    values = [
        parse_numeric(
            r[field],
            field,
            context,
        )
        for r in rows
    ]

    return statistics.median(values)


def write_atomic_tsv(path, fields, rows):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")

    with tmp.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)

    os.replace(tmp, path)


def write_atomic_json(path, obj):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")

    with tmp.open(
        "w",
        encoding="utf-8",
    ) as fh:
        json.dump(
            obj,
            fh,
            indent=2,
            sort_keys=True,
        )
        fh.write("\n")

    os.replace(tmp, path)


def load_input(path):
    with Path(path).open(
        "r",
        encoding="utf-8",
        newline="",
    ) as fh:
        reader = csv.DictReader(
            fh,
            delimiter="\t",
        )

        fields = reader.fieldnames or []
        rows = list(reader)

    if len(rows) != EXPECTED_INPUT_ROWS:
        fail(
            f"Expected {EXPECTED_INPUT_ROWS} rows; "
            f"found {len(rows)}"
        )

    if len(fields) != EXPECTED_INPUT_COLUMNS:
        fail(
            f"Expected {EXPECTED_INPUT_COLUMNS} columns; "
            f"found {len(fields)}"
        )

    required = set(
        IDENTITY_FIELDS
        + BINARY_FEATURES
        + MEDIAN_FEATURES
        + SPMR_FEATURES
        + ANALYSIS_QC_FEATURES
        + [
            "analysis_id",
            "study",
            "condition_key",
            "target",
            "peak_mode",
            "chip_spmr_locus_covered_bp",
            "chip_spmr_locus_covered_fraction",
            "chip_spmr_locus_actual_bp",
            "chip_spmr_pm1000bp_covered_bp",
            "chip_spmr_pm1000bp_covered_fraction",
            "chip_spmr_pm1000bp_actual_window_bp",
            "chip_spmr_pm10000bp_covered_bp",
            "chip_spmr_pm10000bp_covered_fraction",
            "chip_spmr_pm10000bp_actual_window_bp",
            "chip_spmr_pm50000bp_covered_bp",
            "chip_spmr_pm50000bp_covered_fraction",
            "chip_spmr_pm50000bp_actual_window_bp",
        ]
    )

    missing = sorted(required - set(fields))

    if missing:
        fail(
            "Missing required fields: "
            + ", ".join(missing)
        )

    return fields, rows


def validate_spmr_technical_invariants(rows):
    errors = 0

    for row in rows:

        locus_bp = int(
            row["chip_spmr_locus_actual_bp"]
        )

        locus_covered = int(
            row["chip_spmr_locus_covered_bp"]
        )

        locus_fraction = float(
            row[
                "chip_spmr_locus_covered_fraction"
            ]
        )

        expected_locus_bp = (
            int(row["target_end_1based"])
            - int(row["target_start_1based"])
            + 1
        )

        if locus_bp != expected_locus_bp:
            errors += 1

        if locus_covered != locus_bp:
            errors += 1

        if not math.isclose(
            locus_fraction,
            1.0,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            errors += 1

        for half_window, expected_bp in (
            (1000, 2001),
            (10000, 20001),
            (50000, 100001),
        ):
            actual = int(
                row[
                    f"chip_spmr_pm{half_window}bp_"
                    "actual_window_bp"
                ]
            )

            covered = int(
                row[
                    f"chip_spmr_pm{half_window}bp_"
                    "covered_bp"
                ]
            )

            fraction = float(
                row[
                    f"chip_spmr_pm{half_window}bp_"
                    "covered_fraction"
                ]
            )

            if actual != expected_bp:
                errors += 1

            if covered != actual:
                errors += 1

            if not math.isclose(
                fraction,
                1.0,
                rel_tol=0.0,
                abs_tol=1e-12,
            ):
                errors += 1

    if errors:
        fail(
            "SPMR technical invariant errors="
            f"{errors}"
        )

    return errors


def build_output_fields():
    fields = list(IDENTITY_FIELDS)

    fields += [
        "target",
        "study",
        "peak_mode",
        "n_conditions",
        "condition_keys",
        "analysis_ids",
    ]

    for feature in ANALYSIS_QC_FEATURES:
        fields.append(
            f"{feature}_study_median"
        )

    for feature in BINARY_FEATURES:
        fields.extend([
            f"{feature}_n_conditions_positive",
            f"{feature}_fraction_conditions_positive",
            f"{feature}_any_condition",
            f"{feature}_all_conditions",
        ])

    for feature in MEDIAN_FEATURES:
        fields.append(
            f"{feature}_study_median"
        )

    for feature in SPMR_FEATURES:
        fields.extend([
            f"{feature}_study_median",
            f"{feature}_study_min",
            f"{feature}_study_max",
            f"{feature}_study_abs_delta",
        ])

    return fields


def summarize_group(
    locus_id,
    mark,
    study,
    rows,
):
    context = (
        f"{locus_id}/{mark}/{study}"
    )

    if len(rows) != 2:
        fail(
            f"{context}: expected exactly "
            f"2 condition rows; found {len(rows)}"
        )

    condition_keys = sorted({
        r["condition_key"]
        for r in rows
    })

    analysis_ids = sorted({
        r["analysis_id"]
        for r in rows
    })

    if len(condition_keys) != 2:
        fail(
            f"{context}: expected 2 conditions; "
            f"found {len(condition_keys)}"
        )

    if len(analysis_ids) != 2:
        fail(
            f"{context}: expected 2 analyses; "
            f"found {len(analysis_ids)}"
        )

    out = {}

    for field in IDENTITY_FIELDS:
        unique = {
            r[field]
            for r in rows
        }

        if len(unique) != 1:
            fail(
                f"{context}: identity field "
                f"{field} varies within group"
            )

        out[field] = next(iter(unique))

    peak_modes = {
        r["peak_mode"]
        for r in rows
    }

    if len(peak_modes) != 1:
        fail(
            f"{context}: peak_mode varies "
            "within study"
        )

    out["target"] = mark
    out["study"] = study
    out["peak_mode"] = next(iter(peak_modes))
    out["n_conditions"] = "2"

    out["condition_keys"] = ";".join(
        condition_keys
    )

    out["analysis_ids"] = ";".join(
        analysis_ids
    )

    for feature in ANALYSIS_QC_FEATURES:
        out[
            f"{feature}_study_median"
        ] = fmt(
            median_numeric(
                rows,
                feature,
                context,
            )
        )

    for feature in BINARY_FEATURES:

        values = [
            int(
                parse_numeric(
                    r[feature],
                    feature,
                    context,
                )
            )
            for r in rows
        ]

        if any(
            x not in (0, 1)
            for x in values
        ):
            fail(
                f"{context}: {feature} "
                "is not binary"
            )

        n_positive = sum(values)

        out[
            f"{feature}_n_conditions_positive"
        ] = str(n_positive)

        out[
            f"{feature}_fraction_conditions_positive"
        ] = fmt(
            n_positive / len(values)
        )

        out[
            f"{feature}_any_condition"
        ] = str(
            int(n_positive > 0)
        )

        out[
            f"{feature}_all_conditions"
        ] = str(
            int(n_positive == len(values))
        )

    for feature in MEDIAN_FEATURES:

        out[
            f"{feature}_study_median"
        ] = fmt(
            median_numeric(
                rows,
                feature,
                context,
            )
        )

    for feature in SPMR_FEATURES:

        values = sorted([
            parse_numeric(
                r[feature],
                feature,
                context,
            )
            for r in rows
        ])

        out[
            f"{feature}_study_median"
        ] = fmt(
            statistics.median(values)
        )

        out[
            f"{feature}_study_min"
        ] = fmt(
            min(values)
        )

        out[
            f"{feature}_study_max"
        ] = fmt(
            max(values)
        )

        out[
            f"{feature}_study_abs_delta"
        ] = fmt(
            abs(values[1] - values[0])
        )

    return out


def validate_output(rows, fields):
    errors = []

    if len(rows) != EXPECTED_OUTPUT_ROWS:
        errors.append(
            f"rows={len(rows)}"
        )

    loci = {
        r["canonical_locus_id"]
        for r in rows
    }

    marks = {
        r["target"]
        for r in rows
    }

    groups = {
        (
            r["canonical_locus_id"],
            r["target"],
            r["study"],
        )
        for r in rows
    }

    if len(loci) != EXPECTED_LOCI:
        errors.append(
            f"loci={len(loci)}"
        )

    if len(marks) != EXPECTED_MARKS:
        errors.append(
            f"marks={len(marks)}"
        )

    if len(groups) != EXPECTED_OUTPUT_ROWS:
        errors.append(
            f"unique_groups={len(groups)}"
        )

    if any(
        r["n_conditions"] != "2"
        for r in rows
    ):
        errors.append(
            "n_conditions_not_2"
        )

    mark_counts = Counter(
        r["target"]
        for r in rows
    )

    expected_mark_counts = {
        "H3K27ac": 74,
        "H3K4me3": 74,
        "H3K9me3": 74,
        "H3K27me3": 37,
        "H3K36me3": 37,
        "H3K4me1": 37,
    }

    if dict(mark_counts) != expected_mark_counts:
        errors.append(
            "mark_row_distribution_mismatch"
        )

    for row in rows:

        for feature in BINARY_FEATURES:

            frac = float(
                row[
                    f"{feature}_fraction_conditions_positive"
                ]
            )

            n = int(
                row[
                    f"{feature}_n_conditions_positive"
                ]
            )

            any_value = int(
                row[
                    f"{feature}_any_condition"
                ]
            )

            all_value = int(
                row[
                    f"{feature}_all_conditions"
                ]
            )

            if n not in (0, 1, 2):
                errors.append(
                    f"{feature}:invalid_n"
                )

            if not math.isclose(
                frac,
                n / 2.0,
                rel_tol=0.0,
                abs_tol=1e-12,
            ):
                errors.append(
                    f"{feature}:fraction_formula"
                )

            if any_value != int(n > 0):
                errors.append(
                    f"{feature}:any_formula"
                )

            if all_value != int(n == 2):
                errors.append(
                    f"{feature}:all_formula"
                )

    if errors:
        fail(
            "Output validation failed: "
            + "; ".join(errors[:20])
        )

    return {
        "rows": len(rows),
        "columns": len(fields),
        "unique_loci": len(loci),
        "unique_marks": len(marks),
        "unique_locus_mark_study_groups":
            len(groups),
        "groups_with_two_conditions":
            sum(
                r["n_conditions"] == "2"
                for r in rows
            ),
        "H3K27ac_rows":
            mark_counts["H3K27ac"],
        "H3K4me3_rows":
            mark_counts["H3K4me3"],
        "H3K9me3_rows":
            mark_counts["H3K9me3"],
        "H3K27me3_rows":
            mark_counts["H3K27me3"],
        "H3K36me3_rows":
            mark_counts["H3K36me3"],
        "H3K4me1_rows":
            mark_counts["H3K4me1"],
        "validation_errors": 0,
        "status": "PASS",
    }


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument(
        "--input",
        required=True,
    )

    p.add_argument(
        "--output",
        required=True,
    )

    p.add_argument(
        "--qc",
        required=True,
    )

    p.add_argument(
        "--provenance",
        required=True,
    )

    p.add_argument(
        "--sha256s",
        required=True,
    )

    return p.parse_args()


def main():
    args = parse_args()

    output_paths = [
        Path(args.output),
        Path(args.qc),
        Path(args.provenance),
        Path(args.sha256s),
    ]

    existing = [
        str(p)
        for p in output_paths
        if p.exists()
    ]

    if existing:
        fail(
            "Refusing to overwrite outputs: "
            + ", ".join(existing)
        )

    input_fields, rows = load_input(
        args.input
    )

    technical_errors = (
        validate_spmr_technical_invariants(
            rows
        )
    )

    grouped = defaultdict(list)

    for row in rows:
        key = (
            row["canonical_locus_id"],
            row["target"],
            row["study"],
        )

        grouped[key].append(row)

    output_rows = []

    for (
        locus_id,
        mark,
        study,
    ) in sorted(grouped):

        output_rows.append(
            summarize_group(
                locus_id,
                mark,
                study,
                grouped[
                    (
                        locus_id,
                        mark,
                        study,
                    )
                ],
            )
        )

    output_fields = build_output_fields()

    metrics = validate_output(
        output_rows,
        output_fields,
    )

    write_atomic_tsv(
        args.output,
        output_fields,
        output_rows,
    )

    qc_rows = [
        {
            "metric": key,
            "value": str(value),
        }
        for key, value in metrics.items()
    ]

    write_atomic_tsv(
        args.qc,
        ["metric", "value"],
        qc_rows,
    )

    provenance = {
        "workflow": (
            "ChIP-seq locus x histone mark x study "
            "aggregation"
        ),
        "source_table": str(
            Path(args.input)
        ),
        "source_table_sha256":
            sha256_file(args.input),
        "source_rows":
            len(rows),
        "source_columns":
            len(input_fields),
        "grouping": [
            "canonical_locus_id",
            "target",
            "study",
        ],
        "expected_conditions_per_group":
            2,
        "aggregation_policy": {
            "binary_peak_features": (
                "n positive conditions, fraction "
                "positive, any condition, all conditions"
            ),
            "other_peak_features": (
                "median across the two conditions "
                "within each study"
            ),
            "spmr_features": (
                "median, min, max, and absolute "
                "condition delta within study"
            ),
            "analysis_qc_features": (
                "median across the two "
                "condition-specific analyses"
            ),
            "cross_study_aggregation": (
                "NOT performed at this stage"
            ),
        },
        "technical_spmr_fields": {
            "policy": (
                "validated as invariants and excluded "
                "from biological study-summary features"
            ),
            "full_coverage_errors":
                technical_errors,
            "nominal_windows_bp": {
                "pm1000": 2001,
                "pm10000": 20001,
                "pm50000": 100001,
            },
        },
        "validation": metrics,
    }

    write_atomic_json(
        args.provenance,
        provenance,
    )

    sha_tmp = Path(
        args.sha256s
    ).with_name(
        Path(args.sha256s).name + ".tmp"
    )

    with sha_tmp.open(
        "w",
        encoding="utf-8",
    ) as fh:
        for p in (
            Path(args.output),
            Path(args.qc),
            Path(args.provenance),
        ):
            fh.write(
                f"{sha256_file(p)}  "
                f"{p.name}\n"
            )

    os.replace(
        sha_tmp,
        args.sha256s,
    )

    print(
        "CHIP_STUDY_SUMMARY=PASS"
    )

    for key, value in metrics.items():
        print(
            f"{key}={value}"
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(
            f"ERROR: {exc}",
            file=sys.stderr,
        )
        sys.exit(1)
