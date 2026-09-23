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


EXPECTED_INPUT_ROWS = 333
EXPECTED_INPUT_COLUMNS = 75
EXPECTED_OUTPUT_ROWS = 222
EXPECTED_OUTPUT_COLUMNS = 100
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


PEAK_CONTINUOUS_FEATURES = [
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


TWO_STUDY_MARKS = {
    "H3K27ac",
    "H3K4me3",
    "H3K9me3",
}

ONE_STUDY_MARKS = {
    "H3K27me3",
    "H3K36me3",
    "H3K4me1",
}


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
        fail(f"Non-finite value: {value}")

    return format(value, ".12g")


def read_tsv(path):
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

    return fields, rows


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


def mean(values):
    return sum(values) / len(values)


def parse_float(value, context):
    try:
        x = float(value)
    except ValueError:
        fail(
            f"{context}: non-numeric value {value!r}"
        )

    if not math.isfinite(x):
        fail(
            f"{context}: non-finite value {value!r}"
        )

    return x


def rank_percentiles(values):
    """
    Average-rank percentile in [0, 1].

    Ranks are 1-based internally.
    Tied values receive their average rank.

    percentile = (average_rank - 1) / (n - 1)
    """

    n = len(values)

    if n != EXPECTED_LOCI:
        fail(
            f"Expected {EXPECTED_LOCI} values for "
            f"rank harmonization; found {n}"
        )

    order = sorted(
        range(n),
        key=lambda i: values[i],
    )

    ranks = [None] * n

    pos = 0

    while pos < n:

        end = pos + 1

        while (
            end < n
            and values[order[end]]
            == values[order[pos]]
        ):
            end += 1

        average_rank = (
            (pos + 1) + end
        ) / 2.0

        for j in range(pos, end):
            ranks[order[j]] = average_rank

        pos = end

    return [
        (rank - 1.0) / (n - 1.0)
        for rank in ranks
    ]


def required_fields():
    required = set(
        IDENTITY_FIELDS
        + [
            "target",
            "study",
            "peak_mode",
            "n_conditions",
            "condition_keys",
            "analysis_ids",
        ]
    )

    for feature in BINARY_FEATURES:
        required.update([
            f"{feature}_n_conditions_positive",
            f"{feature}_fraction_conditions_positive",
            f"{feature}_any_condition",
            f"{feature}_all_conditions",
        ])

    for feature in PEAK_CONTINUOUS_FEATURES:
        required.add(
            f"{feature}_study_median"
        )

    for feature in SPMR_FEATURES:
        required.add(
            f"{feature}_study_median"
        )

    return required


def build_output_fields():
    fields = list(IDENTITY_FIELDS)

    fields += [
        "target",
        "peak_mode",
        "n_studies",
        "studies",
        "n_conditions_total",
        "source_condition_keys",
        "source_analysis_ids",
    ]

    for feature in BINARY_FEATURES:
        fields += [
            f"{feature}_n_positive_conditions_total",
            f"{feature}_fraction_positive_conditions_raw",
            f"{feature}_fraction_positive_conditions_equal_study_weight",
            f"{feature}_fraction_studies_with_any_condition",
            f"{feature}_fraction_studies_with_all_conditions",
            f"{feature}_any_study",
            f"{feature}_all_studies_have_any_condition",
        ]

    for feature in PEAK_CONTINUOUS_FEATURES:
        fields += [
            f"{feature}_consensus_median_across_studies",
            f"{feature}_cross_study_range",
        ]

    for feature in SPMR_FEATURES:
        fields += [
            f"{feature}_relative_percentile_equal_study_consensus",
            f"{feature}_relative_percentile_min",
            f"{feature}_relative_percentile_max",
            f"{feature}_relative_percentile_cross_study_abs_diff",
        ]

    return fields


def build_percentile_lookup(rows):
    """
    Rank harmonization is performed separately for each
    histone-mark x study x SPMR feature over the 37 loci.
    """

    mark_study = defaultdict(list)

    for row in rows:
        mark_study[
            (
                row["target"],
                row["study"],
            )
        ].append(row)

    lookup = {}

    for (mark, study), rr in sorted(mark_study.items()):

        if len(rr) != EXPECTED_LOCI:
            fail(
                f"{mark}/{study}: expected "
                f"{EXPECTED_LOCI} loci; found {len(rr)}"
            )

        rr = sorted(
            rr,
            key=lambda r: r["canonical_locus_id"],
        )

        locus_ids = [
            r["canonical_locus_id"]
            for r in rr
        ]

        if len(set(locus_ids)) != EXPECTED_LOCI:
            fail(
                f"{mark}/{study}: duplicate locus"
            )

        for feature in SPMR_FEATURES:

            source_field = (
                f"{feature}_study_median"
            )

            values = [
                parse_float(
                    r[source_field],
                    f"{mark}/{study}/{feature}",
                )
                for r in rr
            ]

            percentiles = rank_percentiles(
                values
            )

            for row, percentile in zip(
                rr,
                percentiles,
            ):

                lookup[
                    (
                        row["canonical_locus_id"],
                        mark,
                        study,
                        feature,
                    )
                ] = percentile

    return lookup


def summarize_group(
    locus,
    mark,
    rows,
    percentile_lookup,
):
    context = f"{locus}/{mark}"

    expected_studies = (
        2
        if mark in TWO_STUDY_MARKS
        else 1
    )

    if len(rows) != expected_studies:
        fail(
            f"{context}: expected "
            f"{expected_studies} studies; "
            f"found {len(rows)}"
        )

    rows = sorted(
        rows,
        key=lambda r: r["study"],
    )

    studies = [
        r["study"]
        for r in rows
    ]

    if len(set(studies)) != len(studies):
        fail(
            f"{context}: duplicate study"
        )

    out = {}

    for field in IDENTITY_FIELDS:

        values = {
            r[field]
            for r in rows
        }

        if len(values) != 1:
            fail(
                f"{context}: identity field "
                f"{field} varies across studies"
            )

        out[field] = next(iter(values))

    peak_modes = {
        r["peak_mode"]
        for r in rows
    }

    if len(peak_modes) != 1:
        fail(
            f"{context}: peak_mode differs "
            "between studies"
        )

    out["target"] = mark
    out["peak_mode"] = next(iter(peak_modes))
    out["n_studies"] = str(len(rows))
    out["studies"] = ";".join(studies)

    n_conditions_total = sum(
        int(r["n_conditions"])
        for r in rows
    )

    out["n_conditions_total"] = str(
        n_conditions_total
    )

    condition_keys = sorted({
        condition
        for r in rows
        for condition in r["condition_keys"].split(";")
        if condition
    })

    analysis_ids = sorted({
        analysis
        for r in rows
        for analysis in r["analysis_ids"].split(";")
        if analysis
    })

    out["source_condition_keys"] = ";".join(
        condition_keys
    )

    out["source_analysis_ids"] = ";".join(
        analysis_ids
    )

    # --------------------------------------------------------
    # Peak presence / binary evidence.
    # --------------------------------------------------------

    for feature in BINARY_FEATURES:

        n_positive_by_study = [
            int(
                r[
                    f"{feature}_n_conditions_positive"
                ]
            )
            for r in rows
        ]

        condition_fraction_by_study = [
            parse_float(
                r[
                    f"{feature}_fraction_conditions_positive"
                ],
                context,
            )
            for r in rows
        ]

        any_by_study = [
            int(
                r[
                    f"{feature}_any_condition"
                ]
            )
            for r in rows
        ]

        all_by_study = [
            int(
                r[
                    f"{feature}_all_conditions"
                ]
            )
            for r in rows
        ]

        if any(
            x not in (0, 1)
            for x in any_by_study
        ):
            fail(
                f"{context}: invalid any_condition "
                f"for {feature}"
            )

        if any(
            x not in (0, 1)
            for x in all_by_study
        ):
            fail(
                f"{context}: invalid all_conditions "
                f"for {feature}"
            )

        total_positive = sum(
            n_positive_by_study
        )

        raw_fraction = (
            total_positive
            / n_conditions_total
        )

        equal_study_fraction = mean(
            condition_fraction_by_study
        )

        fraction_studies_any = mean(
            any_by_study
        )

        fraction_studies_all = mean(
            all_by_study
        )

        out[
            f"{feature}_n_positive_conditions_total"
        ] = str(total_positive)

        out[
            f"{feature}_fraction_positive_conditions_raw"
        ] = fmt(raw_fraction)

        out[
            f"{feature}_fraction_positive_conditions_equal_study_weight"
        ] = fmt(equal_study_fraction)

        out[
            f"{feature}_fraction_studies_with_any_condition"
        ] = fmt(fraction_studies_any)

        out[
            f"{feature}_fraction_studies_with_all_conditions"
        ] = fmt(fraction_studies_all)

        out[
            f"{feature}_any_study"
        ] = str(
            int(any(any_by_study))
        )

        out[
            f"{feature}_all_studies_have_any_condition"
        ] = str(
            int(all(any_by_study))
        )

    # --------------------------------------------------------
    # Peak geometry / distance / enriched fraction.
    #
    # These values are in common genomic units and derive
    # from the same peak-calling workflow. We therefore use
    # the median of the study-level medians as consensus.
    # --------------------------------------------------------

    for feature in PEAK_CONTINUOUS_FEATURES:

        field = (
            f"{feature}_study_median"
        )

        values = [
            parse_float(
                r[field],
                f"{context}/{field}",
            )
            for r in rows
        ]

        out[
            f"{feature}_consensus_median_across_studies"
        ] = fmt(
            statistics.median(values)
        )

        if len(values) == 2:
            out[
                f"{feature}_cross_study_range"
            ] = fmt(
                max(values) - min(values)
            )
        else:
            out[
                f"{feature}_cross_study_range"
            ] = "NA"

    # --------------------------------------------------------
    # SPMR.
    #
    # Raw SPMR values are deliberately NOT averaged across
    # studies. Study-specific SPMR remains preserved in the
    # 333-row source table.
    #
    # Here we combine only within-study relative percentiles.
    # --------------------------------------------------------

    for feature in SPMR_FEATURES:

        percentiles = [
            percentile_lookup[
                (
                    locus,
                    mark,
                    r["study"],
                    feature,
                )
            ]
            for r in rows
        ]

        consensus = mean(percentiles)

        out[
            f"{feature}_relative_percentile_equal_study_consensus"
        ] = fmt(consensus)

        out[
            f"{feature}_relative_percentile_min"
        ] = fmt(
            min(percentiles)
        )

        out[
            f"{feature}_relative_percentile_max"
        ] = fmt(
            max(percentiles)
        )

        if len(percentiles) == 2:
            out[
                f"{feature}_relative_percentile_cross_study_abs_diff"
            ] = fmt(
                abs(
                    percentiles[0]
                    - percentiles[1]
                )
            )
        else:
            out[
                f"{feature}_relative_percentile_cross_study_abs_diff"
            ] = "NA"

    return out


def validate_output(rows, fields):
    errors = []

    if len(rows) != EXPECTED_OUTPUT_ROWS:
        errors.append(
            f"rows={len(rows)}"
        )

    if len(fields) != EXPECTED_OUTPUT_COLUMNS:
        errors.append(
            f"columns={len(fields)}"
        )

    loci = {
        r["canonical_locus_id"]
        for r in rows
    }

    marks = {
        r["target"]
        for r in rows
    }

    pairs = {
        (
            r["canonical_locus_id"],
            r["target"],
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

    if len(pairs) != EXPECTED_OUTPUT_ROWS:
        errors.append(
            f"unique_pairs={len(pairs)}"
        )

    mark_counts = Counter(
        r["target"]
        for r in rows
    )

    if set(mark_counts.values()) != {37}:
        errors.append(
            "mark_count_not_37"
        )

    n_study_counts = Counter(
        int(r["n_studies"])
        for r in rows
    )

    if n_study_counts != Counter({
        1: 111,
        2: 111,
    }):
        errors.append(
            f"unexpected_n_study_distribution="
            f"{dict(n_study_counts)}"
        )

    for row in rows:

        mark = row["target"]
        n_studies = int(
            row["n_studies"]
        )

        expected = (
            2
            if mark in TWO_STUDY_MARKS
            else 1
        )

        if n_studies != expected:
            errors.append(
                f"{row['canonical_locus_id']}/"
                f"{mark}:wrong_n_studies"
            )

        expected_conditions = (
            4
            if expected == 2
            else 2
        )

        if (
            int(row["n_conditions_total"])
            != expected_conditions
        ):
            errors.append(
                f"{row['canonical_locus_id']}/"
                f"{mark}:wrong_n_conditions"
            )

        for feature in BINARY_FEATURES:

            for suffix in (
                "fraction_positive_conditions_raw",
                "fraction_positive_conditions_equal_study_weight",
                "fraction_studies_with_any_condition",
                "fraction_studies_with_all_conditions",
            ):

                value = float(
                    row[
                        f"{feature}_{suffix}"
                    ]
                )

                if not (
                    -1e-12
                    <= value
                    <= 1.0 + 1e-12
                ):
                    errors.append(
                        f"{feature}/{suffix}:"
                        "outside_0_1"
                    )

            for suffix in (
                "any_study",
                "all_studies_have_any_condition",
            ):

                value = int(
                    row[
                        f"{feature}_{suffix}"
                    ]
                )

                if value not in (0, 1):
                    errors.append(
                        f"{feature}/{suffix}:"
                        "not_binary"
                    )

        for feature in SPMR_FEATURES:

            consensus = float(
                row[
                    f"{feature}_relative_percentile_equal_study_consensus"
                ]
            )

            minimum = float(
                row[
                    f"{feature}_relative_percentile_min"
                ]
            )

            maximum = float(
                row[
                    f"{feature}_relative_percentile_max"
                ]
            )

            for value in (
                consensus,
                minimum,
                maximum,
            ):
                if not (
                    -1e-12
                    <= value
                    <= 1.0 + 1e-12
                ):
                    errors.append(
                        f"{feature}:percentile_"
                        "outside_0_1"
                    )

            diff_field = (
                f"{feature}_relative_percentile_"
                "cross_study_abs_diff"
            )

            diff = row[diff_field]

            if n_studies == 2:

                if diff == "":
                    errors.append(
                        f"{feature}:missing_diff"
                    )
                else:
                    observed = float(diff)

                    expected_diff = (
                        maximum - minimum
                    )

                    if not math.isclose(
                        observed,
                        expected_diff,
                        rel_tol=1e-9,
                        abs_tol=1e-10,
                    ):
                        errors.append(
                            f"{feature}:diff_formula"
                        )

            else:
                if diff != "NA":
                    errors.append(
                        f"{feature}:single_study_"
                        "diff_should_be_NA"
                    )

    if errors:
        fail(
            "Output validation failed: "
            + "; ".join(errors[:30])
        )

    return {
        "rows": len(rows),
        "columns": len(fields),
        "unique_loci": len(loci),
        "unique_marks": len(marks),
        "unique_locus_mark_pairs":
            len(pairs),
        "one_study_rows":
            n_study_counts[1],
        "two_study_rows":
            n_study_counts[2],
        "H3K27ac_rows":
            mark_counts["H3K27ac"],
        "H3K27me3_rows":
            mark_counts["H3K27me3"],
        "H3K36me3_rows":
            mark_counts["H3K36me3"],
        "H3K4me1_rows":
            mark_counts["H3K4me1"],
        "H3K4me3_rows":
            mark_counts["H3K4me3"],
        "H3K9me3_rows":
            mark_counts["H3K9me3"],
        "validation_errors": 0,
        "status": "PASS",
    }


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
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

    parser.add_argument(
        "--sha256s",
        required=True,
    )

    return parser.parse_args()


def main():
    args = parse_args()

    output_paths = [
        Path(args.output),
        Path(args.qc),
        Path(args.provenance),
        Path(args.sha256s),
    ]

    existing = [
        str(path)
        for path in output_paths
        if path.exists()
    ]

    if existing:
        fail(
            "Refusing to overwrite outputs: "
            + ", ".join(existing)
        )

    fields, rows = read_tsv(
        args.input
    )

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

    missing = sorted(
        required_fields()
        - set(fields)
    )

    if missing:
        fail(
            "Missing required fields: "
            + ", ".join(missing)
        )

    percentile_lookup = (
        build_percentile_lookup(rows)
    )

    grouped = defaultdict(list)

    for row in rows:

        grouped[
            (
                row["canonical_locus_id"],
                row["target"],
            )
        ].append(row)

    if len(grouped) != EXPECTED_OUTPUT_ROWS:
        fail(
            f"Expected {EXPECTED_OUTPUT_ROWS} "
            f"locus-mark groups; found {len(grouped)}"
        )

    output_rows = []

    for (
        locus,
        mark,
    ) in sorted(grouped):

        output_rows.append(
            summarize_group(
                locus,
                mark,
                grouped[(locus, mark)],
                percentile_lookup,
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
            "ChIP-seq locus x histone-mark "
            "cross-study consensus"
        ),
        "source_table":
            str(Path(args.input)),
        "source_table_sha256":
            sha256_file(args.input),
        "source_rows":
            len(rows),
        "source_columns":
            len(fields),
        "grouping": [
            "canonical_locus_id",
            "target",
        ],
        "design": {
            "two_study_marks": sorted(
                TWO_STUDY_MARKS
            ),
            "one_study_marks": sorted(
                ONE_STUDY_MARKS
            ),
            "study_weighting":
                "equal weight per study",
        },
        "peak_feature_policy": {
            "binary_features": (
                "Preserve condition-level prevalence "
                "and study-level support. Equal-study "
                "weighted condition fractions are "
                "reported explicitly."
            ),
            "continuous_peak_features": (
                "Median of study-level medians. "
                "Cross-study range is reported only "
                "when two studies are available."
            ),
        },
        "spmr_policy": {
            "raw_cross_study_average":
                "NOT performed",
            "raw_values_preserved_in":
                str(Path(args.input)),
            "harmonization": (
                "Within each histone-mark x study, "
                "the 37 loci are converted to "
                "average-rank percentiles using "
                "(rank-1)/(n-1), with mid-ranks "
                "for ties."
            ),
            "cross_study_consensus": (
                "Arithmetic mean of within-study "
                "percentiles with equal study weight."
            ),
            "cross_study_disagreement": (
                "Absolute difference between the two "
                "study percentiles when two studies "
                "are available; NA for single-study "
                "marks because cross-study disagreement "
                "is not applicable."
            ),
            "primary_continuous_features": [
                "chip_spmr_pm1000bp_mean",
                "chip_spmr_pm10000bp_mean",
                "chip_spmr_pm50000bp_mean",
            ],
            "secondary_context_feature": {
                "feature":
                    "chip_spmr_locus_mean",
                "reason": (
                    "Benchmark loci have heterogeneous "
                    "interval widths and exact-locus "
                    "SPMR showed more zero/tied values "
                    "and weaker cross-study agreement."
                ),
            },
        },
        "interpretation": (
            "Percentile consensus is a relative "
            "within-mark regulatory-context measure; "
            "it is not an absolute cross-study "
            "quantification of ChIP enrichment."
        ),
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

        for path in (
            Path(args.output),
            Path(args.qc),
            Path(args.provenance),
        ):

            fh.write(
                f"{sha256_file(path)}  "
                f"{path.name}\n"
            )

    os.replace(
        sha_tmp,
        args.sha256s,
    )

    print(
        "CHIP_MARK_CONSENSUS=PASS"
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
