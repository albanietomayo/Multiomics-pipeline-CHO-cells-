#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED_ROWS = 37
EXPECTED_COLUMNS = 588

EXPECTED_MODALITY_COUNTS = {
    "CORE": 13,
    "RNA": 30,
    "ATAC": 34,
    "CHIP": 511,
}

EXPECTED_CONSTANT_COLUMNS = 107
EXPECTED_ALL_MISSING_COLUMNS = 54
EXPECTED_COLUMNS_WITH_MISSING = 55
EXPECTED_BENCHMARK_RELATIVE = 96

MARKS = [
    "H3K27ac",
    "H3K27me3",
    "H3K36me3",
    "H3K4me1",
    "H3K4me3",
    "H3K9me3",
]


CORE_FIELDS = [
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
]


RNA_HIGH_CARDINALITY = {
    "rna_nearest_gene_id",
}

RNA_COORDINATE = {
    "rna_nearest_gene_tss_1based",
}

RNA_CATEGORICAL_BIOLOGICAL = {
    "rna_nearest_gene_strand",
    "rna_locus_relative_to_tss",
    "rna_tss_inside_locus",
}

RNA_NUMERIC_BIOLOGICAL = {
    "rna_min_tss_distance_bp",
    "rna_tss_inside_count",
}

RNA_PROVENANCE = {
    "rna_n_studies",
    "rna_n_experiments",
}

RNA_TECHNICAL = {
    "rna_raw_count_total",
    "rna_experiments_raw_count_positive_n",
    "rna_experiments_raw_count_positive_fraction",
    "rna_raw_all_zero",
}

RNA_SENSITIVITY = {
    "rna_tmm_logcpm_median_without_PRJEB35183",
    "rna_tmm_logcpm_delta_without_PRJEB35183",
    "rna_tmm_logcpm_abs_delta_without_PRJEB35183",
    "rna_tpm_median_without_PRJEB35183",
    "rna_tpm_delta_without_PRJEB35183",
    "rna_tpm_abs_delta_without_PRJEB35183",
    "rna_tpm_relative_change_without_PRJEB35183",
}

RNA_EXPRESSION = {
    "rna_tmm_logcpm_median",
    "rna_tmm_logcpm_q1",
    "rna_tmm_logcpm_q3",
    "rna_tmm_logcpm_iqr",
    "rna_tpm_median",
    "rna_tpm_q1",
    "rna_tpm_q3",
    "rna_tpm_iqr",
    "rna_study_median_tpm_positive_n",
    "rna_study_median_tpm_positive_fraction",
}


CHIP_PROVENANCE_SUFFIXES = {
    "peak_mode",
    "n_studies",
    "studies",
    "n_conditions_total",
    "source_condition_keys",
    "source_analysis_ids",
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

        return (
            reader.fieldnames or [],
            list(reader),
        )


def write_atomic_tsv(
    path,
    fields,
    rows,
):

    path = Path(path)

    tmp = path.with_name(
        path.name + ".tmp"
    )

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

    os.replace(
        tmp,
        path,
    )


def write_atomic_json(
    path,
    obj,
):

    path = Path(path)

    tmp = path.with_name(
        path.name + ".tmp"
    )

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

    os.replace(
        tmp,
        path,
    )


def parse_number(value):

    if value in ("", "NA"):
        return None

    try:
        number = float(value)

        if math.isfinite(number):
            return number

        return None

    except ValueError:
        return None


def profile_field(
    field,
    rows,
):

    values = [
        row[field]
        for row in rows
    ]

    nonmissing = [
        value
        for value in values
        if value not in (
            "",
            "NA",
        )
    ]

    parsed = [
        parse_number(value)
        for value in nonmissing
    ]

    if not nonmissing:
        dtype = "all_missing"

    elif all(
        value is not None
        for value in parsed
    ):
        dtype = "numeric"

    else:
        dtype = "textual"

    unique_nonmissing = len(
        set(nonmissing)
    )

    n_missing = (
        len(values)
        - len(nonmissing)
    )

    return {
        "dtype": dtype,
        "n_missing": n_missing,
        "missing_fraction":
            n_missing / len(values),
        "n_unique_nonmissing":
            unique_nonmissing,
        "constant":
            unique_nonmissing <= 1,
        "all_missing":
            len(nonmissing) == 0,
    }


def chip_family(suffix):

    if suffix in CHIP_PROVENANCE_SUFFIXES:
        return "provenance_metadata"

    if "relative_percentile_" in suffix:
        return "benchmark_relative_spmr"

    if suffix.endswith(
        "_cross_study_range"
    ):
        return "cross_study_variability"

    if suffix.endswith(
        "_consensus_median_across_studies"
    ):
        return "peak_continuous_consensus"

    support_suffixes = (
        "_n_positive_conditions_total",
        "_fraction_positive_conditions_raw",
        "_fraction_positive_conditions_equal_study_weight",
        "_fraction_studies_with_any_condition",
        "_fraction_studies_with_all_conditions",
        "_any_study",
        "_all_studies_have_any_condition",
    )

    if suffix.endswith(
        support_suffixes
    ):
        return "peak_support_summary"

    return "UNCLASSIFIED"


def classify_core(field):

    if field == "canonical_locus_id":
        return (
            "identifier",
            "benchmark_identity",
            "id_only",
            "Stable locus key retained for joins; never a predictor.",
        )

    if field == "benchmark_role":
        return (
            "label",
            "benchmark_outcome",
            "label_only",
            "Benchmark class label; must never be used as an input predictor.",
        )

    if field == "locus_name":
        return (
            "identifier_metadata",
            "benchmark_identity",
            "exclude_metadata",
            "Human-readable locus name; descriptive rather than predictive.",
        )

    if field in {
        "target_seqname",
        "target_start_1based",
        "target_end_1based",
        "locus_width_bp",
    }:
        return (
            "benchmark_coordinate_metadata",
            "benchmark_definition",
            "exclude_coordinate",
            "Part of the curated benchmark locus definition; excluded to avoid coordinate/provenance proxies.",
        )

    if field == "target_assembly":
        return (
            "reference_metadata",
            "reference_genome",
            "exclude_metadata",
            "Reference assembly metadata; constant across this benchmark.",
        )

    if field in {
        "best_evidence_tier",
        "n_observations",
        "locus_semantic_class",
        "mapping_confidence",
        "supporting_studies",
    }:
        return (
            "evidence_metadata",
            "benchmark_evidence",
            "exclude_metadata",
            "Describes benchmark evidence or curation rather than molecular state.",
        )

    fail(
        f"Unclassified CORE field: {field}"
    )


def classify_rna(field):

    if field in RNA_HIGH_CARDINALITY:
        return (
            "biological_annotation",
            "nearest_gene_annotation",
            "exclude_high_cardinality",
            "Gene identifier is retained for interpretation but would act as a high-cardinality locus identifier.",
        )

    if field in RNA_COORDINATE:
        return (
            "biological_annotation",
            "nearest_gene_coordinate",
            "exclude_coordinate",
            "Absolute TSS coordinate is an annotation coordinate, not an expression measurement.",
        )

    if field in RNA_CATEGORICAL_BIOLOGICAL:
        return (
            "biological_context",
            "tss_context",
            "candidate_categorical",
            "Biological locus-to-TSS context; requires explicit categorical/binary encoding before modelling.",
        )

    if field in RNA_NUMERIC_BIOLOGICAL:
        return (
            "biological_context",
            "tss_context",
            "candidate_numeric",
            "Numeric locus-to-TSS context derived independently of the benchmark label.",
        )

    if field in RNA_PROVENANCE:
        return (
            "provenance_metadata",
            "rna_dataset_scope",
            "exclude_metadata",
            "Dataset-scope provenance rather than locus-specific molecular signal.",
        )

    if field in RNA_TECHNICAL:
        return (
            "technical_or_detection_summary",
            "rna_raw_detection",
            "review_technical",
            "Raw-count/detection summary retained for interpretation; not automatically accepted as a portable predictor.",
        )

    if field in RNA_SENSITIVITY:
        return (
            "sensitivity_analysis",
            "rna_PRJEB35183_sensitivity",
            "sensitivity_only",
            "Leave-one-study-out sensitivity metric; retained for robustness analysis rather than baseline prediction.",
        )

    if field in RNA_EXPRESSION:
        return (
            "biological_signal",
            "rna_expression_consensus",
            "candidate_redundancy_review",
            "Expression-context feature; eligible for modelling after redundancy and dimensionality review.",
        )

    fail(
        f"Unclassified RNA field: {field}"
    )


def classify_atac(field):

    if field == "atac_anchor_1based":
        return (
            "benchmark_coordinate_metadata",
            "atac_anchor",
            "exclude_coordinate",
            "Representative benchmark anchor coordinate; not an accessibility measurement.",
        )

    if field.startswith(
        "atac_idr05_"
    ):
        return (
            "sensitivity_analysis",
            "atac_IDR05_sensitivity",
            "sensitivity_only",
            "IDR 0.05 layer retained as a stricter sensitivity analysis; IDR 0.10 is the primary peak layer.",
        )

    if (
        field.startswith(
            "atac_idr10_"
        )
        or field.startswith(
            "atac_spmr_consensus_"
        )
    ):
        return (
            "biological_signal",
            "atac_primary_accessibility",
            "candidate_redundancy_review",
            "Primary accessibility context derived from IDR 0.10 peaks or replicate-consensus SPMR.",
        )

    fail(
        f"Unclassified ATAC field: {field}"
    )


def classify_chip(field):

    if field == "chip_anchor_1based":
        return (
            "benchmark_coordinate_metadata",
            "chip_anchor",
            "exclude_coordinate",
            "Representative benchmark anchor coordinate; not a histone-mark signal.",
            "",
            "",
        )

    if "__" not in field:
        fail(
            f"Unexpected unprefixed ChIP field: {field}"
        )

    mark, suffix = field.split(
        "__",
        1,
    )

    if mark not in MARKS:
        fail(
            f"Unexpected ChIP mark prefix: {mark}"
        )

    family = chip_family(
        suffix
    )

    if family == "UNCLASSIFIED":
        fail(
            f"Unclassified ChIP suffix: {suffix}"
        )

    if family == "provenance_metadata":
        return (
            "provenance_metadata",
            family,
            "exclude_metadata",
            "Mark/study/analysis provenance retained for traceability; not a locus-specific predictor.",
            mark,
            suffix,
        )

    if family == "benchmark_relative_spmr":
        return (
            "benchmark_relative_feature",
            family,
            "exclude_benchmark_relative",
            "Percentile was computed relative to the 37 benchmark loci and is therefore benchmark-composition-dependent.",
            mark,
            suffix,
        )

    if family == "cross_study_variability":
        return (
            "sensitivity_analysis",
            family,
            "sensitivity_only",
            "Between-study variability is retained as robustness/heterogeneity information rather than a baseline biological predictor.",
            mark,
            suffix,
        )

    if family == "peak_support_summary":
        return (
            "biological_signal",
            family,
            "review_redundancy",
            "Peak-support feature is biologically informative but highly related to other support encodings; requires redundancy review.",
            mark,
            suffix,
        )

    if family == "peak_continuous_consensus":
        return (
            "biological_signal",
            family,
            "candidate_redundancy_review",
            "Peak-context consensus is portable to additional loci on the same fixed tracks, subject to redundancy review.",
            mark,
            suffix,
        )

    fail(
        f"Unhandled ChIP family: {family}"
    )


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--master",
        required=True,
    )

    parser.add_argument(
        "--master-provenance",
        required=True,
    )

    parser.add_argument(
        "--manifest",
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
        Path(args.manifest),
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
            + ",".join(existing)
        )

    fields, rows = read_tsv(
        args.master
    )

    if len(rows) != EXPECTED_ROWS:
        fail(
            f"Expected {EXPECTED_ROWS} rows; "
            f"found {len(rows)}"
        )

    if len(fields) != EXPECTED_COLUMNS:
        fail(
            f"Expected {EXPECTED_COLUMNS} columns; "
            f"found {len(fields)}"
        )

    with Path(
        args.master_provenance
    ).open(
        "r",
        encoding="utf-8",
    ) as fh:

        master_prov = json.load(fh)

    blocks = master_prov[
        "feature_blocks"
    ]

    modality_by_field = {}

    for field in CORE_FIELDS:
        modality_by_field[field] = "CORE"

    for modality in (
        "RNA",
        "ATAC",
        "CHIP",
    ):

        for field in blocks[
            modality
        ]:

            if field in modality_by_field:
                fail(
                    "Duplicate modality assignment: "
                    + field
                )

            modality_by_field[
                field
            ] = modality

    if set(modality_by_field) != set(
        fields
    ):
        missing = sorted(
            set(fields)
            - set(modality_by_field)
        )

        extra = sorted(
            set(modality_by_field)
            - set(fields)
        )

        fail(
            "Modality map mismatch. "
            f"missing={missing}; extra={extra}"
        )

    modality_counts = Counter(
        modality_by_field.values()
    )

    if dict(
        modality_counts
    ) != EXPECTED_MODALITY_COUNTS:
        fail(
            "Unexpected modality distribution: "
            + repr(
                dict(modality_counts)
            )
        )

    # Validate RNA rule partition explicitly.
    expected_rna_partition = (
        RNA_HIGH_CARDINALITY
        | RNA_COORDINATE
        | RNA_CATEGORICAL_BIOLOGICAL
        | RNA_NUMERIC_BIOLOGICAL
        | RNA_PROVENANCE
        | RNA_TECHNICAL
        | RNA_SENSITIVITY
        | RNA_EXPRESSION
    )

    if set(
        blocks["RNA"]
    ) != expected_rna_partition:
        fail(
            "RNA classification partition does not "
            "exactly match the frozen RNA feature block"
        )

    # Validate ChIP 510 + 1 structure and 85 suffix definitions.
    chip_fields = blocks[
        "CHIP"
    ]

    mark_prefixes = tuple(
        mark + "__"
        for mark in MARKS
    )

    chip_unprefixed = [
        field
        for field in chip_fields
        if not field.startswith(
            mark_prefixes
        )
    ]

    chip_prefixed = [
        field
        for field in chip_fields
        if field.startswith(
            mark_prefixes
        )
    ]

    if chip_unprefixed != [
        "chip_anchor_1based"
    ]:
        fail(
            "Unexpected ChIP unprefixed fields: "
            + repr(chip_unprefixed)
        )

    if len(chip_prefixed) != 510:
        fail(
            f"Expected 510 mark-prefixed ChIP fields; "
            f"found {len(chip_prefixed)}"
        )

    suffix_maps = defaultdict(
        set
    )

    for field in chip_prefixed:

        mark, suffix = field.split(
            "__",
            1,
        )

        suffix_maps[
            suffix
        ].add(mark)

    if len(suffix_maps) != 85:
        fail(
            f"Expected 85 ChIP suffixes; "
            f"found {len(suffix_maps)}"
        )

    for suffix, marks in suffix_maps.items():

        if marks != set(MARKS):
            fail(
                "Incomplete mark replication for "
                + suffix
            )

    family_counts = Counter(
        chip_family(suffix)
        for suffix in suffix_maps
    )

    expected_family_counts = Counter({
        "provenance_metadata": 6,
        "peak_support_summary": 35,
        "peak_continuous_consensus": 14,
        "cross_study_variability": 14,
        "benchmark_relative_spmr": 16,
    })

    if family_counts != expected_family_counts:
        fail(
            "Unexpected ChIP family counts: "
            + repr(
                dict(family_counts)
            )
        )

    manifest_rows = []

    for index, field in enumerate(
        fields,
        start=1,
    ):

        modality = modality_by_field[
            field
        ]

        p = profile_field(
            field,
            rows,
        )

        mark = ""
        suffix = ""

        if modality == "CORE":

            (
                category,
                family,
                policy,
                reason,
            ) = classify_core(
                field
            )

        elif modality == "RNA":

            (
                category,
                family,
                policy,
                reason,
            ) = classify_rna(
                field
            )

        elif modality == "ATAC":

            (
                category,
                family,
                policy,
                reason,
            ) = classify_atac(
                field
            )

        elif modality == "CHIP":

            (
                category,
                family,
                policy,
                reason,
                mark,
                suffix,
            ) = classify_chip(
                field
            )

        else:
            fail(
                "Unexpected modality: "
                + modality
            )

        benchmark_relative = (
            "relative_percentile_"
            in field
        )

        # Apply data-driven guards only to categories that
        # could otherwise enter or be reviewed for modelling.
        if p["all_missing"]:

            if policy in {
                "candidate_numeric",
                "candidate_categorical",
                "candidate_redundancy_review",
                "review_redundancy",
                "review_technical",
                "sensitivity_only",
            }:

                policy = (
                    "exclude_not_applicable"
                )

                reason = (
                    reason
                    + " This column is all-missing "
                    "across the 37 loci."
                )

        elif p["constant"]:

            if policy in {
                "candidate_numeric",
                "candidate_categorical",
                "candidate_redundancy_review",
                "review_redundancy",
                "review_technical",
            }:

                policy = "exclude_constant"

                reason = (
                    reason
                    + " This column is constant "
                    "across the 37 benchmark loci."
                )

        requires_encoding = (
            p["dtype"] == "textual"
            and policy
            in {
                "candidate_categorical",
                "candidate_redundancy_review",
                "review_redundancy",
                "review_technical",
            }
        )

        manifest_rows.append({
            "column_index":
                str(index),

            "column_name":
                field,

            "modality":
                modality,

            "histone_mark":
                mark,

            "feature_suffix":
                suffix,

            "semantic_category":
                category,

            "feature_family":
                family,

            "dtype":
                p["dtype"],

            "n_missing":
                str(
                    p["n_missing"]
                ),

            "missing_fraction":
                f"{p['missing_fraction']:.12g}",

            "n_unique_nonmissing":
                str(
                    p[
                        "n_unique_nonmissing"
                    ]
                ),

            "constant":
                str(
                    p["constant"]
                ).upper(),

            "all_missing":
                str(
                    p["all_missing"]
                ).upper(),

            "benchmark_relative":
                str(
                    benchmark_relative
                ).upper(),

            "requires_encoding":
                str(
                    requires_encoding
                ).upper(),

            "baseline_model_policy":
                policy,

            "rationale":
                reason,
        })

    # --------------------------------------------------------
    # Global independent validations of manifest content.
    # --------------------------------------------------------

    if len(manifest_rows) != EXPECTED_COLUMNS:
        fail(
            "Manifest row count mismatch"
        )

    if len({
        row["column_name"]
        for row in manifest_rows
    }) != EXPECTED_COLUMNS:
        fail(
            "Manifest column names are not unique"
        )

    if [
        row["column_name"]
        for row in manifest_rows
    ] != fields:
        fail(
            "Manifest does not preserve master column order"
        )

    constant_count = sum(
        row["constant"] == "TRUE"
        for row in manifest_rows
    )

    all_missing_count = sum(
        row["all_missing"] == "TRUE"
        for row in manifest_rows
    )

    columns_with_missing = sum(
        int(
            row["n_missing"]
        ) > 0
        for row in manifest_rows
    )

    benchmark_relative_count = sum(
        row[
            "benchmark_relative"
        ] == "TRUE"
        for row in manifest_rows
    )

    if (
        constant_count
        != EXPECTED_CONSTANT_COLUMNS
    ):
        fail(
            f"Expected {EXPECTED_CONSTANT_COLUMNS} "
            f"constant columns; found "
            f"{constant_count}"
        )

    if (
        all_missing_count
        != EXPECTED_ALL_MISSING_COLUMNS
    ):
        fail(
            f"Expected {EXPECTED_ALL_MISSING_COLUMNS} "
            f"all-missing columns; found "
            f"{all_missing_count}"
        )

    if (
        columns_with_missing
        != EXPECTED_COLUMNS_WITH_MISSING
    ):
        fail(
            f"Expected {EXPECTED_COLUMNS_WITH_MISSING} "
            f"columns with missingness; found "
            f"{columns_with_missing}"
        )

    if (
        benchmark_relative_count
        != EXPECTED_BENCHMARK_RELATIVE
    ):
        fail(
            f"Expected {EXPECTED_BENCHMARK_RELATIVE} "
            "benchmark-relative columns; found "
            f"{benchmark_relative_count}"
        )

    # Benchmark-relative columns must all be ChIP and excluded.
    for row in manifest_rows:

        if (
            row[
                "benchmark_relative"
            ] == "TRUE"
        ):

            if row["modality"] != "CHIP":
                fail(
                    "Non-ChIP benchmark-relative column: "
                    + row["column_name"]
                )

            if (
                row[
                    "baseline_model_policy"
                ]
                != "exclude_benchmark_relative"
            ):
                fail(
                    "Benchmark-relative column not "
                    "excluded: "
                    + row["column_name"]
                )

    if next(
        row
        for row in manifest_rows
        if row["column_name"]
        == "benchmark_role"
    )["baseline_model_policy"] != "label_only":
        fail(
            "benchmark_role policy is not label_only"
        )

    for coordinate in (
        "target_seqname",
        "target_start_1based",
        "target_end_1based",
        "locus_width_bp",
        "rna_nearest_gene_tss_1based",
        "atac_anchor_1based",
        "chip_anchor_1based",
    ):

        row = next(
            item
            for item in manifest_rows
            if item["column_name"]
            == coordinate
        )

        if row[
            "baseline_model_policy"
        ] != "exclude_coordinate":
            fail(
                "Coordinate not excluded: "
                + coordinate
            )

    # support_only is explicitly not treated as negative.
    role_counts = Counter(
        row["benchmark_role"]
        for row in rows
    )

    if role_counts != Counter({
        "positive": 27,
        "negative": 4,
        "support_only": 6,
    }):
        fail(
            "Unexpected benchmark role distribution"
        )

    policy_counts = Counter(
        row[
            "baseline_model_policy"
        ]
        for row in manifest_rows
    )

    category_counts = Counter(
        row[
            "semantic_category"
        ]
        for row in manifest_rows
    )

    modality_manifest_counts = Counter(
        row["modality"]
        for row in manifest_rows
    )

    # --------------------------------------------------------
    # Write manifest.
    # --------------------------------------------------------

    manifest_fields = [
        "column_index",
        "column_name",
        "modality",
        "histone_mark",
        "feature_suffix",
        "semantic_category",
        "feature_family",
        "dtype",
        "n_missing",
        "missing_fraction",
        "n_unique_nonmissing",
        "constant",
        "all_missing",
        "benchmark_relative",
        "requires_encoding",
        "baseline_model_policy",
        "rationale",
    ]

    write_atomic_tsv(
        args.manifest,
        manifest_fields,
        manifest_rows,
    )

    # --------------------------------------------------------
    # QC summary.
    # --------------------------------------------------------

    qc_rows = [
        {
            "metric":
                "master_rows",
            "value":
                str(len(rows)),
        },
        {
            "metric":
                "master_columns",
            "value":
                str(len(fields)),
        },
        {
            "metric":
                "manifest_rows",
            "value":
                str(
                    len(manifest_rows)
                ),
        },
        {
            "metric":
                "unique_manifest_columns",
            "value":
                str(
                    len({
                        row["column_name"]
                        for row
                        in manifest_rows
                    })
                ),
        },
        {
            "metric":
                "constant_columns",
            "value":
                str(constant_count),
        },
        {
            "metric":
                "all_missing_columns",
            "value":
                str(
                    all_missing_count
                ),
        },
        {
            "metric":
                "columns_with_missing",
            "value":
                str(
                    columns_with_missing
                ),
        },
        {
            "metric":
                "benchmark_relative_columns",
            "value":
                str(
                    benchmark_relative_count
                ),
        },
        {
            "metric":
                "benchmark_positive",
            "value":
                str(
                    role_counts["positive"]
                ),
        },
        {
            "metric":
                "benchmark_negative",
            "value":
                str(
                    role_counts["negative"]
                ),
        },
        {
            "metric":
                "benchmark_support_only",
            "value":
                str(
                    role_counts[
                        "support_only"
                    ]
                ),
        },
    ]

    for modality in (
        "CORE",
        "RNA",
        "ATAC",
        "CHIP",
    ):

        qc_rows.append({
            "metric":
                "modality::"
                + modality,
            "value":
                str(
                    modality_manifest_counts[
                        modality
                    ]
                ),
        })

    for policy in sorted(
        policy_counts
    ):

        qc_rows.append({
            "metric":
                "policy::"
                + policy,
            "value":
                str(
                    policy_counts[
                        policy
                    ]
                ),
        })

    for category in sorted(
        category_counts
    ):

        qc_rows.append({
            "metric":
                "category::"
                + category,
            "value":
                str(
                    category_counts[
                        category
                    ]
                ),
        })

    qc_rows.extend([
        {
            "metric":
                "chip_unique_feature_definitions",
            "value":
                "85",
        },
        {
            "metric":
                "chip_provenance_definitions",
            "value":
                "6",
        },
        {
            "metric":
                "chip_peak_support_definitions",
            "value":
                "35",
        },
        {
            "metric":
                "chip_peak_continuous_definitions",
            "value":
                "14",
        },
        {
            "metric":
                "chip_cross_study_variability_definitions",
            "value":
                "14",
        },
        {
            "metric":
                "chip_benchmark_relative_definitions",
            "value":
                "16",
        },
        {
            "metric":
                "validation_errors",
            "value":
                "0",
        },
        {
            "metric":
                "status",
            "value":
                "PASS",
        },
    ])

    write_atomic_tsv(
        args.qc,
        [
            "metric",
            "value",
        ],
        qc_rows,
    )

    # --------------------------------------------------------
    # Provenance.
    # --------------------------------------------------------

    provenance = {
        "workflow":
            "Multiomic feature governance manifest",

        "source_master":
            str(
                Path(args.master)
            ),

        "source_master_sha256":
            sha256_file(
                args.master
            ),

        "source_master_provenance":
            str(
                Path(
                    args.master_provenance
                )
            ),

        "source_master_provenance_sha256":
            sha256_file(
                args.master_provenance
            ),

        "producer_script":
            str(
                Path(sys.argv[0])
            ),

        "producer_script_sha256":
            sha256_file(
                sys.argv[0]
            ),

        "manifest_shape": {
            "rows":
                len(manifest_rows),
            "columns":
                len(
                    manifest_fields
                ),
        },

        "master_shape": {
            "rows":
                len(rows),
            "columns":
                len(fields),
        },

        "classification_scope": (
            "This manifest assigns semantic and "
            "governance roles to every column in the "
            "frozen 37-locus multiomic master. It does "
            "not create a training matrix and does not "
            "perform model fitting or outcome-driven "
            "feature selection."
        ),

        "label_policy": (
            "benchmark_role is retained as the "
            "benchmark label only. It must never be "
            "included among predictors. support_only "
            "records are not negatives."
        ),

        "coordinate_policy": (
            "Absolute benchmark coordinates, locus "
            "width, RNA TSS coordinates, and ATAC/ChIP "
            "representative anchors are excluded from "
            "the baseline predictor space to avoid "
            "coordinate and benchmark-definition "
            "proxies."
        ),

        "rna_policy": (
            "Primary normalized expression and TSS "
            "context are retained as candidate "
            "biological features. Raw-count detection "
            "summaries require technical review. "
            "Features recalculated without PRJEB35183 "
            "are retained as sensitivity-analysis "
            "variables rather than baseline predictors."
        ),

        "atac_policy": (
            "IDR 0.10 and replicate-consensus SPMR "
            "features define the primary accessibility "
            "layer. IDR 0.05 features are retained as "
            "a stricter sensitivity-analysis layer."
        ),

        "chip_policy": (
            "Study/analysis provenance is excluded. "
            "Peak-support and peak-context features "
            "remain biologically interpretable "
            "candidates subject to redundancy review. "
            "Cross-study variability is treated as "
            "sensitivity information."
        ),

        "chip_benchmark_relative_policy": (
            "All 96 ChIP relative_percentile columns "
            "remain excluded from a portable baseline "
            "predictor view because their values depend "
            "on the distribution of the same 37 "
            "benchmark loci."
        ),

        "constant_policy": (
            "Biological or technical features that are "
            "constant across the 37 benchmark loci are "
            "marked exclude_constant for the baseline "
            "benchmark feature view. Their source "
            "values remain preserved in the frozen "
            "master."
        ),

        "missingness_policy": (
            "All-missing cross-study fields that are "
            "not applicable to one-study histone marks "
            "are marked exclude_not_applicable. No "
            "imputation is performed by this manifest."
        ),

        "selection_warning": (
            "Candidate and review policies are "
            "pre-model governance labels, not evidence "
            "that a feature should be used in a final "
            "model. Redundancy, dimensionality, "
            "validation design, and the very small "
            "negative class must be addressed "
            "separately."
        ),

        "benchmark_composition": {
            "positive": 27,
            "negative": 4,
            "support_only": 6,
        },

        "observed_profile": {
            "constant_columns":
                constant_count,
            "all_missing_columns":
                all_missing_count,
            "columns_with_missing":
                columns_with_missing,
            "benchmark_relative_columns":
                benchmark_relative_count,
        },

        "policy_counts":
            dict(
                sorted(
                    policy_counts.items()
                )
            ),

        "semantic_category_counts":
            dict(
                sorted(
                    category_counts.items()
                )
            ),

        "validation": {
            "errors": 0,
            "status": "PASS",
        },
    }

    write_atomic_json(
        args.provenance,
        provenance,
    )

    sha_tmp = Path(
        args.sha256s
    ).with_name(
        Path(args.sha256s).name
        + ".tmp"
    )

    with sha_tmp.open(
        "w",
        encoding="utf-8",
    ) as fh:

        for path in (
            Path(args.manifest),
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
        "MULTIOMIC_FEATURE_MANIFEST_BUILD=PASS"
    )

    print(
        f"master_rows={len(rows)}"
    )

    print(
        f"master_columns={len(fields)}"
    )

    print(
        f"manifest_rows={len(manifest_rows)}"
    )

    print(
        f"constant_columns={constant_count}"
    )

    print(
        f"all_missing_columns={all_missing_count}"
    )

    print(
        f"columns_with_missing={columns_with_missing}"
    )

    print(
        "benchmark_relative_columns="
        f"{benchmark_relative_count}"
    )

    print()

    print("POLICY_COUNTS")

    for policy in sorted(
        policy_counts
    ):

        print(
            f"{policy}="
            f"{policy_counts[policy]}"
        )

    print()

    print("SEMANTIC_CATEGORY_COUNTS")

    for category in sorted(
        category_counts
    ):

        print(
            f"{category}="
            f"{category_counts[category]}"
        )

    print()

    print(
        "validation_errors=0"
    )

    print(
        "status=PASS"
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
