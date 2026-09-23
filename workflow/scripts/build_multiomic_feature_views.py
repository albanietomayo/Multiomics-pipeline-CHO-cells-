#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import os
import sys
from collections import Counter
from pathlib import Path


EXPECTED_MASTER_ROWS = 37
EXPECTED_MASTER_COLUMNS = 588
EXPECTED_MANIFEST_ROWS = 588

EXPECTED_CANDIDATES = 317
EXPECTED_DETERMINISTIC_REMOVALS = 150

EXPECTED_EXTENDED = 167
EXPECTED_EXTENDED_MODALITIES = Counter({
    "RNA": 10,
    "ATAC": 12,
    "CHIP": 145,
})

EXPECTED_CORE = 67
EXPECTED_CORE_MODALITIES = Counter({
    "RNA": 5,
    "ATAC": 8,
    "CHIP": 54,
})


CANDIDATE_POLICIES = {
    "candidate_categorical",
    "candidate_numeric",
    "candidate_redundancy_review",
    "review_redundancy",
}


MARKS = [
    "H3K27ac",
    "H3K27me3",
    "H3K36me3",
    "H3K4me1",
    "H3K4me3",
    "H3K9me3",
]


CORE_RNA = [
    "rna_locus_relative_to_tss",
    "rna_min_tss_distance_bp",
    "rna_tmm_logcpm_median",
    "rna_tpm_median",
    "rna_study_median_tpm_positive_fraction",
]


CORE_ATAC = [
    "atac_idr10_exact_n_unique_peaks",
    "atac_idr10_exact_overlap_fraction_of_locus",

    "atac_idr10_accessible_bp_pm1000bp",
    "atac_spmr_consensus_pm1000bp_mean",

    "atac_idr10_accessible_bp_pm10000bp",
    "atac_spmr_consensus_pm10000bp_mean",

    "atac_idr10_accessible_bp_pm50000bp",
    "atac_spmr_consensus_pm50000bp_mean",
]


CORE_CHIP_SUFFIXES = [
    "chip_exact_overlap_fraction_positive_conditions_equal_study_weight",

    "chip_any_peak_pm1000bp_fraction_positive_conditions_equal_study_weight",
    "chip_any_peak_pm10000bp_fraction_positive_conditions_equal_study_weight",
    "chip_any_peak_pm50000bp_fraction_positive_conditions_equal_study_weight",

    "chip_exact_overlap_fraction_of_locus_consensus_median_across_studies",

    "chip_nearest_peak_distance_to_locus_bp_consensus_median_across_studies",

    "chip_enriched_fraction_pm1000bp_consensus_median_across_studies",
    "chip_enriched_fraction_pm10000bp_consensus_median_across_studies",
    "chip_enriched_fraction_pm50000bp_consensus_median_across_studies",
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


def derive_extended(
    master_fields,
    manifest_rows,
):

    manifest_by_name = {
        row["column_name"]: row
        for row in manifest_rows
    }

    candidates = {
        row["column_name"]: row
        for row in manifest_rows
        if row[
            "baseline_model_policy"
        ] in CANDIDATE_POLICIES
    }

    if len(candidates) != EXPECTED_CANDIDATES:
        fail(
            f"Expected {EXPECTED_CANDIDATES} candidates; "
            f"found {len(candidates)}"
        )

    remove = set()

    # --------------------------------------------------------
    # RNA deterministic redundancy
    # --------------------------------------------------------

    remove.update({
        "rna_tmm_logcpm_iqr",
        "rna_tpm_iqr",
        "rna_study_median_tpm_positive_n",
        "rna_tss_inside_locus",
        "rna_tss_inside_count",
    })

    # --------------------------------------------------------
    # ATAC deterministic redundancy
    # --------------------------------------------------------

    remove.update({
        "atac_idr10_exact_overlap",
        "atac_idr10_exact_overlap_bp",
        "atac_idr10_any_peak_pm1000bp",
        "atac_idr10_any_peak_pm10000bp",
        "atac_idr10_any_peak_pm50000bp",
    })

    # --------------------------------------------------------
    # ChIP deterministic redundancy
    # --------------------------------------------------------

    for name, row in candidates.items():

        if row["modality"] != "CHIP":
            continue

        suffix = row[
            "feature_suffix"
        ]

        if suffix.endswith(
            "_n_positive_conditions_total"
        ):
            remove.add(name)

        if suffix.endswith(
            "_fraction_positive_conditions_raw"
        ):
            remove.add(name)

        if suffix.endswith(
            "_any_study"
        ):
            remove.add(name)

        if suffix.endswith(
            "_all_studies_have_any_condition"
        ):
            remove.add(name)

        if (
            "chip_exact_overlap_bp_"
            in suffix
            and suffix.endswith(
                "_consensus_median_across_studies"
            )
        ):
            remove.add(name)

        if (
            "chip_enriched_bp_"
            in suffix
            and suffix.endswith(
                "_consensus_median_across_studies"
            )
        ):
            remove.add(name)

    applicable_remove = (
        set(candidates)
        & remove
    )

    if (
        len(applicable_remove)
        != EXPECTED_DETERMINISTIC_REMOVALS
    ):
        fail(
            "Expected "
            f"{EXPECTED_DETERMINISTIC_REMOVALS} "
            "deterministic removals; found "
            f"{len(applicable_remove)}"
        )

    extended = [
        field
        for field in master_fields
        if (
            field in candidates
            and field not in applicable_remove
        )
    ]

    if len(extended) != EXPECTED_EXTENDED:
        fail(
            f"Expected {EXPECTED_EXTENDED} "
            f"extended predictors; found "
            f"{len(extended)}"
        )

    modality_counts = Counter(
        manifest_by_name[field][
            "modality"
        ]
        for field in extended
    )

    if (
        modality_counts
        != EXPECTED_EXTENDED_MODALITIES
    ):
        fail(
            "Unexpected extended modality counts: "
            + repr(
                dict(modality_counts)
            )
        )

    return (
        candidates,
        applicable_remove,
        extended,
    )


def derive_core():

    core_chip = []

    for mark in MARKS:

        for suffix in (
            CORE_CHIP_SUFFIXES
        ):

            core_chip.append(
                f"{mark}__{suffix}"
            )

    core = (
        CORE_RNA
        + CORE_ATAC
        + core_chip
    )

    if len(core) != EXPECTED_CORE:
        fail(
            f"Expected {EXPECTED_CORE} core predictors; "
            f"found {len(core)}"
        )

    return core


def validate_view(
    label,
    predictors,
    master_fields,
    manifest_by_name,
):

    if len(set(predictors)) != len(
        predictors
    ):
        fail(
            f"{label}: duplicate predictors"
        )

    missing = [
        field
        for field in predictors
        if field not in master_fields
    ]

    if missing:
        fail(
            f"{label}: predictors missing from master: "
            + ",".join(missing)
        )

    benchmark_relative = []
    constants = []
    all_missing = []
    with_missing = []
    sensitivity = []
    coordinates = []
    wrong_policy = []

    for field in predictors:

        row = manifest_by_name[
            field
        ]

        if row[
            "benchmark_relative"
        ] == "TRUE":
            benchmark_relative.append(
                field
            )

        if row["constant"] == "TRUE":
            constants.append(
                field
            )

        if row["all_missing"] == "TRUE":
            all_missing.append(
                field
            )

        if int(
            row["n_missing"]
        ) > 0:
            with_missing.append(
                field
            )

        if row[
            "baseline_model_policy"
        ] == "sensitivity_only":
            sensitivity.append(
                field
            )

        if row[
            "baseline_model_policy"
        ] == "exclude_coordinate":
            coordinates.append(
                field
            )

        if row[
            "baseline_model_policy"
        ] not in CANDIDATE_POLICIES:
            wrong_policy.append(
                field
            )

    violations = {
        "benchmark_relative":
            benchmark_relative,

        "constant":
            constants,

        "all_missing":
            all_missing,

        "with_missing":
            with_missing,

        "sensitivity":
            sensitivity,

        "coordinate":
            coordinates,

        "noncandidate_policy":
            wrong_policy,
    }

    nonempty = {
        key: values
        for key, values
        in violations.items()
        if values
    }

    if nonempty:
        fail(
            f"{label}: governance violations: "
            + repr(nonempty)
        )


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--master",
        required=True,
    )

    parser.add_argument(
        "--manifest",
        required=True,
    )

    parser.add_argument(
        "--core",
        required=True,
    )

    parser.add_argument(
        "--extended",
        required=True,
    )

    parser.add_argument(
        "--spec",
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
        Path(args.core),
        Path(args.extended),
        Path(args.spec),
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

    (
        master_fields,
        master_rows,
    ) = read_tsv(
        args.master
    )

    (
        manifest_fields,
        manifest_rows,
    ) = read_tsv(
        args.manifest
    )

    if len(
        master_rows
    ) != EXPECTED_MASTER_ROWS:
        fail(
            "Unexpected master row count"
        )

    if len(
        master_fields
    ) != EXPECTED_MASTER_COLUMNS:
        fail(
            "Unexpected master column count"
        )

    if len(
        manifest_rows
    ) != EXPECTED_MANIFEST_ROWS:
        fail(
            "Unexpected manifest row count"
        )

    manifest_by_name = {
        row["column_name"]: row
        for row in manifest_rows
    }

    if len(
        manifest_by_name
    ) != EXPECTED_MANIFEST_ROWS:
        fail(
            "Manifest column names are not unique"
        )

    if list(
        manifest_by_name
    ) != master_fields:
        fail(
            "Manifest does not preserve master column order"
        )

    (
        candidates,
        deterministic_removals,
        extended,
    ) = derive_extended(
        master_fields,
        manifest_rows,
    )

    core = derive_core()

    if not set(core).issubset(
        set(extended)
    ):
        fail(
            "Core is not a subset of extended"
        )

    validate_view(
        "core",
        core,
        master_fields,
        manifest_by_name,
    )

    validate_view(
        "extended",
        extended,
        master_fields,
        manifest_by_name,
    )

    core_modality_counts = Counter(
        manifest_by_name[field][
            "modality"
        ]
        for field in core
    )

    if (
        core_modality_counts
        != EXPECTED_CORE_MODALITIES
    ):
        fail(
            "Unexpected core modality counts: "
            + repr(
                dict(core_modality_counts)
            )
        )

    # --------------------------------------------------------
    # Explicit ID / label safety.
    # --------------------------------------------------------

    for forbidden in (
        "canonical_locus_id",
        "benchmark_role",
        "locus_name",
        "target_seqname",
        "target_start_1based",
        "target_end_1based",
        "locus_width_bp",
    ):

        if forbidden in core:
            fail(
                "Forbidden field in core predictors: "
                + forbidden
            )

        if forbidden in extended:
            fail(
                "Forbidden field in extended predictors: "
                + forbidden
            )

    if (
        "atac_idr10_nearest_peak_distance_bp"
        in core
    ):
        fail(
            "Anchor-based ATAC nearest distance "
            "must not be in core"
        )

    if (
        "atac_idr10_nearest_peak_distance_bp"
        not in extended
    ):
        fail(
            "Expected anchor-based ATAC nearest "
            "distance in extended"
        )

    chip_anchor_distance_core = [
        field
        for field in core
        if "nearest_peak_distance_to_anchor"
        in field
    ]

    if chip_anchor_distance_core:
        fail(
            "Anchor-based ChIP distance found in core"
        )

    chip_locus_distance_core = [
        field
        for field in core
        if "nearest_peak_distance_to_locus"
        in field
    ]

    if len(
        chip_locus_distance_core
    ) != 6:
        fail(
            "Expected six full-locus ChIP distances "
            "in core"
        )

    # --------------------------------------------------------
    # Build tables.
    # --------------------------------------------------------

    id_field = "canonical_locus_id"
    label_field = "benchmark_role"

    core_fields = [
        id_field,
        label_field,
    ] + core

    extended_fields = [
        id_field,
        label_field,
    ] + extended

    core_rows = [
        {
            field: row[field]
            for field in core_fields
        }
        for row in master_rows
    ]

    extended_rows = [
        {
            field: row[field]
            for field
            in extended_fields
        }
        for row in master_rows
    ]

    # No missing predictor values.
    for view_name, rows, predictors in (
        (
            "core",
            core_rows,
            core,
        ),
        (
            "extended",
            extended_rows,
            extended,
        ),
    ):

        missing_cells = []

        for row_index, row in enumerate(
            rows,
            start=1,
        ):

            for field in predictors:

                if row[field] in (
                    "",
                    "NA",
                ):

                    missing_cells.append(
                        (
                            row_index,
                            field,
                        )
                    )

        if missing_cells:
            fail(
                f"{view_name}: missing predictor cells: "
                + repr(
                    missing_cells[:20]
                )
            )

    write_atomic_tsv(
        args.core,
        core_fields,
        core_rows,
    )

    write_atomic_tsv(
        args.extended,
        extended_fields,
        extended_rows,
    )

    # --------------------------------------------------------
    # Build 588-row view specification.
    # --------------------------------------------------------

    core_order = {
        field: index
        for index, field
        in enumerate(
            core,
            start=1,
        )
    }

    extended_order = {
        field: index
        for index, field
        in enumerate(
            extended,
            start=1,
        )
    }

    spec_rows = []

    for row in manifest_rows:

        field = row[
            "column_name"
        ]

        baseline_policy = row[
            "baseline_model_policy"
        ]

        is_candidate = (
            field in candidates
        )

        is_removed = (
            field
            in deterministic_removals
        )

        in_extended = (
            field in extended_order
        )

        in_core = (
            field in core_order
        )

        if in_extended:
            extended_status = (
                "included_predictor"
            )

        elif is_removed:
            extended_status = (
                "excluded_deterministic_redundancy"
            )

        elif baseline_policy == "label_only":
            extended_status = (
                "label_only"
            )

        elif baseline_policy == "id_only":
            extended_status = (
                "identifier_only"
            )

        else:
            extended_status = (
                "excluded_by_feature_governance"
            )

        if in_core:
            core_status = (
                "included_predictor"
            )

        elif in_extended:
            core_status = (
                "extended_only"
            )

        elif baseline_policy == "label_only":
            core_status = (
                "label_only"
            )

        elif baseline_policy == "id_only":
            core_status = (
                "identifier_only"
            )

        elif is_removed:
            core_status = (
                "excluded_deterministic_redundancy"
            )

        else:
            core_status = (
                "excluded_by_feature_governance"
            )

        spec_rows.append({
            "master_column_index":
                row["column_index"],

            "column_name":
                field,

            "modality":
                row["modality"],

            "semantic_category":
                row["semantic_category"],

            "feature_family":
                row["feature_family"],

            "baseline_model_policy":
                baseline_policy,

            "candidate_before_structural_reduction":
                str(
                    is_candidate
                ).upper(),

            "deterministic_redundancy_removed":
                str(
                    is_removed
                ).upper(),

            "extended_predictor":
                str(
                    in_extended
                ).upper(),

            "extended_predictor_order":
                (
                    str(
                        extended_order[
                            field
                        ]
                    )
                    if in_extended
                    else ""
                ),

            "extended_status":
                extended_status,

            "core_predictor":
                str(
                    in_core
                ).upper(),

            "core_predictor_order":
                (
                    str(
                        core_order[
                            field
                        ]
                    )
                    if in_core
                    else ""
                ),

            "core_status":
                core_status,
        })

    spec_fields = [
        "master_column_index",
        "column_name",
        "modality",
        "semantic_category",
        "feature_family",
        "baseline_model_policy",
        "candidate_before_structural_reduction",
        "deterministic_redundancy_removed",
        "extended_predictor",
        "extended_predictor_order",
        "extended_status",
        "core_predictor",
        "core_predictor_order",
        "core_status",
    ]

    if len(spec_rows) != 588:
        fail(
            "Unexpected spec row count"
        )

    write_atomic_tsv(
        args.spec,
        spec_fields,
        spec_rows,
    )

    # --------------------------------------------------------
    # QC summary.
    # --------------------------------------------------------

    role_counts = Counter(
        row[label_field]
        for row in master_rows
    )

    if role_counts != Counter({
        "positive": 27,
        "negative": 4,
        "support_only": 6,
    }):
        fail(
            "Unexpected benchmark-role distribution"
        )

    core_textual = [
        field
        for field in core
        if manifest_by_name[
            field
        ]["dtype"] == "textual"
    ]

    extended_textual = [
        field
        for field in extended
        if manifest_by_name[
            field
        ]["dtype"] == "textual"
    ]

    if core_textual != [
        "rna_locus_relative_to_tss"
    ]:
        fail(
            "Unexpected core textual predictor set"
        )

    if set(
        extended_textual
    ) != {
        "rna_nearest_gene_strand",
        "rna_locus_relative_to_tss",
    }:
        fail(
            "Unexpected extended textual predictor set"
        )

    qc_rows = [
        {
            "metric": "master_rows",
            "value": "37",
        },
        {
            "metric": "master_columns",
            "value": "588",
        },
        {
            "metric": "starting_candidate_predictors",
            "value": "317",
        },
        {
            "metric": "deterministic_redundancy_removals",
            "value": "150",
        },
        {
            "metric": "extended_predictors",
            "value": "167",
        },
        {
            "metric": "extended_table_columns",
            "value": "169",
        },
        {
            "metric": "extended_RNA_predictors",
            "value": "10",
        },
        {
            "metric": "extended_ATAC_predictors",
            "value": "12",
        },
        {
            "metric": "extended_CHIP_predictors",
            "value": "145",
        },
        {
            "metric": "core_predictors",
            "value": "67",
        },
        {
            "metric": "core_table_columns",
            "value": "69",
        },
        {
            "metric": "core_RNA_predictors",
            "value": "5",
        },
        {
            "metric": "core_ATAC_predictors",
            "value": "8",
        },
        {
            "metric": "core_CHIP_predictors",
            "value": "54",
        },
        {
            "metric": "core_subset_of_extended",
            "value": "TRUE",
        },
        {
            "metric": "core_textual_predictors",
            "value": str(
                len(core_textual)
            ),
        },
        {
            "metric": "extended_textual_predictors",
            "value": str(
                len(extended_textual)
            ),
        },
        {
            "metric": "benchmark_relative_predictors",
            "value": "0",
        },
        {
            "metric": "sensitivity_only_predictors",
            "value": "0",
        },
        {
            "metric": "coordinate_predictors",
            "value": "0",
        },
        {
            "metric": "constant_predictors",
            "value": "0",
        },
        {
            "metric": "predictor_columns_with_missingness",
            "value": "0",
        },
        {
            "metric": "benchmark_positive",
            "value": "27",
        },
        {
            "metric": "benchmark_negative",
            "value": "4",
        },
        {
            "metric": "benchmark_support_only",
            "value": "6",
        },
        {
            "metric": "label_used_for_feature_selection",
            "value": "FALSE",
        },
        {
            "metric": "outcome_driven_selection",
            "value": "FALSE",
        },
        {
            "metric": "validation_errors",
            "value": "0",
        },
        {
            "metric": "status",
            "value": "PASS",
        },
    ]

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
            "Multiomic core and extended feature views",

        "source_master":
            str(
                Path(args.master)
            ),

        "source_master_sha256":
            sha256_file(
                args.master
            ),

        "source_feature_manifest":
            str(
                Path(args.manifest)
            ),

        "source_feature_manifest_sha256":
            sha256_file(
                args.manifest
            ),

        "producer_script":
            str(
                Path(sys.argv[0])
            ),

        "producer_script_sha256":
            sha256_file(
                sys.argv[0]
            ),

        "table_roles": {
            "canonical_locus_id":
                "identifier_only",

            "benchmark_role":
                "label_only",

            "predictors":
                "all remaining columns in each feature-view table",
        },

        "core_view": {
            "rows": 37,
            "table_columns": 69,
            "predictors": 67,
            "RNA_predictors": 5,
            "ATAC_predictors": 8,
            "CHIP_predictors": 54,
            "textual_predictors": core_textual,
        },

        "extended_view": {
            "rows": 37,
            "table_columns": 169,
            "predictors": 167,
            "RNA_predictors": 10,
            "ATAC_predictors": 12,
            "CHIP_predictors": 145,
            "textual_predictors": extended_textual,
        },

        "selection_principle": (
            "Feature membership was pre-specified from "
            "biological semantics, feature governance, "
            "and deterministic/algebraic redundancy. "
            "benchmark_role was not used to select or "
            "rank predictors."
        ),

        "structural_reduction": (
            "The extended view starts from the 317 "
            "governed biological candidate/review "
            "columns and removes 150 representations "
            "demonstrated to be deterministic or "
            "algebraically redundant. No outcome "
            "association was used."
        ),

        "core_principle": (
            "The core view is a compact subset of the "
            "extended view emphasizing central RNA "
            "expression/TSS context, primary IDR<=0.10 "
            "ATAC accessibility, and full-locus or "
            "regional ChIP peak context across all six "
            "histone marks."
        ),

        "anchor_policy": (
            "Point-anchor-specific variables are "
            "excluded from the core view. ATAC nearest "
            "peak distance is anchor-based and is "
            "retained only in the extended view. ChIP "
            "core nearest-peak distances use the full "
            "canonical locus interval. Fixed +/-1 kb, "
            "+/-10 kb and +/-50 kb regional summaries "
            "remain anchor-centred context features."
        ),

        "benchmark_relative_policy": (
            "No ChIP relative-percentile feature is "
            "included because those values depend on "
            "the composition of the same 37 benchmark "
            "loci."
        ),

        "sensitivity_policy": (
            "RNA leave-one-study-out sensitivity "
            "features, ATAC IDR<=0.05 features, and "
            "ChIP cross-study variability features "
            "are not included in either baseline "
            "feature view."
        ),

        "missingness_policy": (
            "No predictor column with missing values "
            "is included. No imputation is performed."
        ),

        "encoding_policy": (
            "Categorical biological predictors are "
            "preserved as raw categorical values in "
            "these feature views. No one-hot encoding "
            "or other model-specific transformation "
            "is applied at this stage."
        ),

        "scaling_policy": (
            "No scaling, centering, normalization, "
            "or feature transformation is fitted at "
            "this stage. Any model-specific transform "
            "must be fitted only within the training "
            "data of the eventual modelling design."
        ),

        "label_policy": (
            "benchmark_role is provided only as a "
            "separate label column. support_only is "
            "preserved as a distinct benchmark class "
            "and must not be silently recoded as "
            "negative."
        ),

        "scope_warning": (
            "These are model-facing benchmark feature "
            "views, not evidence that the 37 loci are "
            "sufficient for supervised model training "
            "or inferential performance estimation. "
            "The benchmark contains only four negative "
            "loci."
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

    # --------------------------------------------------------
    # Output checksum manifest.
    # --------------------------------------------------------

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
            Path(args.core),
            Path(args.extended),
            Path(args.spec),
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
        "MULTIOMIC_FEATURE_VIEWS_BUILD=PASS"
    )

    print(
        "master_rows=37"
    )

    print(
        "master_columns=588"
    )

    print(
        "starting_candidate_predictors=317"
    )

    print(
        "deterministic_redundancy_removals=150"
    )

    print(
        "extended_predictors=167"
    )

    print(
        "extended_table_shape=37x169"
    )

    print(
        "extended_RNA=10"
    )

    print(
        "extended_ATAC=12"
    )

    print(
        "extended_CHIP=145"
    )

    print(
        "core_predictors=67"
    )

    print(
        "core_table_shape=37x69"
    )

    print(
        "core_RNA=5"
    )

    print(
        "core_ATAC=8"
    )

    print(
        "core_CHIP=54"
    )

    print(
        "core_subset_of_extended=YES"
    )

    print(
        "benchmark_relative_predictors=0"
    )

    print(
        "sensitivity_predictors=0"
    )

    print(
        "coordinate_predictors=0"
    )

    print(
        "constant_predictors=0"
    )

    print(
        "predictor_columns_with_missingness=0"
    )

    print(
        "label_used_for_feature_selection=NO"
    )

    print(
        "outcome_driven_selection=NO"
    )

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
