#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import os
import sys
from collections import Counter
from pathlib import Path


EXPECTED_SHAPES = {
    "RNA": (37, 43),
    "ATAC": (37, 47),
    "CHIP": (37, 524),
}

EXPECTED_NONCORE = {
    "RNA": 30,
    "ATAC": 34,
    "CHIP": 511,
}

EXPECTED_OUTPUT_ROWS = 37
EXPECTED_OUTPUT_COLUMNS = 588
EXPECTED_ROUNDTRIP_COMPARISONS = 22718
EXPECTED_CHIP_RELATIVE_PERCENTILES = 96


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


EXPECTED_ROLES = Counter({
    "positive": 27,
    "negative": 4,
    "support_only": 6,
})


EXPECTED_SEMANTICS = Counter({
    "point_coordinate": 22,
    "interbase_cut_boundary": 4,
    "paired_nick_interval": 3,
    "mapped_interval": 8,
})


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


def validate_table(
    name,
    fields,
    rows,
):

    expected_rows, expected_cols = (
        EXPECTED_SHAPES[name]
    )

    if len(rows) != expected_rows:
        fail(
            f"{name}: expected {expected_rows} rows; "
            f"found {len(rows)}"
        )

    if len(fields) != expected_cols:
        fail(
            f"{name}: expected {expected_cols} columns; "
            f"found {len(fields)}"
        )

    missing = [
        field
        for field in CORE_FIELDS
        if field not in fields
    ]

    if missing:
        fail(
            f"{name}: missing core fields: "
            + ",".join(missing)
        )

    loci = [
        row["canonical_locus_id"]
        for row in rows
    ]

    if len(set(loci)) != 37:
        fail(
            f"{name}: expected 37 unique loci; "
            f"found {len(set(loci))}"
        )

    noncore = [
        field
        for field in fields
        if field not in CORE_FIELDS
    ]

    if len(noncore) != EXPECTED_NONCORE[name]:
        fail(
            f"{name}: expected "
            f"{EXPECTED_NONCORE[name]} non-core fields; "
            f"found {len(noncore)}"
        )

    roles = Counter(
        row["benchmark_role"]
        for row in rows
    )

    semantics = Counter(
        row["locus_semantic_class"]
        for row in rows
    )

    assemblies = Counter(
        row["target_assembly"]
        for row in rows
    )

    if roles != EXPECTED_ROLES:
        fail(
            f"{name}: unexpected role distribution"
        )

    if semantics != EXPECTED_SEMANTICS:
        fail(
            f"{name}: unexpected semantic distribution"
        )

    if assemblies != Counter({
        "GCF_003668045.3": 37,
    }):
        fail(
            f"{name}: unexpected assembly distribution"
        )

    return noncore


def build_master(
    tables,
    noncore_fields,
):

    locus_maps = {}

    for name in (
        "RNA",
        "ATAC",
        "CHIP",
    ):

        locus_maps[name] = {
            row["canonical_locus_id"]: row
            for row in tables[name]["rows"]
        }

    locus_sets = {
        name: set(mapping)
        for name, mapping in locus_maps.items()
    }

    if not (
        locus_sets["RNA"]
        == locus_sets["ATAC"]
        == locus_sets["CHIP"]
    ):
        fail(
            "The three modalities do not contain "
            "the same locus set"
        )

    # Explicit namespace collision check.
    for a, b in (
        ("RNA", "ATAC"),
        ("RNA", "CHIP"),
        ("ATAC", "CHIP"),
    ):

        overlap = (
            set(noncore_fields[a])
            & set(noncore_fields[b])
        )

        if overlap:
            fail(
                f"Non-core namespace collision "
                f"{a}/{b}: "
                + ",".join(
                    sorted(overlap)
                )
            )

    output_fields = list(
        CORE_FIELDS
    )

    for name in (
        "RNA",
        "ATAC",
        "CHIP",
    ):
        output_fields.extend(
            noncore_fields[name]
        )

    if len(output_fields) != EXPECTED_OUTPUT_COLUMNS:
        fail(
            "Expected "
            f"{EXPECTED_OUTPUT_COLUMNS} master columns; "
            f"found {len(output_fields)}"
        )

    if len(set(output_fields)) != len(output_fields):
        fail(
            "Duplicate columns in master schema"
        )

    output_rows = []

    identity_comparisons = 0
    identity_errors = 0

    # Preserve RNA benchmark ordering explicitly.
    locus_order = [
        row["canonical_locus_id"]
        for row in tables["RNA"]["rows"]
    ]

    for locus in locus_order:

        rna = locus_maps["RNA"][locus]
        atac = locus_maps["ATAC"][locus]
        chip = locus_maps["CHIP"][locus]

        out = {}

        for field in CORE_FIELDS:

            reference = rna[field]

            for other in (
                atac,
                chip,
            ):

                identity_comparisons += 1

                if other[field] != reference:
                    identity_errors += 1

            out[field] = reference

        for field in noncore_fields["RNA"]:
            out[field] = rna[field]

        for field in noncore_fields["ATAC"]:
            out[field] = atac[field]

        for field in noncore_fields["CHIP"]:
            out[field] = chip[field]

        output_rows.append(out)

    expected_identity_comparisons = (
        37
        * len(CORE_FIELDS)
        * 2
    )

    if (
        identity_comparisons
        != expected_identity_comparisons
    ):
        fail(
            "Unexpected identity comparison count: "
            f"{identity_comparisons}"
        )

    if identity_errors:
        fail(
            f"Cross-modal identity errors="
            f"{identity_errors}"
        )

    return (
        output_fields,
        output_rows,
        locus_maps,
        identity_comparisons,
    )


def roundtrip_validate(
    tables,
    noncore_fields,
    master_rows,
):

    master_by_locus = {
        row["canonical_locus_id"]: row
        for row in master_rows
    }

    comparisons_by_modality = {}
    errors_by_modality = {}

    total_comparisons = 0
    total_errors = 0

    for name in (
        "RNA",
        "ATAC",
        "CHIP",
    ):

        fields = tables[name]["fields"]
        source_rows = tables[name]["rows"]

        comparisons = 0
        errors = 0

        for source in source_rows:

            locus = source[
                "canonical_locus_id"
            ]

            if locus not in master_by_locus:
                errors += 1
                continue

            master = master_by_locus[locus]

            reconstructed = {}

            for field in CORE_FIELDS:
                reconstructed[field] = (
                    master[field]
                )

            for field in noncore_fields[name]:
                reconstructed[field] = (
                    master[field]
                )

            if set(reconstructed) != set(fields):
                errors += 1
                continue

            for field in fields:

                comparisons += 1

                if (
                    reconstructed[field]
                    != source[field]
                ):
                    errors += 1

        comparisons_by_modality[name] = (
            comparisons
        )

        errors_by_modality[name] = errors

        total_comparisons += comparisons
        total_errors += errors

    if (
        total_comparisons
        != EXPECTED_ROUNDTRIP_COMPARISONS
    ):
        fail(
            "Expected "
            f"{EXPECTED_ROUNDTRIP_COMPARISONS} "
            "roundtrip comparisons; found "
            f"{total_comparisons}"
        )

    if total_errors:
        fail(
            f"Roundtrip errors={total_errors}"
        )

    return {
        "RNA_roundtrip_comparisons":
            comparisons_by_modality["RNA"],
        "RNA_roundtrip_errors":
            errors_by_modality["RNA"],

        "ATAC_roundtrip_comparisons":
            comparisons_by_modality["ATAC"],
        "ATAC_roundtrip_errors":
            errors_by_modality["ATAC"],

        "CHIP_roundtrip_comparisons":
            comparisons_by_modality["CHIP"],
        "CHIP_roundtrip_errors":
            errors_by_modality["CHIP"],

        "total_roundtrip_comparisons":
            total_comparisons,
        "total_roundtrip_errors":
            total_errors,
    }


def validate_master(
    fields,
    rows,
):

    if len(rows) != EXPECTED_OUTPUT_ROWS:
        fail(
            f"Expected {EXPECTED_OUTPUT_ROWS} master rows; "
            f"found {len(rows)}"
        )

    if len(fields) != EXPECTED_OUTPUT_COLUMNS:
        fail(
            f"Expected {EXPECTED_OUTPUT_COLUMNS} master "
            f"columns; found {len(fields)}"
        )

    loci = [
        row["canonical_locus_id"]
        for row in rows
    ]

    if len(set(loci)) != 37:
        fail(
            "Master locus uniqueness failed"
        )

    roles = Counter(
        row["benchmark_role"]
        for row in rows
    )

    semantics = Counter(
        row["locus_semantic_class"]
        for row in rows
    )

    if roles != EXPECTED_ROLES:
        fail(
            "Master role distribution failed"
        )

    if semantics != EXPECTED_SEMANTICS:
        fail(
            "Master semantic distribution failed"
        )

    relative_percentiles = [
        field
        for field in fields
        if "relative_percentile_" in field
    ]

    if (
        len(relative_percentiles)
        != EXPECTED_CHIP_RELATIVE_PERCENTILES
    ):
        fail(
            "Expected "
            f"{EXPECTED_CHIP_RELATIVE_PERCENTILES} "
            "benchmark-relative percentile columns; "
            f"found {len(relative_percentiles)}"
        )

    return {
        "output_rows": len(rows),
        "output_columns": len(fields),
        "output_unique_loci":
            len(set(loci)),
        "benchmark_positive":
            roles["positive"],
        "benchmark_negative":
            roles["negative"],
        "benchmark_support_only":
            roles["support_only"],
        "semantic_point_coordinate":
            semantics["point_coordinate"],
        "semantic_interbase_cut_boundary":
            semantics[
                "interbase_cut_boundary"
            ],
        "semantic_paired_nick_interval":
            semantics[
                "paired_nick_interval"
            ],
        "semantic_mapped_interval":
            semantics["mapped_interval"],
        "chip_benchmark_relative_percentile_columns":
            len(relative_percentiles),
    }


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--rna",
        required=True,
    )

    parser.add_argument(
        "--atac",
        required=True,
    )

    parser.add_argument(
        "--chip",
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
            + ",".join(existing)
        )

    input_paths = {
        "RNA": Path(args.rna),
        "ATAC": Path(args.atac),
        "CHIP": Path(args.chip),
    }

    tables = {}

    noncore_fields = {}

    for name in (
        "RNA",
        "ATAC",
        "CHIP",
    ):

        fields, rows = read_tsv(
            input_paths[name]
        )

        tables[name] = {
            "fields": fields,
            "rows": rows,
        }

        noncore_fields[name] = (
            validate_table(
                name,
                fields,
                rows,
            )
        )

    (
        output_fields,
        output_rows,
        locus_maps,
        identity_comparisons,
    ) = build_master(
        tables,
        noncore_fields,
    )

    output_metrics = validate_master(
        output_fields,
        output_rows,
    )

    roundtrip = roundtrip_validate(
        tables,
        noncore_fields,
        output_rows,
    )

    metrics = {
        "RNA_source_rows":
            len(tables["RNA"]["rows"]),
        "RNA_source_columns":
            len(tables["RNA"]["fields"]),
        "RNA_noncore_columns":
            len(noncore_fields["RNA"]),

        "ATAC_source_rows":
            len(tables["ATAC"]["rows"]),
        "ATAC_source_columns":
            len(tables["ATAC"]["fields"]),
        "ATAC_noncore_columns":
            len(noncore_fields["ATAC"]),

        "CHIP_source_rows":
            len(tables["CHIP"]["rows"]),
        "CHIP_source_columns":
            len(tables["CHIP"]["fields"]),
        "CHIP_noncore_columns":
            len(noncore_fields["CHIP"]),

        "core_identity_columns":
            len(CORE_FIELDS),
        "cross_modal_identity_comparisons":
            identity_comparisons,

        **output_metrics,
        **roundtrip,

        "validation_errors": 0,
        "status": "PASS",
    }

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
        [
            "metric",
            "value",
        ],
        qc_rows,
    )

    provenance = {
        "workflow":
            "37-locus RNA + ATAC + ChIP master integration",

        "producer_script":
            str(Path(sys.argv[0])),

        "producer_script_sha256":
            sha256_file(sys.argv[0]),

        "sources": {
            name: {
                "path":
                    str(input_paths[name]),
                "sha256":
                    sha256_file(
                        input_paths[name]
                    ),
                "rows":
                    len(
                        tables[name]["rows"]
                    ),
                "columns":
                    len(
                        tables[name]["fields"]
                    ),
                "noncore_columns":
                    len(
                        noncore_fields[name]
                    ),
            }
            for name in (
                "RNA",
                "ATAC",
                "CHIP",
            )
        },

        "output_shape": {
            "rows":
                len(output_rows),
            "columns":
                len(output_fields),
        },

        "join_key":
            "canonical_locus_id",

        "output_order_policy": (
            "The output preserves the validated RNA "
            "37-locus row order. ATAC and ChIP are "
            "joined by canonical_locus_id; no "
            "integration step relies on row position."
        ),

        "shared_identity_fields":
            CORE_FIELDS,

        "feature_blocks": {
            "RNA":
                noncore_fields["RNA"],
            "ATAC":
                noncore_fields["ATAC"],
            "CHIP":
                noncore_fields["CHIP"],
        },

        "namespace_policy": (
            "Original validated source column names "
            "are preserved. Integration is permitted "
            "only because the three non-core feature "
            "namespaces are disjoint."
        ),

        "aggregation_policy": (
            "No new biological aggregation, scaling, "
            "imputation or feature transformation is "
            "performed. This table is a deterministic "
            "keyed join of three frozen 37-locus "
            "context tables."
        ),

        "lossless_validation": (
            "Each of the three source tables was "
            "reconstructed from the master table and "
            "every source field was compared as an "
            "exact string."
        ),

        "master_table_scope": (
            "This file is a complete multi-omic "
            "benchmark context table, not a finalized "
            "machine-learning feature matrix."
        ),

        "benchmark_relative_feature_warning": (
            "The 96 ChIP columns containing "
            "'relative_percentile_' are relative to "
            "the distribution of the 37 benchmark "
            "loci used to construct the ChIP "
            "harmonization layer. They are retained "
            "for descriptive traceability in this "
            "master table but must not automatically "
            "be treated as portable model predictors."
        ),

        "label_and_metadata_warning": (
            "benchmark_role and other evidence, "
            "provenance, study or technical metadata "
            "are retained in the master table and "
            "must be separated from predictors when "
            "a model-facing feature view is built."
        ),

        "validation":
            metrics,
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
        "MULTIOMIC_MASTER_BUILD=PASS"
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
