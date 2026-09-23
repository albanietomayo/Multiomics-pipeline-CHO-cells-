#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED_SOURCE_ROWS = 222
EXPECTED_SOURCE_COLUMNS = 100
EXPECTED_OUTPUT_ROWS = 37
EXPECTED_OUTPUT_COLUMNS = 524
EXPECTED_LOCI = 37
EXPECTED_MARKS = 6
EXPECTED_NA_CELLS = 1998


MARK_ORDER = [
    "H3K27ac",
    "H3K27me3",
    "H3K36me3",
    "H3K4me1",
    "H3K4me3",
    "H3K9me3",
]


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


EXPECTED_ROLE_COUNTS = {
    "positive": 27,
    "negative": 4,
    "support_only": 6,
}


EXPECTED_SEMANTIC_COUNTS = {
    "point_coordinate": 22,
    "interbase_cut_boundary": 4,
    "paired_nick_interval": 3,
    "mapped_interval": 8,
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


def prefixed(mark, field):
    return f"{mark}__{field}"


def validate_source(fields, rows):

    if len(rows) != EXPECTED_SOURCE_ROWS:
        fail(
            f"Expected {EXPECTED_SOURCE_ROWS} input rows; "
            f"found {len(rows)}"
        )

    if len(fields) != EXPECTED_SOURCE_COLUMNS:
        fail(
            f"Expected {EXPECTED_SOURCE_COLUMNS} input columns; "
            f"found {len(fields)}"
        )

    missing_identity = [
        field
        for field in IDENTITY_FIELDS
        if field not in fields
    ]

    if missing_identity:
        fail(
            "Missing identity fields: "
            + ", ".join(missing_identity)
        )

    if "target" not in fields:
        fail("Missing target field")

    marks = Counter(
        row["target"]
        for row in rows
    )

    expected_marks = Counter({
        mark: EXPECTED_LOCI
        for mark in MARK_ORDER
    })

    if marks != expected_marks:
        fail(
            "Unexpected mark distribution: "
            + repr(dict(marks))
        )

    keys = [
        (
            row["canonical_locus_id"],
            row["target"],
        )
        for row in rows
    ]

    if len(set(keys)) != EXPECTED_SOURCE_ROWS:
        fail(
            "Duplicate locus x mark key detected"
        )

    loci = {
        row["canonical_locus_id"]
        for row in rows
    }

    if len(loci) != EXPECTED_LOCI:
        fail(
            f"Expected {EXPECTED_LOCI} loci; "
            f"found {len(loci)}"
        )

    na_cells = sum(
        value == "NA"
        for row in rows
        for value in row.values()
    )

    if na_cells != EXPECTED_NA_CELLS:
        fail(
            f"Expected {EXPECTED_NA_CELLS} NA cells; "
            f"found {na_cells}"
        )

    return {
        "unique_loci": len(loci),
        "mark_counts": dict(marks),
        "source_na_cells": na_cells,
    }


def pivot(fields, rows):

    mark_specific_fields = [
        field
        for field in fields
        if (
            field not in IDENTITY_FIELDS
            and field != "target"
        )
    ]

    if len(mark_specific_fields) != 85:
        fail(
            "Expected 85 mark-specific fields; "
            f"found {len(mark_specific_fields)}"
        )

    output_fields = list(
        IDENTITY_FIELDS
    )

    for mark in MARK_ORDER:
        for field in mark_specific_fields:
            output_fields.append(
                prefixed(mark, field)
            )

    if len(output_fields) != EXPECTED_OUTPUT_COLUMNS:
        fail(
            f"Expected {EXPECTED_OUTPUT_COLUMNS} output "
            f"columns; found {len(output_fields)}"
        )

    by_locus = defaultdict(list)

    for row in rows:
        by_locus[
            row["canonical_locus_id"]
        ].append(row)

    if len(by_locus) != EXPECTED_LOCI:
        fail(
            f"Expected {EXPECTED_LOCI} locus groups; "
            f"found {len(by_locus)}"
        )

    output_rows = []

    identity_errors = 0
    missing_mark_errors = 0

    for locus in sorted(by_locus):

        locus_rows = by_locus[locus]

        mark_map = {
            row["target"]: row
            for row in locus_rows
        }

        if set(mark_map) != set(MARK_ORDER):
            missing_mark_errors += 1
            continue

        out = {}

        for field in IDENTITY_FIELDS:

            values = {
                row[field]
                for row in locus_rows
            }

            if len(values) != 1:
                identity_errors += 1
                fail(
                    f"{locus}: identity field "
                    f"{field} differs across marks"
                )

            out[field] = next(
                iter(values)
            )

        for mark in MARK_ORDER:

            source_row = mark_map[mark]

            for field in mark_specific_fields:
                out[
                    prefixed(mark, field)
                ] = source_row[field]

        output_rows.append(out)

    if identity_errors:
        fail(
            f"identity_errors={identity_errors}"
        )

    if missing_mark_errors:
        fail(
            f"missing_mark_errors={missing_mark_errors}"
        )

    return (
        output_fields,
        mark_specific_fields,
        output_rows,
    )


def roundtrip_validate(
    source_fields,
    source_rows,
    mark_specific_fields,
    output_rows,
):
    """
    Reconstruct all 222 source rows from the 37-row wide table.

    The target column is recovered from the column-block prefix.
    Every one of the 100 fields is then compared as an exact
    string to the source table.
    """

    source_by_key = {
        (
            row["canonical_locus_id"],
            row["target"],
        ): row
        for row in source_rows
    }

    comparisons = 0
    errors = 0

    reconstructed_rows = 0

    for wide in output_rows:

        locus = wide[
            "canonical_locus_id"
        ]

        for mark in MARK_ORDER:

            key = (
                locus,
                mark,
            )

            if key not in source_by_key:
                errors += 1
                continue

            reconstructed = {}

            for field in IDENTITY_FIELDS:
                reconstructed[field] = (
                    wide[field]
                )

            reconstructed["target"] = mark

            for field in mark_specific_fields:
                reconstructed[field] = (
                    wide[
                        prefixed(mark, field)
                    ]
                )

            source = source_by_key[key]

            if set(reconstructed) != set(
                source_fields
            ):
                errors += 1
                continue

            reconstructed_rows += 1

            for field in source_fields:

                comparisons += 1

                if (
                    reconstructed[field]
                    != source[field]
                ):
                    errors += 1

    if reconstructed_rows != EXPECTED_SOURCE_ROWS:
        fail(
            "Roundtrip reconstructed "
            f"{reconstructed_rows} rows; expected "
            f"{EXPECTED_SOURCE_ROWS}"
        )

    expected_comparisons = (
        EXPECTED_SOURCE_ROWS
        * EXPECTED_SOURCE_COLUMNS
    )

    if comparisons != expected_comparisons:
        fail(
            f"Roundtrip comparisons={comparisons}; "
            f"expected {expected_comparisons}"
        )

    if errors:
        fail(
            f"Roundtrip validation errors={errors}"
        )

    return {
        "roundtrip_reconstructed_rows":
            reconstructed_rows,
        "roundtrip_field_comparisons":
            comparisons,
        "roundtrip_errors":
            errors,
    }


def validate_output(
    fields,
    rows,
):

    if len(rows) != EXPECTED_OUTPUT_ROWS:
        fail(
            f"Expected {EXPECTED_OUTPUT_ROWS} output rows; "
            f"found {len(rows)}"
        )

    if len(fields) != EXPECTED_OUTPUT_COLUMNS:
        fail(
            f"Expected {EXPECTED_OUTPUT_COLUMNS} output "
            f"columns; found {len(fields)}"
        )

    loci = [
        row["canonical_locus_id"]
        for row in rows
    ]

    if len(set(loci)) != EXPECTED_LOCI:
        fail(
            "Duplicate or missing locus in wide output"
        )

    role_counts = Counter(
        row["benchmark_role"]
        for row in rows
    )

    if dict(role_counts) != EXPECTED_ROLE_COUNTS:
        fail(
            "Unexpected benchmark-role distribution: "
            + repr(dict(role_counts))
        )

    semantic_counts = Counter(
        row["locus_semantic_class"]
        for row in rows
    )

    if dict(semantic_counts) != EXPECTED_SEMANTIC_COUNTS:
        fail(
            "Unexpected locus-semantic distribution: "
            + repr(dict(semantic_counts))
        )

    output_na_cells = sum(
        value == "NA"
        for row in rows
        for value in row.values()
    )

    if output_na_cells != EXPECTED_NA_CELLS:
        fail(
            f"Expected {EXPECTED_NA_CELLS} NA cells in "
            f"wide output; found {output_na_cells}"
        )

    return {
        "rows": len(rows),
        "columns": len(fields),
        "unique_loci": len(set(loci)),
        "benchmark_positive": role_counts["positive"],
        "benchmark_negative": role_counts["negative"],
        "benchmark_support_only":
            role_counts["support_only"],
        "semantic_point_coordinate":
            semantic_counts["point_coordinate"],
        "semantic_interbase_cut_boundary":
            semantic_counts["interbase_cut_boundary"],
        "semantic_paired_nick_interval":
            semantic_counts["paired_nick_interval"],
        "semantic_mapped_interval":
            semantic_counts["mapped_interval"],
        "output_na_cells": output_na_cells,
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

    source_fields, source_rows = read_tsv(
        args.input
    )

    source_metrics = validate_source(
        source_fields,
        source_rows,
    )

    (
        output_fields,
        mark_specific_fields,
        output_rows,
    ) = pivot(
        source_fields,
        source_rows,
    )

    output_metrics = validate_output(
        output_fields,
        output_rows,
    )

    roundtrip_metrics = roundtrip_validate(
        source_fields,
        source_rows,
        mark_specific_fields,
        output_rows,
    )

    metrics = {
        "source_rows":
            len(source_rows),
        "source_columns":
            len(source_fields),
        "source_unique_loci":
            source_metrics["unique_loci"],
        "source_mark_specific_fields":
            len(mark_specific_fields),
        "source_na_cells":
            source_metrics["source_na_cells"],
        "output_rows":
            output_metrics["rows"],
        "output_columns":
            output_metrics["columns"],
        "output_unique_loci":
            output_metrics["unique_loci"],
        "output_na_cells":
            output_metrics["output_na_cells"],
        "benchmark_positive":
            output_metrics["benchmark_positive"],
        "benchmark_negative":
            output_metrics["benchmark_negative"],
        "benchmark_support_only":
            output_metrics["benchmark_support_only"],
        "semantic_point_coordinate":
            output_metrics[
                "semantic_point_coordinate"
            ],
        "semantic_interbase_cut_boundary":
            output_metrics[
                "semantic_interbase_cut_boundary"
            ],
        "semantic_paired_nick_interval":
            output_metrics[
                "semantic_paired_nick_interval"
            ],
        "semantic_mapped_interval":
            output_metrics[
                "semantic_mapped_interval"
            ],
        **roundtrip_metrics,
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
            "ChIP-seq benchmark locus wide-context pivot",
        "source_table":
            str(Path(args.input)),
        "source_table_sha256":
            sha256_file(args.input),
        "producer_script":
            str(Path(sys.argv[0])),
        "producer_script_sha256":
            sha256_file(sys.argv[0]),
        "source_shape": {
            "rows": len(source_rows),
            "columns": len(source_fields),
        },
        "output_shape": {
            "rows": len(output_rows),
            "columns": len(output_fields),
        },
        "identity_fields":
            IDENTITY_FIELDS,
        "mark_order":
            MARK_ORDER,
        "mark_specific_fields":
            mark_specific_fields,
        "column_naming": (
            "<histone_mark>__<mark_level_field>"
        ),
        "target_handling": (
            "The source target column is encoded in the "
            "histone-mark column prefix and is therefore "
            "not duplicated as six target columns."
        ),
        "aggregation_policy": (
            "No additional biological aggregation is "
            "performed. The transformation is a "
            "deterministic wide pivot from one row per "
            "locus x histone mark to one row per locus."
        ),
        "cross_mark_policy": (
            "Histone marks remain separate feature "
            "blocks. No averaging, summation, scoring "
            "or other combination across marks is "
            "performed."
        ),
        "missing_value_policy": (
            "Literal NA values from the validated "
            "mark-level source are preserved unchanged."
        ),
        "roundtrip_validation": (
            "The 37-row wide table was expanded back "
            "to the 222 source locus x mark records and "
            "all 100 source fields were compared as "
            "exact strings."
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

    print("CHIP_LOCUS_CONTEXT_PIVOT=PASS")

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
