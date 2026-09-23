#!/usr/bin/env python3

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path


WINDOWS = (1000, 10000, 50000)

EXPECTED_INPUT_ROWS = 666
EXPECTED_INPUT_COLUMNS = 46
EXPECTED_LOCI = 37
EXPECTED_ANALYSES = 18

SPMR_COLUMNS = [
    "chip_spmr_locus_mean",
    "chip_spmr_locus_covered_bp",
    "chip_spmr_locus_covered_fraction",
    "chip_spmr_locus_actual_bp",

    "chip_spmr_pm1000bp_mean",
    "chip_spmr_pm1000bp_covered_bp",
    "chip_spmr_pm1000bp_covered_fraction",
    "chip_spmr_pm1000bp_actual_window_bp",

    "chip_spmr_pm10000bp_mean",
    "chip_spmr_pm10000bp_covered_bp",
    "chip_spmr_pm10000bp_covered_fraction",
    "chip_spmr_pm10000bp_actual_window_bp",

    "chip_spmr_pm50000bp_mean",
    "chip_spmr_pm50000bp_covered_bp",
    "chip_spmr_pm50000bp_covered_fraction",
    "chip_spmr_pm50000bp_actual_window_bp",
]


def fail(message):
    raise RuntimeError(message)


def sha256_file(path, chunk_size=1024 * 1024 * 8):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def atomic_text_writer(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    return tmp


def fmt_float(value):
    if not math.isfinite(value):
        fail(f"Non-finite value generated: {value}")
    if value < 0:
        fail(f"Negative value generated: {value}")
    return format(value, ".12g")


def load_fai(path):
    sizes = {}

    with Path(path).open("r", encoding="utf-8") as fh:
        for line_number, line in enumerate(fh, 1):
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 2:
                fail(
                    f"Malformed FAI line {line_number}: "
                    f"expected >=2 fields"
                )

            chrom = fields[0]

            try:
                size = int(fields[1])
            except ValueError:
                fail(
                    f"Malformed FAI size at line {line_number}: "
                    f"{fields[1]}"
                )

            if size <= 0:
                fail(
                    f"Invalid FAI chromosome size at line "
                    f"{line_number}: {size}"
                )

            if chrom in sizes:
                fail(f"Duplicate FAI contig: {chrom}")

            sizes[chrom] = size

    if not sizes:
        fail("FAI is empty")

    return sizes


def load_peak_table(path):
    with Path(path).open(
        "r",
        encoding="utf-8",
        newline="",
    ) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)

    if len(fields) != EXPECTED_INPUT_COLUMNS:
        fail(
            f"Expected {EXPECTED_INPUT_COLUMNS} input columns; "
            f"found {len(fields)}"
        )

    if len(rows) != EXPECTED_INPUT_ROWS:
        fail(
            f"Expected {EXPECTED_INPUT_ROWS} rows; "
            f"found {len(rows)}"
        )

    required = {
        "canonical_locus_id",
        "target_seqname",
        "target_start_1based",
        "target_end_1based",
        "chip_anchor_1based",
        "analysis_id",
    }

    missing = sorted(required - set(fields))

    if missing:
        fail(
            "Missing required peak-table columns: "
            + ", ".join(missing)
        )

    loci = {
        row["canonical_locus_id"]
        for row in rows
    }

    analyses = {
        row["analysis_id"]
        for row in rows
    }

    pairs = {
        (
            row["canonical_locus_id"],
            row["analysis_id"],
        )
        for row in rows
    }

    if len(loci) != EXPECTED_LOCI:
        fail(
            f"Expected {EXPECTED_LOCI} loci; "
            f"found {len(loci)}"
        )

    if len(analyses) != EXPECTED_ANALYSES:
        fail(
            f"Expected {EXPECTED_ANALYSES} analyses; "
            f"found {len(analyses)}"
        )

    if len(pairs) != EXPECTED_INPUT_ROWS:
        fail(
            "Locus-analysis pairs are not unique: "
            f"{len(pairs)} unique pairs"
        )

    return fields, rows


def load_inventory(path):
    with Path(path).open(
        "r",
        encoding="utf-8",
        newline="",
    ) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)

    required = {
        "analysis_id",
        "artifact_type",
        "relative_path",
        "size_bytes",
        "sha256",
        "record_count",
    }

    missing = sorted(required - set(fields))

    if missing:
        fail(
            "Missing inventory columns: "
            + ", ".join(missing)
        )

    spmr = [
        row
        for row in rows
        if row["artifact_type"] == "treat_pileup_spmr"
    ]

    if len(spmr) != EXPECTED_ANALYSES:
        fail(
            f"Expected {EXPECTED_ANALYSES} SPMR artifacts; "
            f"found {len(spmr)}"
        )

    by_analysis = {}

    for row in spmr:
        analysis_id = row["analysis_id"]

        if analysis_id in by_analysis:
            fail(
                f"Duplicate SPMR artifact for {analysis_id}"
            )

        by_analysis[analysis_id] = row

    return by_analysis


def build_regions(rows, chrom_sizes):
    """
    Build the canonical locus plus anchor-centred windows for
    every locus-analysis row.

    Coordinates used internally are 0-based, half-open.
    """

    by_analysis = defaultdict(
        lambda: defaultdict(list)
    )

    region_meta = {}

    for row_index, row in enumerate(rows):
        analysis_id = row["analysis_id"]
        chrom = row["target_seqname"]

        if chrom not in chrom_sizes:
            fail(
                f"{row['canonical_locus_id']}: "
                f"contig {chrom} absent from FAI"
            )

        chrom_size = chrom_sizes[chrom]

        try:
            start1 = int(row["target_start_1based"])
            end1 = int(row["target_end_1based"])
            anchor1 = int(row["chip_anchor_1based"])
        except ValueError as exc:
            fail(
                f"Invalid coordinate for "
                f"{row['canonical_locus_id']}: {exc}"
            )

        if start1 < 1:
            fail(
                f"{row['canonical_locus_id']}: "
                f"start < 1"
            )

        if end1 < start1:
            fail(
                f"{row['canonical_locus_id']}: "
                f"end < start"
            )

        if end1 > chrom_size:
            fail(
                f"{row['canonical_locus_id']}: "
                f"end exceeds chromosome size"
            )

        expected_anchor1 = (start1 + end1) // 2

        if anchor1 != expected_anchor1:
            fail(
                f"{row['canonical_locus_id']}: "
                f"anchor mismatch; observed={anchor1}, "
                f"expected={expected_anchor1}"
            )

        start0 = start1 - 1
        end0 = end1
        anchor0 = anchor1 - 1

        locus_bp = end0 - start0

        if locus_bp != end1 - start1 + 1:
            fail(
                f"{row['canonical_locus_id']}: "
                "canonical locus width mismatch"
            )

        regions = [
            (
                "locus",
                start0,
                end0,
                locus_bp,
            )
        ]

        for half_window in WINDOWS:
            qstart = max(
                0,
                anchor0 - half_window,
            )

            qend = min(
                chrom_size,
                anchor0 + half_window + 1,
            )

            actual_bp = qend - qstart

            if actual_bp <= 0:
                fail(
                    f"{row['canonical_locus_id']}: "
                    f"invalid window +/-{half_window}"
                )

            regions.append(
                (
                    f"pm{half_window}bp",
                    qstart,
                    qend,
                    actual_bp,
                )
            )

        for label, start, end, actual_bp in regions:
            key = (row_index, label)

            region_meta[key] = {
                "weighted_sum": 0.0,
                "covered_bp": 0,
                "actual_bp": actual_bp,
            }

            by_analysis[analysis_id][chrom].append(
                {
                    "row_index": row_index,
                    "label": label,
                    "start": start,
                    "end": end,
                }
            )

    for analysis_id in by_analysis:
        for chrom in by_analysis[analysis_id]:
            by_analysis[analysis_id][chrom].sort(
                key=lambda q: (
                    q["start"],
                    q["end"],
                    q["row_index"],
                    q["label"],
                )
            )

    return by_analysis, region_meta


def verify_artifact(
    analysis_id,
    inventory_row,
    chip_root,
):
    path = (
        Path(chip_root)
        / inventory_row["relative_path"]
    )

    if not path.is_file():
        fail(
            f"{analysis_id}: missing SPMR file: {path}"
        )

    observed_size = path.stat().st_size

    try:
        expected_size = int(
            inventory_row["size_bytes"]
        )
    except ValueError:
        fail(
            f"{analysis_id}: invalid inventory size_bytes"
        )

    if observed_size != expected_size:
        fail(
            f"{analysis_id}: file-size mismatch; "
            f"observed={observed_size}, "
            f"expected={expected_size}"
        )

    observed_sha = sha256_file(path)
    expected_sha = inventory_row["sha256"]

    if observed_sha != expected_sha:
        fail(
            f"{analysis_id}: SHA256 mismatch; "
            f"observed={observed_sha}, "
            f"expected={expected_sha}"
        )

    return path, observed_sha


def stream_spmr_track(
    analysis_id,
    path,
    queries_by_chrom,
    region_meta,
    chrom_sizes,
    expected_record_count=None,
):
    """
    Stream one gzip-compressed bedGraph once.

    BedGraph intervals are required to be:
      - 0-based half-open
      - sorted within chromosome
      - non-overlapping within chromosome
      - finite and non-negative

    Only regions required for the 37 benchmark loci are
    accumulated.
    """

    record_count = 0

    current_chrom = None
    previous_start = None
    previous_end = None

    closed_chroms = set()

    query_state = {}

    for chrom, queries in queries_by_chrom.items():
        query_state[chrom] = {
            "queries": queries,
            "next_index": 0,
            "active": [],
        }

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
    ) as fh:

        for line_number, line in enumerate(fh, 1):
            stripped = line.strip()

            if not stripped:
                continue

            if (
                stripped.startswith("#")
                or stripped.startswith("track")
                or stripped.startswith("browser")
            ):
                continue

            fields = stripped.split()

            if len(fields) < 4:
                fail(
                    f"{analysis_id}: malformed bedGraph "
                    f"line {line_number}: <4 fields"
                )

            chrom = fields[0]

            try:
                start = int(fields[1])
                end = int(fields[2])
                value = float(fields[3])
            except ValueError:
                fail(
                    f"{analysis_id}: non-numeric bedGraph "
                    f"value at line {line_number}"
                )

            if chrom not in chrom_sizes:
                fail(
                    f"{analysis_id}: bedGraph contig "
                    f"{chrom} absent from FAI "
                    f"at line {line_number}"
                )

            if start < 0:
                fail(
                    f"{analysis_id}: negative bedGraph "
                    f"start at line {line_number}"
                )

            if end <= start:
                fail(
                    f"{analysis_id}: end <= start "
                    f"at line {line_number}"
                )

            if end > chrom_sizes[chrom]:
                fail(
                    f"{analysis_id}: interval exceeds "
                    f"{chrom} length at line {line_number}"
                )

            if not math.isfinite(value):
                fail(
                    f"{analysis_id}: non-finite SPMR "
                    f"at line {line_number}"
                )

            if value < 0:
                fail(
                    f"{analysis_id}: negative SPMR "
                    f"at line {line_number}"
                )

            if chrom != current_chrom:
                if current_chrom is not None:
                    closed_chroms.add(current_chrom)

                if chrom in closed_chroms:
                    fail(
                        f"{analysis_id}: chromosome "
                        f"{chrom} reappeared after its block "
                        f"ended at line {line_number}"
                    )

                current_chrom = chrom
                previous_start = None
                previous_end = None

            if previous_start is not None:
                if start < previous_start:
                    fail(
                        f"{analysis_id}: coordinate order "
                        f"decreased at line {line_number}"
                    )

                if start < previous_end:
                    fail(
                        f"{analysis_id}: overlapping "
                        f"bedGraph segments at line "
                        f"{line_number}"
                    )

            previous_start = start
            previous_end = end

            record_count += 1

            state = query_state.get(chrom)

            if state is None:
                continue

            queries = state["queries"]
            next_index = state["next_index"]
            active = state["active"]

            while (
                next_index < len(queries)
                and queries[next_index]["start"] < end
            ):
                active.append(
                    queries[next_index]
                )
                next_index += 1

            if active:
                active = [
                    q
                    for q in active
                    if q["end"] > start
                ]

                for query in active:
                    overlap_start = max(
                        start,
                        query["start"],
                    )

                    overlap_end = min(
                        end,
                        query["end"],
                    )

                    if overlap_end <= overlap_start:
                        continue

                    overlap_bp = (
                        overlap_end - overlap_start
                    )

                    key = (
                        query["row_index"],
                        query["label"],
                    )

                    region_meta[key][
                        "weighted_sum"
                    ] += value * overlap_bp

                    region_meta[key][
                        "covered_bp"
                    ] += overlap_bp

            state["next_index"] = next_index
            state["active"] = active

    if expected_record_count not in (
        None,
        "",
        "NA",
        "na",
    ):
        try:
            expected = int(expected_record_count)
        except ValueError:
            fail(
                f"{analysis_id}: invalid inventory "
                f"record_count={expected_record_count}"
            )

        if record_count != expected:
            fail(
                f"{analysis_id}: bedGraph record-count "
                f"mismatch; observed={record_count}, "
                f"expected={expected}"
            )

    return record_count


def feature_values(row_index, region_meta):
    values = {}

    mapping = [
        (
            "locus",
            "chip_spmr_locus_mean",
            "chip_spmr_locus_covered_bp",
            "chip_spmr_locus_covered_fraction",
            "chip_spmr_locus_actual_bp",
        ),
        (
            "pm1000bp",
            "chip_spmr_pm1000bp_mean",
            "chip_spmr_pm1000bp_covered_bp",
            "chip_spmr_pm1000bp_covered_fraction",
            "chip_spmr_pm1000bp_actual_window_bp",
        ),
        (
            "pm10000bp",
            "chip_spmr_pm10000bp_mean",
            "chip_spmr_pm10000bp_covered_bp",
            "chip_spmr_pm10000bp_covered_fraction",
            "chip_spmr_pm10000bp_actual_window_bp",
        ),
        (
            "pm50000bp",
            "chip_spmr_pm50000bp_mean",
            "chip_spmr_pm50000bp_covered_bp",
            "chip_spmr_pm50000bp_covered_fraction",
            "chip_spmr_pm50000bp_actual_window_bp",
        ),
    ]

    for (
        label,
        mean_col,
        covered_col,
        fraction_col,
        actual_col,
    ) in mapping:

        stats = region_meta[
            (row_index, label)
        ]

        actual_bp = stats["actual_bp"]
        covered_bp = stats["covered_bp"]

        if actual_bp <= 0:
            fail(
                f"row={row_index}, {label}: "
                f"actual_bp <= 0"
            )

        if covered_bp < 0:
            fail(
                f"row={row_index}, {label}: "
                f"covered_bp < 0"
            )

        if covered_bp > actual_bp:
            fail(
                f"row={row_index}, {label}: "
                f"covered_bp={covered_bp} > "
                f"actual_bp={actual_bp}"
            )

        mean_value = (
            stats["weighted_sum"] / actual_bp
        )

        covered_fraction = (
            covered_bp / actual_bp
        )

        if not (
            0.0 <= covered_fraction <= 1.0
        ):
            fail(
                f"row={row_index}, {label}: "
                f"invalid covered_fraction="
                f"{covered_fraction}"
            )

        values[mean_col] = fmt_float(
            mean_value
        )

        values[covered_col] = str(
            covered_bp
        )

        values[fraction_col] = fmt_float(
            covered_fraction
        )

        values[actual_col] = str(
            actual_bp
        )

    return values


def write_output(
    path,
    original_fields,
    original_rows,
    region_meta,
):
    output_fields = (
        list(original_fields)
        + SPMR_COLUMNS
    )

    tmp = atomic_text_writer(path)

    with tmp.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=output_fields,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        for row_index, original in enumerate(
            original_rows
        ):
            out = dict(original)
            out.update(
                feature_values(
                    row_index,
                    region_meta,
                )
            )
            writer.writerow(out)

    os.replace(tmp, path)

    return output_fields


def validate_written_output(
    path,
    original_fields,
    original_rows,
):
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

    expected_fields = (
        list(original_fields)
        + SPMR_COLUMNS
    )

    if fields != expected_fields:
        fail(
            "Output column order differs from expected"
        )

    if len(rows) != EXPECTED_INPUT_ROWS:
        fail(
            f"Output rows={len(rows)}, expected "
            f"{EXPECTED_INPUT_ROWS}"
        )

    inherited_mismatches = 0

    for old, new in zip(
        original_rows,
        rows,
    ):
        for field in original_fields:
            if old[field] != new[field]:
                inherited_mismatches += 1

    loci = {
        r["canonical_locus_id"]
        for r in rows
    }

    analyses = {
        r["analysis_id"]
        for r in rows
    }

    pairs = {
        (
            r["canonical_locus_id"],
            r["analysis_id"],
        )
        for r in rows
    }

    if len(loci) != EXPECTED_LOCI:
        fail(
            f"Output unique loci={len(loci)}"
        )

    if len(analyses) != EXPECTED_ANALYSES:
        fail(
            f"Output unique analyses={len(analyses)}"
        )

    if len(pairs) != EXPECTED_INPUT_ROWS:
        fail(
            f"Output unique pairs={len(pairs)}"
        )

    if inherited_mismatches != 0:
        fail(
            "Inherited peak-feature mismatch count="
            f"{inherited_mismatches}"
        )

    for row in rows:
        start1 = int(
            row["target_start_1based"]
        )
        end1 = int(
            row["target_end_1based"]
        )

        expected_locus_bp = (
            end1 - start1 + 1
        )

        if int(
            row["chip_spmr_locus_actual_bp"]
        ) != expected_locus_bp:
            fail(
                f"{row['canonical_locus_id']}: "
                "SPMR locus width mismatch"
            )

        for half_window in WINDOWS:
            col = (
                f"chip_spmr_pm{half_window}bp_"
                f"actual_window_bp"
            )

            value = int(row[col])

            if value <= 0:
                fail(
                    f"{row['canonical_locus_id']}: "
                    f"invalid {col}={value}"
                )

    return {
        "rows": len(rows),
        "columns": len(fields),
        "unique_loci": len(loci),
        "unique_analyses": len(analyses),
        "unique_locus_analysis_pairs": len(
            pairs
        ),
        "inherited_peak_feature_mismatches":
            inherited_mismatches,
    }


def write_qc(path, metrics):
    tmp = atomic_text_writer(path)

    with tmp.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as fh:
        writer = csv.writer(
            fh,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            ["metric", "value"]
        )

        for key, value in metrics.items():
            writer.writerow(
                [key, value]
            )

    os.replace(tmp, path)


def write_json(path, obj):
    tmp = atomic_text_writer(path)

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


def write_sha256s(path, files):
    tmp = atomic_text_writer(path)

    with tmp.open(
        "w",
        encoding="utf-8",
    ) as fh:
        for file_path in files:
            file_path = Path(file_path)
            fh.write(
                f"{sha256_file(file_path)}  "
                f"{file_path.name}\n"
            )

    os.replace(tmp, path)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Append condition-specific MACS3 SPMR "
            "features to the validated ChIP peak "
            "locus-analysis table."
        )
    )

    parser.add_argument(
        "--peak-table",
        required=True,
    )

    parser.add_argument(
        "--inventory",
        required=True,
    )

    parser.add_argument(
        "--chip-root",
        required=True,
    )

    parser.add_argument(
        "--fai",
        required=True,
    )

    parser.add_argument(
        "--output-table",
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
        Path(args.output_table),
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
            "Refusing to overwrite existing outputs: "
            + ", ".join(existing)
        )

    peak_fields, peak_rows = load_peak_table(
        args.peak_table
    )

    inventory = load_inventory(
        args.inventory
    )

    chrom_sizes = load_fai(
        args.fai
    )

    analysis_ids = {
        row["analysis_id"]
        for row in peak_rows
    }

    if analysis_ids != set(inventory):
        fail(
            "Analysis IDs in peak table and frozen "
            "SPMR inventory do not match exactly"
        )

    (
        queries_by_analysis,
        region_meta,
    ) = build_regions(
        peak_rows,
        chrom_sizes,
    )

    artifact_audit = {}
    track_record_counts = {}

    for number, analysis_id in enumerate(
        sorted(analysis_ids),
        1,
    ):
        inventory_row = inventory[
            analysis_id
        ]

        print(
            f"[{number}/{EXPECTED_ANALYSES}] "
            f"{analysis_id}: verifying frozen "
            f"SPMR artifact",
            flush=True,
        )

        path, observed_sha = verify_artifact(
            analysis_id,
            inventory_row,
            args.chip_root,
        )

        print(
            f"[{number}/{EXPECTED_ANALYSES}] "
            f"{analysis_id}: streaming "
            f"{path.name}",
            flush=True,
        )

        record_count = stream_spmr_track(
            analysis_id,
            path,
            queries_by_analysis[
                analysis_id
            ],
            region_meta,
            chrom_sizes,
            inventory_row.get(
                "record_count",
                "",
            ),
        )

        artifact_audit[analysis_id] = {
            "relative_path":
                inventory_row[
                    "relative_path"
                ],
            "size_bytes":
                int(
                    inventory_row[
                        "size_bytes"
                    ]
                ),
            "sha256": observed_sha,
            "sha256_matches_frozen_inventory":
                True,
        }

        track_record_counts[
            analysis_id
        ] = record_count

    output_fields = write_output(
        args.output_table,
        peak_fields,
        peak_rows,
        region_meta,
    )

    validation = validate_written_output(
        args.output_table,
        peak_fields,
        peak_rows,
    )

    nominal_width_checks = {}

    with Path(args.output_table).open(
        "r",
        encoding="utf-8",
        newline="",
    ) as fh:
        output_rows = list(
            csv.DictReader(
                fh,
                delimiter="\t",
            )
        )

    for half_window in WINDOWS:
        col = (
            f"chip_spmr_pm{half_window}bp_"
            f"actual_window_bp"
        )

        values = sorted({
            int(row[col])
            for row in output_rows
        })

        nominal_width_checks[
            col
        ] = values

    qc_metrics = {
        **validation,
        "input_columns":
            EXPECTED_INPUT_COLUMNS,
        "appended_spmr_columns":
            len(SPMR_COLUMNS),
        "source_spmr_tracks":
            len(artifact_audit),
        "source_sha256_mismatches":
            0,
        "reference_contigs":
            len(chrom_sizes),
        "pm1000bp_actual_window_bp_values":
            ",".join(
                map(
                    str,
                    nominal_width_checks[
                        "chip_spmr_pm1000bp_actual_window_bp"
                    ],
                )
            ),
        "pm10000bp_actual_window_bp_values":
            ",".join(
                map(
                    str,
                    nominal_width_checks[
                        "chip_spmr_pm10000bp_actual_window_bp"
                    ],
                )
            ),
        "pm50000bp_actual_window_bp_values":
            ",".join(
                map(
                    str,
                    nominal_width_checks[
                        "chip_spmr_pm50000bp_actual_window_bp"
                    ],
                )
            ),
        "status": "PASS",
    }

    write_qc(
        args.qc,
        qc_metrics,
    )

    provenance = {
        "workflow": (
            "condition-specific ChIP-seq "
            "SPMR locus integration"
        ),
        "input_peak_table": str(
            Path(args.peak_table)
        ),
        "input_peak_table_sha256":
            sha256_file(
                Path(args.peak_table)
            ),
        "input_inventory": str(
            Path(args.inventory)
        ),
        "input_inventory_sha256":
            sha256_file(
                Path(args.inventory)
            ),
        "reference_fai": str(
            Path(args.fai)
        ),
        "reference_fai_sha256":
            sha256_file(
                Path(args.fai)
            ),
        "coordinate_semantics": {
            "benchmark": "1-based closed",
            "bedgraph": "0-based half-open",
            "canonical_locus": (
                "start0 = target_start_1based - 1; "
                "end0 = target_end_1based"
            ),
            "anchor": (
                "anchor1 = floor((start1 + end1)/2); "
                "anchor0 = anchor1 - 1"
            ),
            "windows": (
                "anchor +/- N bp including anchor; "
                "nominal width = 2*N + 1"
            ),
        },
        "spmr_semantics": {
            "mean": (
                "sum(signal_value * overlap_bp) "
                "/ actual_region_bp"
            ),
            "uncovered_positions": (
                "contribute zero to the full-region mean"
            ),
            "covered_bp": (
                "bp represented by valid non-overlapping "
                "bedGraph segments within region"
            ),
        },
        "windows_bp": list(WINDOWS),
        "output_columns": output_fields,
        "source_artifacts":
            artifact_audit,
        "track_record_counts":
            track_record_counts,
        "validation":
            qc_metrics,
    }

    write_json(
        args.provenance,
        provenance,
    )

    write_sha256s(
        args.sha256s,
        [
            args.output_table,
            args.qc,
            args.provenance,
        ],
    )

    print()
    print("CHIP_SPMR_INTEGRATION=PASS")
    print(
        f"rows={validation['rows']}"
    )
    print(
        f"columns={validation['columns']}"
    )
    print(
        "unique_locus_analysis_pairs="
        f"{validation['unique_locus_analysis_pairs']}"
    )
    print(
        "inherited_peak_feature_mismatches="
        f"{validation['inherited_peak_feature_mismatches']}"
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
