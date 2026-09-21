#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_benchmark_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    required = {
        "canonical_locus_id",
        "benchmark_role",
        "target_assembly",
        "target_seqname",
        "target_start",
        "target_end",
    }

    missing = required - set(reader.fieldnames or [])
    if missing:
        raise ValueError(
            f"Missing benchmark columns: {sorted(missing)}"
        )

    seen = set()

    for row in rows:
        locus_id = row["canonical_locus_id"]

        if locus_id in seen:
            raise ValueError(f"Duplicate locus: {locus_id}")

        seen.add(locus_id)

        start = int(row["target_start"])
        end = int(row["target_end"])

        if start < 1 or end < start:
            raise ValueError(
                f"Invalid coordinates for {locus_id}: "
                f"{start}-{end}"
            )

    return rows


def read_benchmark_bed(path):
    rows = {}

    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 4:
                raise ValueError(
                    f"BED line {line_number}: expected >=4 columns"
                )

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            locus_id = fields[3]

            if locus_id in rows:
                raise ValueError(
                    f"Duplicate BED locus ID: {locus_id}"
                )

            rows[locus_id] = (chrom, start, end)

    return rows


def validate_benchmark_coordinates(tsv_rows, bed_rows):

    if len(tsv_rows) != len(bed_rows):
        raise ValueError(
            "Benchmark TSV/BED contain different numbers of loci"
        )

    for row in tsv_rows:

        locus_id = row["canonical_locus_id"]

        if locus_id not in bed_rows:
            raise ValueError(
                f"{locus_id} missing from benchmark BED"
            )

        expected = (
            row["target_seqname"],
            int(row["target_start"]) - 1,
            int(row["target_end"]),
        )

        observed = bed_rows[locus_id]

        if expected != observed:
            raise ValueError(
                f"TSV/BED mismatch for {locus_id}: "
                f"expected={expected}, observed={observed}"
            )


def read_idr(path):

    raw_count = 0
    unique = set()

    with open(path, encoding="utf-8") as handle:

        for line_number, line in enumerate(handle, 1):

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 3:
                raise ValueError(
                    f"{path}:{line_number}: fewer than 3 columns"
                )

            chrom = fields[0]

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: invalid BED coordinates"
                ) from exc

            if start < 0 or end <= start:
                raise ValueError(
                    f"{path}:{line_number}: invalid interval "
                    f"{chrom}:{start}-{end}"
                )

            raw_count += 1
            unique.add((chrom, start, end))

    by_chrom = defaultdict(list)

    for chrom, start, end in unique:
        by_chrom[chrom].append((start, end))

    for chrom in by_chrom:
        by_chrom[chrom].sort()

    return raw_count, unique, by_chrom


def merge_intervals(intervals):

    if not intervals:
        return []

    merged = []

    current_start, current_end = intervals[0]

    for start, end in intervals[1:]:

        if start <= current_end:
            current_end = max(current_end, end)

        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end

    merged.append((current_start, current_end))

    return merged


def overlapping(intervals, qstart, qend):
    return [
        (start, end)
        for start, end in intervals
        if start < qend and end > qstart
    ]


def overlap_bp(intervals, qstart, qend):

    total = 0

    for start, end in merge_intervals(
        overlapping(intervals, qstart, qend)
    ):
        total += max(
            0,
            min(end, qend) - max(start, qstart)
        )

    return total


def nearest_distance(intervals, anchor):

    if not intervals:
        return None

    best = None

    for start, end in intervals:

        if start <= anchor < end:
            return 0

        if anchor < start:
            distance = start - anchor
        else:
            distance = anchor - (end - 1)

        if best is None or distance < best:
            best = distance

    return best


def write_unique_bed(path, unique):

    rows = sorted(
        unique,
        key=lambda x: (x[0], x[1], x[2])
    )

    with open(path, "w", encoding="utf-8") as handle:
        for chrom, start, end in rows:
            handle.write(
                f"{chrom}\t{start}\t{end}\n"
            )


def add_peak_features(
    out,
    prefix,
    chrom_intervals,
    locus_start,
    locus_end,
    anchor,
    windows,
):

    intervals = chrom_intervals

    exact = overlapping(
        intervals,
        locus_start,
        locus_end,
    )

    out[f"{prefix}_exact_overlap"] = int(bool(exact))
    out[f"{prefix}_exact_n_unique_peaks"] = len(exact)

    out[f"{prefix}_exact_overlap_bp"] = overlap_bp(
        intervals,
        locus_start,
        locus_end,
    )

    anchor_hits = overlapping(
        intervals,
        anchor,
        anchor + 1,
    )

    out[f"{prefix}_anchor_in_peak"] = int(
        bool(anchor_hits)
    )

    distance = nearest_distance(
        intervals,
        anchor,
    )

    out[f"{prefix}_nearest_peak_distance_bp"] = (
        "" if distance is None else distance
    )

    for half_window in windows:

        qstart = max(0, anchor - half_window)
        qend = anchor + half_window + 1

        hits = overlapping(
            intervals,
            qstart,
            qend,
        )

        label = f"pm{half_window}bp"

        out[f"{prefix}_n_unique_peaks_{label}"] = len(hits)

        out[f"{prefix}_accessible_bp_{label}"] = overlap_bp(
            intervals,
            qstart,
            qend,
        )


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--benchmark-tsv", required=True)
    parser.add_argument("--benchmark-bed", required=True)

    parser.add_argument("--idr10", required=True)
    parser.add_argument("--idr05", required=True)

    parser.add_argument("--windows", required=True)

    parser.add_argument("--features", required=True)
    parser.add_argument("--anchors-bed", required=True)
    parser.add_argument("--idr10-unique-bed", required=True)
    parser.add_argument("--idr05-unique-bed", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--provenance", required=True)

    args = parser.parse_args()

    windows = sorted(
        {
            int(x)
            for x in args.windows.split(",")
            if x.strip()
        }
    )

    if not windows or any(x <= 0 for x in windows):
        raise ValueError(
            "Windows must contain positive integer half-window sizes"
        )

    benchmark = read_benchmark_tsv(
        args.benchmark_tsv
    )

    benchmark_bed = read_benchmark_bed(
        args.benchmark_bed
    )

    validate_benchmark_coordinates(
        benchmark,
        benchmark_bed,
    )

    raw10, unique10, by_chrom10 = read_idr(
        args.idr10
    )

    raw05, unique05, by_chrom05 = read_idr(
        args.idr05
    )

    # IDR <= 0.05 should be a genomic subset of IDR <= 0.10.
    if not unique05.issubset(unique10):
        raise ValueError(
            "IDR <= 0.05 unique intervals are not a subset "
            "of IDR <= 0.10 intervals"
        )

    for output in (
        args.features,
        args.anchors_bed,
        args.idr10_unique_bed,
        args.idr05_unique_bed,
        args.summary,
        args.provenance,
    ):
        Path(output).parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    write_unique_bed(
        args.idr10_unique_bed,
        unique10,
    )

    write_unique_bed(
        args.idr05_unique_bed,
        unique05,
    )

    output_rows = []

    with open(
        args.anchors_bed,
        "w",
        encoding="utf-8",
    ) as anchors:

        for row in benchmark:

            start1 = int(row["target_start"])
            end1 = int(row["target_end"])

            start0 = start1 - 1
            end0 = end1

            # Integer centre of the 1-based closed experimental locus.
            anchor1 = (start1 + end1) // 2
            anchor0 = anchor1 - 1

            locus_id = row["canonical_locus_id"]
            chrom = row["target_seqname"]

            anchors.write(
                f"{chrom}\t"
                f"{anchor0}\t"
                f"{anchor0 + 1}\t"
                f"{locus_id}\t"
                f"0\t.\n"
            )

            out = {
                "canonical_locus_id": locus_id,
                "locus_name": row.get("locus_name", ""),
                "benchmark_role": row["benchmark_role"],
                "best_evidence_tier": row.get(
                    "best_evidence_tier", ""
                ),
                "target_assembly": row["target_assembly"],
                "target_seqname": chrom,
                "target_start_1based": start1,
                "target_end_1based": end1,
                "locus_length_bp": end1 - start1 + 1,
                "anchor_1based": anchor1,
            }

            add_peak_features(
                out,
                "atac_idr10",
                by_chrom10.get(chrom, []),
                start0,
                end0,
                anchor0,
                windows,
            )

            add_peak_features(
                out,
                "atac_idr05",
                by_chrom05.get(chrom, []),
                start0,
                end0,
                anchor0,
                windows,
            )

            output_rows.append(out)

    fieldnames = list(output_rows[0].keys())

    with open(
        args.features,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(output_rows)

    roles = Counter(
        row["benchmark_role"]
        for row in benchmark
    )

    summary = [
        ("benchmark_loci", len(benchmark)),
        ("benchmark_positive", roles.get("positive", 0)),
        ("benchmark_negative", roles.get("negative", 0)),
        ("idr10_raw_records", raw10),
        ("idr10_unique_intervals", len(unique10)),
        ("idr10_exact_duplicate_records", raw10 - len(unique10)),
        ("idr05_raw_records", raw05),
        ("idr05_unique_intervals", len(unique05)),
        ("idr05_exact_duplicate_records", raw05 - len(unique05)),
        (
            "gold_loci_anchor_in_idr10_peak",
            sum(
                int(row["atac_idr10_anchor_in_peak"])
                for row in output_rows
            ),
        ),
        (
            "gold_loci_anchor_in_idr05_peak",
            sum(
                int(row["atac_idr05_anchor_in_peak"])
                for row in output_rows
            ),
        ),
    ]

    with open(
        args.summary,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(["metric", "value"])
        writer.writerows(summary)

    provenance = {
        "benchmark_tsv": {
            "path": str(
                Path(args.benchmark_tsv).resolve()
            ),
            "sha256": sha256(args.benchmark_tsv),
        },
        "benchmark_bed": {
            "path": str(
                Path(args.benchmark_bed).resolve()
            ),
            "sha256": sha256(args.benchmark_bed),
        },
        "idr10": {
            "path": str(Path(args.idr10).resolve()),
            "sha256": sha256(args.idr10),
            "raw_records": raw10,
            "unique_intervals": len(unique10),
        },
        "idr05": {
            "path": str(Path(args.idr05).resolve()),
            "sha256": sha256(args.idr05),
            "raw_records": raw05,
            "unique_intervals": len(unique05),
        },
        "anchor_definition": (
            "integer centre of 1-based closed benchmark interval"
        ),
        "anchor_half_windows_bp": windows,
        "coordinate_system": {
            "benchmark_tsv": "1-based closed",
            "internal_and_bed": "0-based half-open",
        },
    }

    Path(args.provenance).write_text(
        json.dumps(
            provenance,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        f"Benchmark loci: {len(benchmark)}"
    )
    print(
        f"IDR10 raw/unique: {raw10}/{len(unique10)}"
    )
    print(
        f"IDR05 raw/unique: {raw05}/{len(unique05)}"
    )
    print(
        f"Features written: {args.features}"
    )


if __name__ == "__main__":
    main()
