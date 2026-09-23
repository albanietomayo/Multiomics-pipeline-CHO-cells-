#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

import pyBigWig


def sha256(path):
    h = hashlib.sha256()

    with open(path, "rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


def summarize_window(bw, chrom, start, end):
    """
    Summarise BigWig signal across the COMPLETE requested interval.

    BigWig positions without an explicitly stored interval are treated
    as zero signal. Coordinates are 0-based half-open.
    """

    chroms = bw.chroms()

    if chrom not in chroms:
        raise ValueError(
            f"Contig absent from BigWig: {chrom}"
        )

    chrom_length = chroms[chrom]

    clipped_start = max(0, start)
    clipped_end = min(chrom_length, end)

    if clipped_end <= clipped_start:
        raise ValueError(
            f"Empty interval after clipping: "
            f"{chrom}:{start}-{end}"
        )

    window_bp = clipped_end - clipped_start

    intervals = bw.intervals(
        chrom,
        clipped_start,
        clipped_end
    )

    if intervals is None:
        intervals = []

    weighted_sum = 0.0
    covered_bp = 0
    max_signal = 0.0

    for interval_start, interval_end, value in intervals:

        value = float(value)

        if not math.isfinite(value):
            raise ValueError(
                f"Non-finite BigWig value at "
                f"{chrom}:{interval_start}-{interval_end}"
            )

        if value < 0:
            raise ValueError(
                f"Unexpected negative SPMR signal at "
                f"{chrom}:{interval_start}-{interval_end}: "
                f"{value}"
            )

        ov_start = max(
            clipped_start,
            interval_start
        )

        ov_end = min(
            clipped_end,
            interval_end
        )

        if ov_end <= ov_start:
            continue

        bp = ov_end - ov_start

        weighted_sum += value * bp

        if value > 0:
            covered_bp += bp

        if value > max_signal:
            max_signal = value

    mean_signal = (
        weighted_sum / window_bp
    )

    return {
        "mean": mean_signal,
        "max": max_signal,
        "covered_bp": covered_bp,
        "covered_fraction": (
            covered_bp / window_bp
        ),
        "actual_window_bp": window_bp,
    }


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--features",
        required=True
    )

    parser.add_argument(
        "--rep1",
        required=True
    )

    parser.add_argument(
        "--rep2",
        required=True
    )

    parser.add_argument(
        "--rep1-name",
        default="SRR12774931"
    )

    parser.add_argument(
        "--rep2-name",
        default="SRR12774932"
    )

    parser.add_argument(
        "--windows",
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )

    parser.add_argument(
        "--summary",
        required=True
    )

    parser.add_argument(
        "--provenance",
        required=True
    )

    args = parser.parse_args()

    windows = sorted({
        int(x)
        for x in args.windows.split(",")
        if x.strip()
    })

    if not windows:
        raise ValueError(
            "No windows supplied"
        )

    if any(x <= 0 for x in windows):
        raise ValueError(
            "Window half-sizes must be positive"
        )

    with open(
        args.features,
        encoding="utf-8",
        newline=""
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        rows = list(reader)
        original_fields = list(
            reader.fieldnames or []
        )

    required = {
        "canonical_locus_id",
        "benchmark_role",
        "target_seqname",
        "anchor_1based",
    }

    missing = required - set(
        original_fields
    )

    if missing:
        raise ValueError(
            f"Missing feature columns: "
            f"{sorted(missing)}"
        )

    if not rows:
        raise ValueError(
            "Feature table contains no benchmark loci"
        )

    locus_ids = [
        row["canonical_locus_id"]
        for row in rows
    ]

    if len(locus_ids) != len(set(locus_ids)):
        raise ValueError(
            "Duplicate canonical_locus_id"
        )

    bw1 = pyBigWig.open(args.rep1)
    bw2 = pyBigWig.open(args.rep2)

    if not bw1.isBigWig():
        raise ValueError(
            f"Not a BigWig: {args.rep1}"
        )

    if not bw2.isBigWig():
        raise ValueError(
            f"Not a BigWig: {args.rep2}"
        )

    try:

        benchmark_contigs = {
            row["target_seqname"]
            for row in rows
        }

        for name, bw in (
            (args.rep1_name, bw1),
            (args.rep2_name, bw2),
        ):

            missing_contigs = (
                benchmark_contigs
                - set(bw.chroms())
            )

            if missing_contigs:
                raise ValueError(
                    f"{name}: benchmark contigs "
                    f"missing from BigWig: "
                    f"{sorted(missing_contigs)}"
                )

        for row in rows:

            chrom = row["target_seqname"]

            # anchor_1based is a 1-based genomic position.
            anchor0 = (
                int(row["anchor_1based"])
                - 1
            )

            for half_window in windows:

                # Inclusive biological interpretation:
                # anchor +/- N bp -> 2N+1 positions.
                start0 = (
                    anchor0
                    - half_window
                )

                end0 = (
                    anchor0
                    + half_window
                    + 1
                )

                s1 = summarize_window(
                    bw1,
                    chrom,
                    start0,
                    end0
                )

                s2 = summarize_window(
                    bw2,
                    chrom,
                    start0,
                    end0
                )

                label = (
                    f"pm{half_window}bp"
                )

                for rep_name, summary in (
                    (args.rep1_name, s1),
                    (args.rep2_name, s2),
                ):

                    prefix = (
                        f"atac_spmr_"
                        f"{rep_name}_"
                        f"{label}"
                    )

                    row[
                        f"{prefix}_mean"
                    ] = summary["mean"]

                    row[
                        f"{prefix}_max"
                    ] = summary["max"]

                    row[
                        f"{prefix}_covered_bp"
                    ] = summary["covered_bp"]

                    row[
                        f"{prefix}_covered_fraction"
                    ] = summary["covered_fraction"]

                    row[
                        f"{prefix}_actual_window_bp"
                    ] = summary["actual_window_bp"]

                consensus_prefix = (
                    f"atac_spmr_consensus_"
                    f"{label}"
                )

                row[
                    f"{consensus_prefix}_mean"
                ] = statistics.mean([
                    s1["mean"],
                    s2["mean"]
                ])

                row[
                    f"{consensus_prefix}_mean_of_rep_max"
                ] = statistics.mean([
                    s1["max"],
                    s2["max"]
                ])

                row[
                    f"{consensus_prefix}_mean_covered_fraction"
                ] = statistics.mean([
                    s1["covered_fraction"],
                    s2["covered_fraction"]
                ])

    finally:
        bw1.close()
        bw2.close()

    Path(args.output).parent.mkdir(
        parents=True,
        exist_ok=True
    )

    fieldnames = list(rows[0].keys())

    with open(
        args.output,
        "w",
        encoding="utf-8",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writeheader()
        writer.writerows(rows)

    summary_rows = []

    for half_window in windows:

        col = (
            "atac_spmr_consensus_"
            f"pm{half_window}bp_mean"
        )

        for role in sorted({
            row["benchmark_role"]
            for row in rows
        }):

            values = [
                float(row[col])
                for row in rows
                if row["benchmark_role"] == role
            ]

            summary_rows.append({
                "half_window_bp":
                    half_window,
                "full_window_bp":
                    2 * half_window + 1,
                "benchmark_role":
                    role,
                "n":
                    len(values),
                "median_consensus_mean_spmr":
                    statistics.median(values),
                "min_consensus_mean_spmr":
                    min(values),
                "max_consensus_mean_spmr":
                    max(values),
            })

    summary_fields = [
        "half_window_bp",
        "full_window_bp",
        "benchmark_role",
        "n",
        "median_consensus_mean_spmr",
        "min_consensus_mean_spmr",
        "max_consensus_mean_spmr",
    ]

    with open(
        args.summary,
        "w",
        encoding="utf-8",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=summary_fields,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writeheader()
        writer.writerows(
            summary_rows
        )

    provenance = {
        "pybigwig_version":
            pyBigWig.__version__,

        "input_features": {
            "path":
                str(
                    Path(
                        args.features
                    ).resolve()
                ),
            "sha256":
                sha256(
                    args.features
                ),
        },

        "replicate_bigwigs": {
            args.rep1_name: {
                "path":
                    str(
                        Path(
                            args.rep1
                        ).resolve()
                    ),
                "sha256":
                    sha256(
                        args.rep1
                    ),
            },

            args.rep2_name: {
                "path":
                    str(
                        Path(
                            args.rep2
                        ).resolve()
                    ),
                "sha256":
                    sha256(
                        args.rep2
                    ),
            },
        },

        "anchor_half_windows_bp":
            windows,

        "window_definition":
            "anchor +/- N bp; 2N+1 genomic bases before chromosome-boundary clipping",

        "missing_bigwig_positions":
            "treated as zero signal",

        "consensus_signal":
            "arithmetic mean of replicate SPMR window means",

        "n_loci":
            len(rows),
    }

    Path(
        args.provenance
    ).write_text(
        json.dumps(
            provenance,
            indent=2,
            sort_keys=True
        ) + "\n",
        encoding="utf-8"
    )

    print(
        "Loci processed =",
        len(rows)
    )

    print(
        "Windows =",
        ",".join(
            str(x)
            for x in windows
        )
    )

    print(
        "Output =",
        args.output
    )


if __name__ == "__main__":
    main()
