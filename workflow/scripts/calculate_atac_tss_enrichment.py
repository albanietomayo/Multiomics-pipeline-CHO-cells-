#!/usr/bin/env python3

import argparse
import bisect
from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam


def load_tss_bed(path, reference_lengths, window):
    """
    Load unique TSS positions from BED6.

    TSS too close to sequence boundaries are excluded because a complete
    +/- window cannot be evaluated around them.
    """
    positions = defaultdict(list)
    strands = defaultdict(list)

    total_tss = 0
    usable_tss = 0
    skipped_boundary = 0
    skipped_missing_reference = 0

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 6:
                raise ValueError("TSS BED must contain at least 6 columns.")

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            strand = fields[5]

            total_tss += 1

            if end - start != 1:
                raise ValueError(
                    f"TSS record is not 1 bp: {chrom}:{start}-{end}"
                )

            if strand not in {"+", "-"}:
                raise ValueError(
                    f"Unsupported TSS strand '{strand}' at {chrom}:{start}"
                )

            if chrom not in reference_lengths:
                skipped_missing_reference += 1
                continue

            tss = start
            chrom_length = reference_lengths[chrom]

            if tss - window < 0 or tss + window > chrom_length:
                skipped_boundary += 1
                continue

            positions[chrom].append(tss)
            strands[chrom].append(strand)
            usable_tss += 1

    for chrom in positions:
        ordered = sorted(zip(positions[chrom], strands[chrom]))
        positions[chrom] = [x[0] for x in ordered]
        strands[chrom] = [x[1] for x in ordered]

    return (
        positions,
        strands,
        total_tss,
        usable_tss,
        skipped_boundary,
        skipped_missing_reference,
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--bam", required=True)
    parser.add_argument("--tss-bed", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--summary", required=True)

    parser.add_argument("--window", type=int, required=True)
    parser.add_argument("--bin-size", type=int, required=True)
    parser.add_argument("--background", type=int, required=True)

    parser.add_argument("--forward-shift", type=int, required=True)
    parser.add_argument("--reverse-shift", type=int, required=True)
    parser.add_argument("--min-mapq", type=int, required=True)

    args = parser.parse_args()

    if (2 * args.window) % args.bin_size != 0:
        raise ValueError(
            "The full TSS window must be divisible by bin size."
        )

    if args.background % args.bin_size != 0:
        raise ValueError(
            "Background size must be divisible by bin size."
        )

    n_bins = (2 * args.window) // args.bin_size
    background_bins = args.background // args.bin_size

    if 2 * background_bins >= n_bins:
        raise ValueError("Background regions occupy the complete profile.")

    profile = np.zeros(n_bins, dtype=np.int64)

    profile_path = Path(args.profile)
    summary_path = Path(args.summary)

    profile_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        reference_lengths = dict(zip(bam.references, bam.lengths))

        (
            tss_positions,
            tss_strands,
            total_tss,
            usable_tss,
            skipped_boundary,
            skipped_missing_reference,
        ) = load_tss_bed(
            args.tss_bed,
            reference_lengths,
            args.window,
        )

        reads_seen = 0
        reads_used = 0
        insertion_events = 0

        for read in bam.fetch(until_eof=True):
            reads_seen += 1

            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.is_duplicate
                or read.is_qcfail
                or read.mapping_quality < args.min_mapq
            ):
                continue

            chrom = bam.get_reference_name(read.reference_id)

            if chrom not in tss_positions:
                continue

            if read.is_reverse:
                if read.reference_end is None:
                    continue
                insertion = read.reference_end + args.reverse_shift
            else:
                insertion = read.reference_start + args.forward_shift

            if insertion < 0 or insertion >= reference_lengths[chrom]:
                continue

            reads_used += 1

            positions = tss_positions[chrom]
            strands = tss_strands[chrom]

            left = bisect.bisect_left(
                positions,
                insertion - args.window + 1,
            )
            right = bisect.bisect_right(
                positions,
                insertion + args.window,
            )

            for idx in range(left, right):
                tss = positions[idx]
                strand = strands[idx]

                if strand == "+":
                    relative = insertion - tss
                else:
                    relative = tss - insertion

                if -args.window <= relative < args.window:
                    bin_index = (
                        relative + args.window
                    ) // args.bin_size

                    profile[bin_index] += 1
                    insertion_events += 1

    left_background = profile[:background_bins]
    right_background = profile[-background_bins:]

    background_values = np.concatenate(
        [left_background, right_background]
    )

    background_mean = float(np.mean(background_values))

    if background_mean <= 0:
        raise ValueError(
            "Background signal is zero; TSS enrichment cannot be normalized."
        )

    normalized = profile.astype(float) / background_mean
    tss_enrichment_score = float(np.max(normalized))

    with open(profile_path, "w", encoding="utf-8") as handle:
        handle.write(
            "relative_start\trelative_end\t"
            "insertion_count\tnormalized_signal\n"
        )

        for i in range(n_bins):
            relative_start = -args.window + i * args.bin_size
            relative_end = relative_start + args.bin_size

            handle.write(
                f"{relative_start}\t"
                f"{relative_end}\t"
                f"{int(profile[i])}\t"
                f"{normalized[i]:.6f}\n"
            )

    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"total_tss\t{total_tss}\n")
        handle.write(f"usable_tss\t{usable_tss}\n")
        handle.write(
            f"tss_skipped_boundary\t{skipped_boundary}\n"
        )
        handle.write(
            "tss_skipped_missing_reference\t"
            f"{skipped_missing_reference}\n"
        )
        handle.write(f"reads_seen\t{reads_seen}\n")
        handle.write(f"reads_used\t{reads_used}\n")
        handle.write(
            f"insertion_events_in_tss_windows\t{insertion_events}\n"
        )
        handle.write(f"window_bp\t{args.window}\n")
        handle.write(f"bin_size_bp\t{args.bin_size}\n")
        handle.write(f"background_bp_per_side\t{args.background}\n")
        handle.write(
            f"background_mean\t{background_mean:.6f}\n"
        )
        handle.write(
            f"tss_enrichment_score\t{tss_enrichment_score:.6f}\n"
        )


if __name__ == "__main__":
    main()
