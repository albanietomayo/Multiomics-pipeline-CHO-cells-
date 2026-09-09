#!/usr/bin/env python3

import argparse
from pathlib import Path


def parse_transcript_tss(gtf_path):
    """
    Extract unique transcript-level TSS coordinates from a GTF file.

    GTF coordinates are 1-based inclusive.
    BED coordinates are 0-based half-open.
    """
    transcript_count = 0
    tss_records = set()

    with open(gtf_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                raise ValueError(
                    f"Expected 9 GTF columns, found {len(fields)}."
                )

            chrom, _, feature, start, end, _, strand, _, _ = fields

            if feature != "transcript":
                continue

            transcript_count += 1

            start = int(start)
            end = int(end)

            if strand == "+":
                bed_start = start - 1
                bed_end = start
            elif strand == "-":
                bed_start = end - 1
                bed_end = end
            else:
                raise ValueError(
                    f"Unsupported strand '{strand}' for transcript on {chrom}."
                )

            if bed_start < 0:
                raise ValueError(
                    f"Invalid BED coordinate generated for {chrom}: "
                    f"{bed_start}-{bed_end}"
                )

            tss_records.add((chrom, bed_start, bed_end, strand))

    if transcript_count == 0:
        raise ValueError("No transcript features were found in the GTF.")

    return transcript_count, tss_records


def write_bed(tss_records, output_path):
    """
    Write unique TSS positions as BED6.
    """
    ordered = sorted(
        tss_records,
        key=lambda record: (record[0], record[1], record[2], record[3]),
    )

    with open(output_path, "w", encoding="utf-8") as handle:
        for index, (chrom, start, end, strand) in enumerate(
            ordered, start=1
        ):
            name = f"TSS_{index:06d}"
            handle.write(
                f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n"
            )


def write_summary(transcript_count, unique_tss_count, output_path):
    """
    Write structured provenance statistics for the TSS resource.
    """
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"transcript_features\t{transcript_count}\n")
        handle.write(f"unique_tss\t{unique_tss_count}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    gtf_path = Path(args.gtf)
    bed_path = Path(args.bed)
    summary_path = Path(args.summary)

    bed_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    transcript_count, tss_records = parse_transcript_tss(gtf_path)

    write_bed(tss_records, bed_path)
    write_summary(
        transcript_count,
        len(tss_records),
        summary_path,
    )


if __name__ == "__main__":
    main()
