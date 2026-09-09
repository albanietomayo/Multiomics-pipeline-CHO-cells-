#!/usr/bin/env python3

import argparse
import csv
import gzip
from pathlib import Path


def classify_header(header):
    """
    Classify one FASTQ header as originating from mate R1, mate R2,
    or as having no recognizable mate suffix.
    """

    header = header.strip()

    if header.endswith("/1"):
        return "R1"

    if header.endswith("/2"):
        return "R2"

    return "OTHER"


def inspect_fastq(fastq_file):
    """
    Inspect all records in one local gzipped FASTQ file.

    Returns complete mate-provenance counts. This intentionally scans
    the entire validated local FASTQ rather than a sample.
    """

    counts = {
        "R1": 0,
        "R2": 0,
        "OTHER": 0,
    }

    records = 0

    with gzip.open(
        fastq_file,
        "rt",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        while True:
            header = handle.readline()

            if not header:
                break

            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            if not sequence or not plus or not quality:
                raise ValueError(
                    f"Incomplete FASTQ record detected in {fastq_file}"
                )

            if not header.startswith("@"):
                raise ValueError(
                    f"Invalid FASTQ header in {fastq_file}: "
                    f"{header.rstrip()}"
                )

            mate = classify_header(header)

            counts[mate] += 1
            records += 1

    return records, counts


def determine_provenance(records, counts):
    """
    Convert complete mate counts into a categorical provenance label.
    """

    if records == 0:
        return "EMPTY"

    r1 = counts["R1"]
    r2 = counts["R2"]
    other = counts["OTHER"]

    if r1 > 0 and r2 == 0 and other == 0:
        return "R1_ONLY"

    if r2 > 0 and r1 == 0 and other == 0:
        return "R2_ONLY"

    if r1 > 0 and r2 > 0 and other == 0:
        return "MIXED_R1_R2"

    return "UNKNOWN"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Inspect a local gzipped UNPAIRED FASTQ and determine "
            "whether its reads originate from R1, R2, both, or an "
            "unrecognized header convention."
        )
    )

    parser.add_argument(
        "--fastq",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--run-accession",
        required=True,
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )

    args = parser.parse_args()

    if not args.fastq.is_file():
        raise FileNotFoundError(
            f"FASTQ not found: {args.fastq}"
        )

    records, counts = inspect_fastq(args.fastq)

    provenance = determine_provenance(
        records,
        counts,
    )

    row = {
        "run_accession": args.run_accession,
        "fastq": str(args.fastq),
        "records": records,
        "r1_reads": counts["R1"],
        "r2_reads": counts["R2"],
        "other_reads": counts["OTHER"],
        "r1_fraction": (
            counts["R1"] / records
            if records else 0.0
        ),
        "r2_fraction": (
            counts["R2"] / records
            if records else 0.0
        ),
        "other_fraction": (
            counts["OTHER"] / records
            if records else 0.0
        ),
        "mate_provenance": provenance,
    }

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with args.output.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=row.keys(),
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerow(row)

    print(f"Run: {args.run_accession}")
    print(f"FASTQ: {args.fastq}")
    print(f"Records: {records}")
    print(f"R1 reads: {counts['R1']}")
    print(f"R2 reads: {counts['R2']}")
    print(f"Other reads: {counts['OTHER']}")
    print(f"Mate provenance: {provenance}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
