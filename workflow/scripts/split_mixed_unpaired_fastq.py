#!/usr/bin/env python3

import argparse
import csv
import gzip
from pathlib import Path


def classify_header(header):
    """
    Classify a FASTQ header as originating from R1 or R2.

    The complete header line is inspected because ENA/SRA FASTQ
    records may encode the mate identifier in the description field,
    for example:

        @SRRXXXX.12345 12345/1
        @SRRXXXX.12346 12346/2
    """

    header = header.rstrip(b"\r\n")

    if header.endswith(b"/1"):
        return "R1"

    if header.endswith(b"/2"):
        return "R2"

    return "OTHER"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Split a MIXED_R1_R2 gzipped UNPAIRED FASTQ into "
            "separate R1- and R2-derived FASTQ files."
        )
    )

    parser.add_argument(
        "--input-fastq",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--run-accession",
        required=True,
    )

    parser.add_argument(
        "--r1-output",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--r2-output",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
    )

    args = parser.parse_args()

    if not args.input_fastq.is_file():
        raise FileNotFoundError(
            f"Input FASTQ not found: {args.input_fastq}"
        )

    for output in (
        args.r1_output,
        args.r2_output,
        args.summary,
    ):
        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    r1_tmp = Path(str(args.r1_output) + ".tmp")
    r2_tmp = Path(str(args.r2_output) + ".tmp")
    summary_tmp = Path(str(args.summary) + ".tmp")

    for tmp in (r1_tmp, r2_tmp, summary_tmp):
        tmp.unlink(missing_ok=True)

    total = 0
    r1_records = 0
    r2_records = 0
    other_records = 0

    try:
        with gzip.open(args.input_fastq, "rb") as source, \
             gzip.open(r1_tmp, "wb", compresslevel=6) as r1_handle, \
             gzip.open(r2_tmp, "wb", compresslevel=6) as r2_handle:

            while True:
                header = source.readline()

                if not header:
                    break

                sequence = source.readline()
                plus = source.readline()
                quality = source.readline()

                if not sequence or not plus or not quality:
                    raise ValueError(
                        f"Incomplete FASTQ record in "
                        f"{args.input_fastq}"
                    )

                if not header.startswith(b"@"):
                    raise ValueError(
                        f"Invalid FASTQ header in "
                        f"{args.input_fastq}: "
                        f"{header!r}"
                    )

                if not plus.startswith(b"+"):
                    raise ValueError(
                        f"Invalid FASTQ separator in "
                        f"{args.input_fastq}"
                    )

                record = (
                    header
                    + sequence
                    + plus
                    + quality
                )

                mate = classify_header(header)

                total += 1

                if mate == "R1":
                    r1_handle.write(record)
                    r1_records += 1

                elif mate == "R2":
                    r2_handle.write(record)
                    r2_records += 1

                else:
                    other_records += 1

        if total == 0:
            raise ValueError(
                f"Input FASTQ is empty: {args.input_fastq}"
            )

        if other_records > 0:
            raise ValueError(
                f"Cannot safely split {args.run_accession}: "
                f"{other_records} reads have unrecognized "
                "mate provenance."
            )

        if r1_records == 0 or r2_records == 0:
            raise ValueError(
                f"{args.run_accession} is not MIXED_R1_R2: "
                f"R1={r1_records}, R2={r2_records}."
            )

        classification = "MIXED_R1_R2"

        row = {
            "run_accession": args.run_accession,
            "input_fastq": str(args.input_fastq),
            "total_records": total,
            "r1_records": r1_records,
            "r2_records": r2_records,
            "other_records": other_records,
            "r1_fraction": r1_records / total,
            "r2_fraction": r2_records / total,
            "mate_provenance": classification,
        }

        with summary_tmp.open(
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

        r1_tmp.replace(args.r1_output)
        r2_tmp.replace(args.r2_output)
        summary_tmp.replace(args.summary)

    except Exception:
        r1_tmp.unlink(missing_ok=True)
        r2_tmp.unlink(missing_ok=True)
        summary_tmp.unlink(missing_ok=True)
        raise

    print(f"Run: {args.run_accession}")
    print(f"Total records: {total}")
    print(f"R1 records: {r1_records}")
    print(f"R2 records: {r2_records}")
    print(f"Other records: {other_records}")
    print("Mate provenance: MIXED_R1_R2")
    print(f"R1 output: {args.r1_output}")
    print(f"R2 output: {args.r2_output}")
    print(f"Summary: {args.summary}")


if __name__ == "__main__":
    main()
