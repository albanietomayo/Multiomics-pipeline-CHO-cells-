#!/usr/bin/env python3

import argparse
import gzip
import re


def normalize_read_id(header):
    """Return FASTQ read identifier without @ and terminal /1 or /2."""
    token = header.strip().split()[0]
    if token.startswith("@"):
        token = token[1:]
    return re.sub(r"/[12]$", "", token)


def read_ids(path):
    """Yield normalized read IDs from a gzipped FASTQ."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        record = []
        for line in handle:
            record.append(line)
            if len(record) == 4:
                yield normalize_read_id(record[0])
                record = []

        if record:
            raise ValueError(f"Incomplete FASTQ record in {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--unpaired", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    print("Loading paired read identifiers...")

    r1_ids = set(read_ids(args.r1))
    r2_ids = set(read_ids(args.r2))
    paired_ids = r1_ids | r2_ids

    print("Checking unpaired reads...")

    unpaired_total = 0
    overlap_r1 = 0
    overlap_r2 = 0
    overlap_any = 0

    for read_id in read_ids(args.unpaired):
        unpaired_total += 1

        in_r1 = read_id in r1_ids
        in_r2 = read_id in r2_ids

        overlap_r1 += int(in_r1)
        overlap_r2 += int(in_r2)
        overlap_any += int(read_id in paired_ids)

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("metric\tvalue\n")
        out.write(f"r1_unique_ids\t{len(r1_ids)}\n")
        out.write(f"r2_unique_ids\t{len(r2_ids)}\n")
        out.write(f"paired_union_unique_ids\t{len(paired_ids)}\n")
        out.write(f"unpaired_records\t{unpaired_total}\n")
        out.write(f"unpaired_overlap_r1\t{overlap_r1}\n")
        out.write(f"unpaired_overlap_r2\t{overlap_r2}\n")
        out.write(f"unpaired_overlap_any_paired\t{overlap_any}\n")

    print(f"UNPAIRED records: {unpaired_total}")
    print(f"Overlap with R1: {overlap_r1}")
    print(f"Overlap with R2: {overlap_r2}")
    print(f"Overlap with paired union: {overlap_any}")


if __name__ == "__main__":
    main()
