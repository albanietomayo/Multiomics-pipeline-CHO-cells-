#!/usr/bin/env python3
"""Validate the ATAC BigWig header against reference and retained signal."""

import argparse
from pathlib import Path

from clip_atac_bedgraph import read_sizes
from file_hash import sha256_file


def signal_chromosomes(path):
    chromosomes = set()
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("track", "browser", "#")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"Malformed clipped bedGraph line {number}")
            try:
                start, end = int(fields[1]), int(fields[2])
                float(fields[3])
            except ValueError as exc:
                raise ValueError(f"Invalid clipped bedGraph line {number}") from exc
            if start < 0 or start >= end:
                raise ValueError(f"Invalid clipped bedGraph coordinates at line {number}")
            chromosomes.add(fields[0])
    if not chromosomes:
        raise ValueError("Clipped bedGraph contains no retained signal")
    return chromosomes


def validate_header(header, sizes, signal):
    """Enforce the pilot contract: header chromosomes equal signal chromosomes."""
    header = dict(header)
    header_chromosomes = set(header)
    unknown = header_chromosomes - set(sizes)
    if unknown:
        raise ValueError("BigWig contains chromosomes outside the reference: " + ",".join(sorted(unknown)))
    mismatched = sorted(c for c in header_chromosomes if header[c] != sizes[c])
    if mismatched:
        raise ValueError("BigWig chromosome length mismatch: " + ",".join(mismatched))
    missing = set(signal) - header_chromosomes
    if missing:
        raise ValueError("BigWig is missing signal chromosomes: " + ",".join(sorted(missing)))
    extra = header_chromosomes - set(signal)
    if extra:
        raise ValueError("BigWig header contains chromosomes without retained signal: " + ",".join(sorted(extra)))


def validate_bigwig(bigwig, chrom_sizes, clipped_bedgraph):
    import pyBigWig

    sizes = read_sizes(chrom_sizes)
    signal = signal_chromosomes(clipped_bedgraph)
    unknown_signal = signal - set(sizes)
    if unknown_signal:
        raise ValueError("Clipped signal contains chromosomes outside the reference: " + ",".join(sorted(unknown_signal)))
    with pyBigWig.open(str(bigwig)) as handle:
        header = handle.chroms()
    if not header:
        raise ValueError("BigWig has an empty chromosome header")
    validate_header(header, sizes, signal)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bigwig", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--clipped-bedgraph", required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()
    validate_bigwig(args.bigwig, args.chrom_sizes, args.clipped_bedgraph)
    Path(args.sha256).write_text(
        f"{sha256_file(args.bigwig)}  {Path(args.bigwig).name}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
