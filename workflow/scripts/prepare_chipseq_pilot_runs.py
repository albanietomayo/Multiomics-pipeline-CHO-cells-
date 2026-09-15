#!/usr/bin/env python3
"""Prepare or verify the FASTQ selection for the documented ChIP pilot."""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--output", type=Path)
    action.add_argument("--check-manifest", type=Path)
    args = parser.parse_args()
    pilot = json.loads(args.pilot.read_text())
    actual_sha = hashlib.sha256(args.samples.read_bytes()).hexdigest()
    if actual_sha != pilot["provenance_sha256"]["samples"]:
        raise ValueError("Sample catalog differs from the documented pilot")
    if pilot["schema_version"] != 1:
        raise ValueError("Unsupported pilot schema")
    expected = {pilot["ip_run_accession"]: "ip", pilot["input_run_accession"]: "input"}
    declared = pilot["runs"]
    if (len(expected) != 2 or len(declared) != 2 or
            {r["run_accession"]: r["role"] for r in declared} != expected):
        raise ValueError("The pilot must contain exactly one IP and one distinct input")
    with args.samples.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    samples = {r["run_accession"]: r for r in rows}
    if len(samples) != len(rows):
        raise ValueError("Duplicate run accessions in sample catalog")
    physical = {}
    for item in declared:
        run = item["run_accession"]
        row = samples[run]
        if (not re.fullmatch(r"[DES]RR[0-9]+", run) or
                row["omics"] != "ChIP-seq" or row["library_layout"] != "SINGLE" or
                row["instrument_platform"] != "ILLUMINA" or
                row["study_accession"] != pilot["study_accession"]):
            raise ValueError(f"Unexpected pilot metadata: {run}")
        urls, sizes, md5s = [row[k].split(";") for k in ("fastq_ftp", "fastq_bytes", "fastq_md5")]
        if len(urls) != 1 or len(sizes) != 1 or len(md5s) != 1:
            raise ValueError(f"Expected one physical FASTQ: {run}")
        filename = urls[0].rsplit("/", 1)[-1]
        size = int(sizes[0])
        if (not re.fullmatch(r"[A-Za-z0-9_.-]+[.]fastq[.]gz", filename) or
                not re.fullmatch(r"[a-fA-F0-9]{32}", md5s[0]) or
                size <= 0 or size != item["fastq_bytes"]):
            raise ValueError(f"Invalid FASTQ metadata: {run}")
        physical[run] = (filename, size)
    if sum(size for _, size in physical.values()) != pilot["combined_fastq_bytes"]:
        raise ValueError("Total FASTQ size differs from documented pilot")
    if args.check_manifest:
        with args.check_manifest.open(encoding="utf-8", newline="") as handle:
            manifest = list(csv.DictReader(handle, delimiter="\t"))
        if len(manifest) != 2 or {r["run_accession"] for r in manifest} != set(expected):
            raise ValueError("Manifest must contain exactly the two pilot runs")
        for row in manifest:
            filename, size = physical[row["run_accession"]]
            if (row["fastq_role"] != "SINGLE" or row["filename"] != filename or
                    int(row["fastq_bytes"]) != size):
                raise ValueError("Unexpected physical FASTQ in pilot manifest")
        print("[OK] Manifest verified: two SINGLE FASTQ files, only the documented pilot.")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text("run_accession\n" + "".join(r + "\n" for r in sorted(expected)))
        print("[OK] Pilot selection generated.")


if __name__ == "__main__":
    main()
