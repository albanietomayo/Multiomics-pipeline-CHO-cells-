#!/usr/bin/env python3
"""Verify every planned ENA record and write a portable XML inventory."""
import argparse
import csv
import io
import json
from pathlib import Path
from fetch_chipseq_ena_xml import digest, validate_xml, write_atomic


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--snapshot", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    from chipseq_metadata_snapshot import summarize
    summarize(args.plan, args.snapshot, args.output)
    print("[OK] Verified sanitized scientific metadata; historical inventory preserved")


if __name__ == "__main__":
    main()
