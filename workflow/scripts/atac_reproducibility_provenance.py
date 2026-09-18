#!/usr/bin/env python3
"""Collect group-level ATAC reproducibility provenance with streaming hashes."""

import argparse
import json
from pathlib import Path

from file_hash import sha256_file


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", required=True)
    parser.add_argument("--eligibility", required=True)
    parser.add_argument("--input", action="append", nargs="+", default=[])
    parser.add_argument("--parameters-json", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    eligibility = json.loads(Path(args.eligibility).read_text(encoding="utf-8"))
    if not eligibility.get("eligible") or eligibility.get("replicate_group") != args.group:
        raise ValueError("Replicate-group eligibility record is invalid")
    inputs = [path for group in args.input for path in group]
    value = {
        "schema_version": 1,
        "replicate_group": args.group,
        "members": eligibility.get("members", []),
        "replicate_type": eligibility.get("replicate_type"),
        "eligibility": {
            "path": str(Path(args.eligibility)),
            "sha256": sha256_file(args.eligibility),
        },
        "parameters": json.loads(args.parameters_json),
        "artifacts": {
            Path(path).name: {"path": str(Path(path)), "sha256": sha256_file(path)}
            for path in inputs
        },
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
