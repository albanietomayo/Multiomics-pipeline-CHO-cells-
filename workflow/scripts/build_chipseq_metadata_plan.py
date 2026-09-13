#!/usr/bin/env python3
"""Build a deterministic ENA metadata request plan for ChIP-seq."""

import argparse
import csv
import hashlib
import io
import json
import os
import re
import tempfile
from pathlib import Path


def atomic_write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="",
            dir=path.parent, delete=False
        ) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def tsv_text(fields, rows):
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    raw = args.samples.read_bytes()
    reader = csv.DictReader(
        io.StringIO(raw.decode("utf-8-sig"), newline=""),
        delimiter="\t"
    )

    fields = [
        "run_accession", "study_accession", "experiment_accession",
        "sample_accession", "omics", "experiment_title",
        "sample_title", "library_name", "library_layout",
        "instrument_platform"
    ]
    missing = sorted(set(fields) - set(reader.fieldnames or []))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    runs = []
    seen_runs = set()
    requests = set()

    for row in reader:
        if (row["omics"] or "").strip().casefold() != "chip-seq":
            continue
        record = {field: (row[field] or "").strip() for field in fields}
        run = record["run_accession"]

        if not re.fullmatch(r"[SED]RR[0-9]+", run):
            raise ValueError(f"Invalid run accession: {run!r}")
        if run in seen_runs:
            raise ValueError(f"Duplicated run accession: {run}")
        seen_runs.add(run)

        for kind, field, pattern in (
            ("experiment", "experiment_accession", r"[SED]RX[0-9]+"),
            ("sample", "sample_accession", r"(?:SAM[END][A-Z]*|[SED]RS)[0-9]+"),
        ):
            accession = record[field]
            if not re.fullmatch(pattern, accession):
                raise ValueError(
                    f"{run}: missing or invalid {field}: {accession!r}"
                )
            requests.add((kind, accession))

        runs.append(record)

    if not runs:
        raise ValueError("No ChIP-seq runs found; no outputs written.")

    runs.sort(key=lambda row: row["run_accession"])
    request_rows = [
        {
            "record_type": kind,
            "accession": accession,
            "source_url": (
                f"https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
            ),
        }
        for kind, accession in sorted(requests)
    ]

    provenance = {
        "schema_version": 1,
        "samples_sha256": hashlib.sha256(raw).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "chipseq_runs": len(runs),
        "studies": len({row["study_accession"] for row in runs}),
        "unique_experiments": sum(kind == "experiment" for kind, _ in requests),
        "unique_samples": sum(kind == "sample" for kind, _ in requests),
    }

    atomic_write(args.outdir / "chipseq_runs.tsv", tsv_text(fields, runs))
    atomic_write(
        args.outdir / "ena_request_plan.tsv",
        tsv_text(["record_type", "accession", "source_url"], request_rows)
    )
    atomic_write(
        args.outdir / "request_plan_provenance.json",
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )

    print(json.dumps(provenance, indent=2, sort_keys=True))
    print(f"[OK] Outputs written to: {args.outdir}")


if __name__ == "__main__":
    main()
