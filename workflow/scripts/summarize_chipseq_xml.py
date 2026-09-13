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

    raw = args.plan.read_bytes()
    if (args.snapshot / "request_plan.tsv").read_bytes() != raw:
        raise ValueError("Snapshot and current request plan differ")

    rows = list(csv.DictReader(
        io.StringIO(raw.decode("utf-8-sig")), delimiter="\t"
    ))
    if not rows or len({r["accession"] for r in rows}) != len(rows):
        raise ValueError("Empty plan or duplicate accessions")

    output = []
    for row in sorted(rows, key=lambda r: r["accession"]):
        accession = row["accession"]
        relative = Path("records") / accession / "record.xml"
        data = (args.snapshot / relative).read_bytes()
        receipt = json.loads(
            (args.snapshot / relative.parent / "receipt.json").read_text()
        )

        for key in ("record_type", "accession", "source_url"):
            if receipt[key] != row[key]:
                raise ValueError(f"{accession}: receipt mismatch for {key}")

        if digest(data) != receipt["sha256"]:
            raise ValueError(f"{accession}: SHA-256 mismatch")
        validate_xml(data, row["record_type"], accession)

        output.append(dict(
            row,
            xml_path=relative.as_posix(),
            sha256=receipt["sha256"],
            retrieved_at_utc=receipt["retrieved_at_utc"],
            fetch_script_sha256=receipt["fetch_script_sha256"],
            plan_sha256=digest(raw),
            status="verified",
        ))

    fields = [
        "record_type", "accession", "source_url", "xml_path", "sha256",
        "retrieved_at_utc", "fetch_script_sha256", "plan_sha256", "status",
    ]
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_atomic(args.output, buffer.getvalue().encode())
    print(f"[OK] Verified XML records: {len(output)}")


if __name__ == "__main__":
    main()
