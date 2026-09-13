#!/usr/bin/env python3
"""Retrieve, validate and cache ENA experiment/sample XML records."""

import argparse
import csv
import hashlib
import io
import json
import os
import re
import shutil
import tempfile
import time
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_atomic(path, data):
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(data)
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def validate_xml(data, kind, accession):
    root = ET.fromstring(data)
    tag = kind.upper()
    records = [root] if root.tag == tag else root.findall(tag)
    if len(records) != 1:
        raise ValueError(f"Expected exactly one {tag} record")
    record = records[0]
    identifiers = {record.get("accession", "")}
    identifiers.update(
        (element.text or "").strip()
        for element in record.findall("./IDENTIFIERS/*")
    )
    if accession not in identifiers:
        raise ValueError(f"Response does not identify {accession}")


def retrieve(row, directory):
    accession = row["accession"]
    kind = row["record_type"]
    url = row["source_url"]
    destination = directory / accession

    if destination.exists() and not destination.is_dir():
        raise ValueError(f"Cache destination is not a directory: {destination}")
    if destination.is_dir() and any(destination.iterdir()):
        if not all((destination / name).is_file()
                   for name in ("record.xml", "receipt.json")):
            raise ValueError(f"Incomplete cache: {destination}")
        data = (destination / "record.xml").read_bytes()
        receipt = json.loads((destination / "receipt.json").read_text())
        for key in ("accession", "record_type", "source_url"):
            if receipt[key] != row[key]:
                raise ValueError(f"Cached receipt mismatch: {key}")
        if digest(data) != receipt["sha256"]:
            raise ValueError("Cached XML failed SHA-256 verification")
        validate_xml(data, kind, accession)
        return receipt, "cached"

    for attempt in range(1, 4):
        try:
            time.sleep(0.4)
            request = urllib.request.Request(
                url, headers={"User-Agent": "CHO-TFM-ChIP-metadata/1.0"}
            )
            with urllib.request.urlopen(request, timeout=45) as response:
                data = response.read(10_000_001)
            if len(data) > 10_000_000:
                raise ValueError("XML response exceeds the 10 MB limit")
            validate_xml(data, kind, accession)
            break
        except Exception:
            if attempt == 3:
                raise
            time.sleep(attempt * 2)

    receipt = dict(row)
    receipt.update(
        sha256=digest(data),
        retrieved_at_utc=datetime.now(timezone.utc).isoformat(),
        fetch_script_sha256=digest(Path(__file__).read_bytes()),
    )
    temporary = Path(tempfile.mkdtemp(prefix=".pending-", dir=directory))
    try:
        (temporary / "record.xml").write_bytes(data)
        (temporary / "receipt.json").write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.rename(temporary, destination)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
    return receipt, "downloaded"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--accessions", nargs="+")
    args = parser.parse_args()

    raw = args.plan.read_bytes()
    rows = list(csv.DictReader(
        io.StringIO(raw.decode("utf-8-sig")), delimiter="\t"
    ))
    seen = set()

    for row in rows:
        kind, accession = row["record_type"], row["accession"]
        pattern = {
            "experiment": r"[SED]RX[0-9]+",
            "sample": r"(?:SAM[END][A-Z]*|[SED]RS)[0-9]+",
        }.get(kind)
        if pattern is None or not re.fullmatch(pattern, accession):
            raise ValueError(f"Invalid record: {row}")
        if accession in seen:
            raise ValueError(f"Duplicate accession: {accession}")
        seen.add(accession)
        expected_url = (
            f"https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
        )
        if row["source_url"] != expected_url:
            raise ValueError(f"Unexpected ENA URL for {accession}")

    if not rows:
        raise ValueError("The request plan is empty")

    selected = set(args.accessions) if args.accessions else seen
    if selected - seen:
        raise ValueError(
            f"Accessions absent from plan: {sorted(selected - seen)}"
        )
    records = sorted(
        (r for r in rows if r["accession"] in selected),
        key=lambda r: r["accession"],
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    saved_plan = args.outdir / "request_plan.tsv"
    if saved_plan.exists() and saved_plan.read_bytes() != raw:
        raise ValueError(
            "This snapshot belongs to a different plan; use a new outdir"
        )
    if not saved_plan.exists():
        write_atomic(saved_plan, raw)

    directory = args.outdir / "records"
    directory.mkdir(exist_ok=True)
    report = []

    for number, row in enumerate(records, 1):
        result = dict(
            row, status="error", sha256="", retrieved_at_utc="", error=""
        )
        try:
            receipt, status = retrieve(row, directory)
            result.update(
                status=status,
                sha256=receipt["sha256"],
                retrieved_at_utc=receipt["retrieved_at_utc"],
            )
        except Exception as error:
            result["error"] = str(error)
        report.append(result)
        print(
            f"[{number}/{len(records)}] {row['accession']}: "
            f"{result['status']} {result['error']}",
            flush=True,
        )

    selection_id = digest("\n".join(sorted(selected)).encode())[:16]
    fields = [
        "record_type", "accession", "source_url", "status",
        "sha256", "retrieved_at_utc", "error",
    ]
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(report)
    write_atomic(
        args.outdir / f"fetch_{selection_id}.tsv",
        buffer.getvalue().encode(),
    )
    failures = sum(r["status"] == "error" for r in report)
    print(
        f"RESULT: selected={len(records)} total_plan={len(rows)} "
        f"verified={len(records)-failures} errors={failures}"
    )
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
