#!/usr/bin/env python3
"""Replay only the scientific ENA fields; never fetch or write raw XML."""
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from fetch_chipseq_ena_xml import write_atomic

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "snapshots/chipseq/sanitized_metadata_001"
INVENTORY = ROOT / "snapshots/chipseq/xml_validation_001/xml_inventory.tsv"
FIELDS = {"sample": {"sample_alias", "sample_title", "sample_description"},
          "experiment": {"experiment_title", "library_name"}}
KEYS = {"record_type", "accession", "source_url", "historical_source_sha256",
        "retrieved_at_utc", "fetch_script_sha256", "fields"}


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def catalog():
    proof = json.loads((SOURCE / "provenance.json").read_text())
    data = (SOURCE / "metadata.json").read_bytes()
    if hashlib.sha256(data).hexdigest() != proof["metadata_sha256"]:
        raise ValueError("Sanitized metadata checksum mismatch")
    value = json.loads(data)
    if set(value) != {"schema_version", "records"} or value["schema_version"] != 1:
        raise ValueError("Unsupported sanitized snapshot")
    records = {}
    for row in value["records"]:
        if set(row) != KEYS or row["record_type"] not in FIELDS:
            raise ValueError("Unexpected metadata/contact fields")
        if set(row["fields"]) != FIELDS[row["record_type"]]:
            raise ValueError("Unexpected scientific fields")
        if any(not isinstance(v, str) or re.search(r"[\w.%+-]+@[\w.-]+\.[A-Za-z]{2,}", v)
               for v in row["fields"].values()):
            raise ValueError("Invalid scientific field/contact identifier")
        if row["accession"] in records:
            raise ValueError("Duplicate metadata accession")
        records[row["accession"]] = row
    if len(records) != proof["record_count"]:
        raise ValueError("Sanitized record count mismatch")
    return records


def materialize(plan, outdir, accession):
    if Path(plan).read_bytes() != (SOURCE / "ena_request_plan.tsv").read_bytes():
        raise ValueError("Active and frozen ENA query plans differ")
    record = catalog()[accession]
    with Path(plan).open() as handle:
        requested = next(r for r in csv.DictReader(handle, delimiter="\t") if r["accession"] == accession)
    if any(record[k] != requested[k] for k in ("accession", "record_type", "source_url")):
        raise ValueError("Scientific metadata query mismatch")
    root = Path(outdir); dest = root / "records" / accession
    dest.mkdir(parents=True, exist_ok=True)
    raw = encoded(record)
    receipt = {k: record[k] for k in ("accession", "record_type", "source_url", "retrieved_at_utc", "fetch_script_sha256")}
    receipt.update(sha256=record["historical_source_sha256"], structured_sha256=hashlib.sha256(raw).hexdigest(),
                   scope="sanitized_historical_scientific_metadata")
    for path, data in ((root / "request_plan.tsv", Path(plan).read_bytes()),
                       (dest / "record.json", raw), (dest / "receipt.json", encoded(receipt))):
        if path.exists() and path.read_bytes() != data:
            raise ValueError("Sanitized metadata destination collision")
        write_atomic(path, data)


def verified_record(snapshot, item, accession, kind):
    expected = catalog()[accession]
    raw = (Path(snapshot) / "records" / accession / "record.json").read_bytes()
    receipt = json.loads((Path(snapshot) / "records" / accession / "receipt.json").read_text())
    if (json.loads(raw) != expected or hashlib.sha256(raw).hexdigest() != receipt["structured_sha256"]
            or item["status"] != "verified" or item["record_type"] != kind
            or expected["record_type"] != kind or item["sha256"] != expected["historical_source_sha256"]
            or receipt["sha256"] != expected["historical_source_sha256"]):
        raise ValueError("Sanitized scientific metadata/inventory mismatch")
    return expected


def record_element(snapshot, inventory, accession, kind):
    # Keep reviewed label extraction unchanged by recreating its minimal input
    # interface in memory. No XML is fetched, written or published.
    row = verified_record(snapshot, inventory[accession], accession, kind)
    fields = row["fields"]; element = ET.Element(kind.upper(), accession=accession)
    if kind == "sample":
        element.set("alias", fields["sample_alias"])
        ET.SubElement(element, "TITLE").text = fields["sample_title"]
        ET.SubElement(element, "DESCRIPTION").text = fields["sample_description"]
    else:
        ET.SubElement(element, "TITLE").text = fields["experiment_title"]
        library = ET.SubElement(ET.SubElement(element, "DESIGN"), "LIBRARY_DESCRIPTOR")
        ET.SubElement(library, "LIBRARY_NAME").text = fields["library_name"]
    return element


def summarize(plan, snapshot, output):
    raw = Path(plan).read_bytes()
    if raw != (SOURCE / "ena_request_plan.tsv").read_bytes() or (Path(snapshot) / "request_plan.tsv").read_bytes() != raw:
        raise ValueError("Metadata query differs from frozen definition")
    with INVENTORY.open() as handle:
        inventory = list(csv.DictReader(handle, delimiter="\t"))
    if {r["accession"] for r in inventory} != set(catalog()):
        raise ValueError("Historical inventory/catalog mismatch")
    for row in inventory:
        verified_record(snapshot, row, row["accession"], row["record_type"])
    # Inventory hashes describe historical XML, not a published XML body.
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    write_atomic(Path(output), INVENTORY.read_bytes())
