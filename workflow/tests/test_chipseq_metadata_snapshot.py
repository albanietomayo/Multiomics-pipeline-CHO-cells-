"""Sanitized metadata reproduction must not require or publish raw XML."""
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow/scripts"))
spec = importlib.util.spec_from_file_location("metadata_snapshot_test", ROOT / "workflow/scripts/chipseq_metadata_snapshot.py")
M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)


class SanitizedMetadataTests(unittest.TestCase):
    def test_catalog_contains_only_consumed_fields(self):
        values = M.catalog()
        self.assertEqual(len(values), 232)
        for row in values.values():
            self.assertEqual(set(row), M.KEYS)
            self.assertEqual(set(row["fields"]), M.FIELDS[row["record_type"]])

    def test_materialization_never_writes_xml(self):
        accession = next(iter(M.catalog()))
        with tempfile.TemporaryDirectory() as temp:
            M.materialize(M.SOURCE / "ena_request_plan.tsv", temp, accession)
            self.assertFalse(list(Path(temp).rglob("*.xml")))
            row = M.catalog()[accession]
            item = {"status": "verified", "record_type": row["record_type"], "sha256": row["historical_source_sha256"]}
            self.assertEqual(M.verified_record(temp, item, accession, row["record_type"]), row)

    def test_modified_materialized_label_fails(self):
        accession = next(iter(M.catalog()))
        with tempfile.TemporaryDirectory() as temp:
            M.materialize(M.SOURCE / "ena_request_plan.tsv", temp, accession)
            f = Path(temp) / "records" / accession / "record.json"
            f.write_text("{}")
            row = M.catalog()[accession]
            item = {"status": "verified", "record_type": row["record_type"], "sha256": row["historical_source_sha256"]}
            with self.assertRaisesRegex(ValueError, "mismatch"):
                M.verified_record(temp, item, accession, row["record_type"])

    def test_extra_contact_field_rejected_even_with_updated_digest(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); value = json.loads((M.SOURCE / "metadata.json").read_text())
            value["records"][0]["submitter_name"] = "unneeded contact"
            data = M.encoded(value); (root / "metadata.json").write_bytes(data)
            (root / "provenance.json").write_text(json.dumps({"metadata_sha256": hashlib.sha256(data).hexdigest(), "record_count": 232}))
            with mock.patch.object(M, "SOURCE", root):
                with self.assertRaisesRegex(ValueError, "contact fields"):
                    M.catalog()
