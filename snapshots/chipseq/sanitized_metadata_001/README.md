# Minimal historical scientific metadata

This snapshot retains only the fields consumed by ChIP-seq label/condition
planning: sample alias, title and description; experiment title and library
name. Each accession retains its query URL, retrieval date and historical
source checksum. Submitter identities, contacts and unrelated XML fields are
excluded. The source inventory and validation summary remain in
`../xml_validation_001/`; its XML paths describe historical records, not
published XML bodies.

The portable workflow materializes JSON, verifies its structured checksum and
reconstructs only the reviewed planning functions' input interface in memory.
Biological label extraction, eligibility and pairings are unchanged. Both
generated planning tables must match their frozen full-record equivalents.

`provenance.json` records the historical archive checksum and exact query
definition. Raw XML and `ena_xml_snapshot.tar.gz` must never be published in
Git, GitHub release assets or Zenodo. The historical archive remains private;
it is not required to execute metadata replay.
