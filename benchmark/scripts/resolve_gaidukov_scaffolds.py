#!/usr/bin/env python3

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

if len(sys.argv) != 4:
    raise SystemExit(
        "Usage: resolve_gaidukov_scaffolds.py "
        "gaidukov_table1.tsv assembly_report.txt output.tsv"
    )

table_path = Path(sys.argv[1])
report_path = Path(sys.argv[2])
output_path = Path(sys.argv[3])

# ------------------------------------------------------------
# Parse NCBI assembly report
# ------------------------------------------------------------

records = []

with report_path.open() as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue

        f = line.rstrip("\n").split("\t")

        if len(f) < 9:
            raise SystemExit(
                "ERROR: malformed NCBI assembly-report row"
            )

        records.append(
            {
                "sequence_name": f[0],
                "sequence_role": f[1],
                "genbank_accession": f[4],
                "relationship": f[5],
                "refseq_accession": f[6],
                "assembly_unit": f[7],
                "sequence_length": int(f[8]),
            }
        )

if not records:
    raise SystemExit(
        "ERROR: no sequences parsed from assembly report"
    )

# CriGri_1.0 uses names such as scaffold329,
# scaffold934, scaffold1924, etc.
pattern = re.compile(r"^scaffold(\d+)$", re.I)

by_number = defaultdict(list)

for rec in records:
    m = pattern.match(rec["sequence_name"])

    if m:
        by_number[int(m.group(1))].append(rec)

print(
    f"Assembly sequences parsed: {len(records)}"
)
print(
    f"Numbered scaffolds parsed: {len(by_number)}"
)

# ------------------------------------------------------------
# Load Gaidukov Table 1
# ------------------------------------------------------------

with table_path.open(newline="") as fh:
    loci = list(
        csv.DictReader(
            fh,
            delimiter="\t",
        )
    )

if len(loci) != 21:
    raise SystemExit(
        f"ERROR: expected 21 Gaidukov loci, found {len(loci)}"
    )

resolved = []

for locus in loci:

    number = int(
        locus["source_scaffold_number"]
    )

    candidates = by_number.get(number, [])

    if len(candidates) != 1:
        raise SystemExit(
            f"ERROR: scaffold{number} resolves to "
            f"{len(candidates)} assembly-report records"
        )

    rec = candidates[0]

    position = int(
        locus["source_position"]
    )

    if not (
        1 <= position <= rec["sequence_length"]
    ):
        raise SystemExit(
            f"ERROR: {locus['source_locus_id']} "
            f"coordinate {position} outside "
            f"scaffold{number} length "
            f"{rec['sequence_length']}"
        )

    refseq = rec["refseq_accession"]

    if not refseq.startswith("NW_"):
        raise SystemExit(
            f"ERROR: unexpected RefSeq accession "
            f"for scaffold{number}: {refseq}"
        )

    resolved.append(
        {
            **locus,
            "source_assembly_name": "CriGri_1.0",
            "source_assembly_accession": "GCF_000223135.1",
            "source_sequence_name": rec["sequence_name"],
            "source_seqname": refseq,
            "source_genbank_accession": rec[
                "genbank_accession"
            ],
            "source_start": position,
            "source_end": position,
            "source_sequence_length": rec[
                "sequence_length"
            ],
        }
    )

# ------------------------------------------------------------
# Final QC
# ------------------------------------------------------------

if len(resolved) != 21:
    raise SystemExit(
        "ERROR: did not resolve all 21 loci"
    )

if len({
    x["source_locus_id"]
    for x in resolved
}) != 21:
    raise SystemExit(
        "ERROR: duplicated source_locus_id"
    )

fields = list(resolved[0].keys())

output_path.parent.mkdir(
    parents=True,
    exist_ok=True,
)

with output_path.open(
    "w",
    newline="",
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(resolved)

print()
print("GAIDUKOV SCAFFOLD RESOLUTION: PASS")
print(f"Physical loci:          {len(resolved)}")
print(
    "Unique source scaffolds:",
    len({
        x["source_seqname"]
        for x in resolved
    }),
)

print()
print(
    "paper_scaffold\tRefSeq_accession\t"
    "position\tlocus"
)

for row in resolved:
    print(
        f"scaf{row['source_scaffold_number']}\t"
        f"{row['source_seqname']}\t"
        f"{row['source_position']}\t"
        f"{row['source_locus_id']}"
    )
