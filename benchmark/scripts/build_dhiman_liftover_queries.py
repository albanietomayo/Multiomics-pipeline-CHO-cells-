#!/usr/bin/env python3

import csv
import sys
from collections import defaultdict
from pathlib import Path

if len(sys.argv) != 5:
    raise SystemExit(
        "Usage: build_dhiman_liftover_queries.py "
        "loci.tsv source.fa flank output_prefix"
    )

table = Path(sys.argv[1])
fasta = Path(sys.argv[2])
flank = int(sys.argv[3])
prefix = Path(sys.argv[4])

if flank != 2000:
    raise SystemExit("ERROR: expected flank=2000")

# ------------------------------------------------------------
# Build requested source-coordinate windows.
#
# Published RefSeq-style coordinates are interpreted as
# 1-based closed for sequence extraction. The original values
# themselves remain preserved in the source tables.
# ------------------------------------------------------------

queries = []
by_seq = defaultdict(list)

with table.open(newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

if len(rows) != 13:
    raise SystemExit(
        f"ERROR: expected 13 source loci, found {len(rows)}"
    )

for row in rows:
    locus = row["source_locus_id"]
    seq = row["source_seqname"]
    start = int(row["source_start"])
    end = int(row["source_end"])

    boundaries = [("START", start)]

    if end != start:
        boundaries.append(("END", end))

    for boundary, position in boundaries:
        extract_start = max(1, position - flank)
        extract_end = position + flank

        query_id = f"{locus}|{boundary}"

        rec = {
            "query_id": query_id,
            "source_locus_id": locus,
            "boundary": boundary,
            "source_seqname": seq,
            "source_position": position,
            "extract_start": extract_start,
            "extract_end_requested": extract_end,
            "boundary_offset0": position - extract_start,
            "parts": [],
            "observed_seq_end": 0,
        }

        queries.append(rec)
        by_seq[seq].append(rec)

if len(queries) != 21:
    raise SystemExit(
        f"ERROR: expected 21 boundary queries, found {len(queries)}"
    )

# ------------------------------------------------------------
# Stream through the 2.37-Gb FASTA once.
# Only bases overlapping the 21 requested windows are retained.
# ------------------------------------------------------------

wanted = set(by_seq)

current = None
position = 1
seen = set()

with fasta.open() as fh:
    for line in fh:
        if line.startswith(">"):
            current = line[1:].split()[0]
            position = 1

            if current in wanted:
                seen.add(current)

            continue

        sequence = line.strip().upper()

        if not sequence:
            continue

        chunk_start = position
        chunk_end = position + len(sequence) - 1

        if current in by_seq:
            for rec in by_seq[current]:

                left = max(
                    chunk_start,
                    rec["extract_start"],
                )

                right = min(
                    chunk_end,
                    rec["extract_end_requested"],
                )

                if left <= right:
                    a = left - chunk_start
                    b = right - chunk_start + 1

                    rec["parts"].append(sequence[a:b])
                    rec["observed_seq_end"] = max(
                        rec["observed_seq_end"],
                        right,
                    )

        position = chunk_end + 1

missing_seq = sorted(wanted - seen)

if missing_seq:
    raise SystemExit(
        "ERROR: requested FASTA sequences absent: "
        + ",".join(missing_seq)
    )

fasta_out = prefix.with_suffix(".fa")
meta_out = prefix.with_suffix(".metadata.tsv")

fasta_out.parent.mkdir(parents=True, exist_ok=True)

with fasta_out.open("w") as fout, \
     meta_out.open("w", newline="") as mout:

    fields = [
        "query_id",
        "source_locus_id",
        "boundary",
        "source_seqname",
        "source_position",
        "extract_start",
        "extract_end",
        "query_length",
        "boundary_offset0",
    ]

    writer = csv.DictWriter(
        mout,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()

    for rec in queries:
        sequence = "".join(rec["parts"])

        if not sequence:
            raise SystemExit(
                f"ERROR: empty sequence for {rec['query_id']}"
            )

        offset = rec["boundary_offset0"]

        if offset < 0 or offset >= len(sequence):
            raise SystemExit(
                f"ERROR: boundary absent from query "
                f"{rec['query_id']}"
            )

        n_fraction = sequence.count("N") / len(sequence)

        if n_fraction > 0.10:
            raise SystemExit(
                f"ERROR: >10% N in {rec['query_id']}"
            )

        fout.write(f">{rec['query_id']}\n")

        for i in range(0, len(sequence), 80):
            fout.write(sequence[i:i + 80] + "\n")

        writer.writerow(
            {
                "query_id": rec["query_id"],
                "source_locus_id": rec["source_locus_id"],
                "boundary": rec["boundary"],
                "source_seqname": rec["source_seqname"],
                "source_position": rec["source_position"],
                "extract_start": rec["extract_start"],
                "extract_end": rec["observed_seq_end"],
                "query_length": len(sequence),
                "boundary_offset0": offset,
            }
        )

print("QUERY EXTRACTION: PASS")
print(f"Source loci:       {len(rows)}")
print(f"Boundary queries:  {len(queries)}")
print(f"Flank each side:   {flank}")
print(f"FASTA:             {fasta_out}")
print(f"Metadata:          {meta_out}")
