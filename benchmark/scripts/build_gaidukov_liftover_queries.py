#!/usr/bin/env python3

import csv
import sys
from collections import defaultdict
from pathlib import Path

if len(sys.argv) != 5:
    raise SystemExit(
        "Usage: build_gaidukov_liftover_queries.py "
        "loci.tsv source.fa flank output_prefix"
    )

table = Path(sys.argv[1])
fasta = Path(sys.argv[2])
flank = int(sys.argv[3])
prefix = Path(sys.argv[4])

if flank != 2000:
    raise SystemExit("ERROR: expected flank=2000")

SEARCH_FLANK = 10000
MAX_N_FRACTION = 0.10
MIN_QUERY_LENGTH = 1000
MIN_SIDE = 250

with table.open(newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

if len(rows) != 21:
    raise SystemExit(
        f"ERROR: expected 21 loci, found {len(rows)}"
    )

# ------------------------------------------------------------
# Extract ±10 kb around every published point in one FASTA pass.
# This gives us room to construct a gap-free local query when
# the standard ±2 kb window intersects an assembly gap.
# ------------------------------------------------------------

requests = []
by_seq = defaultdict(list)

for row in rows:

    start = int(row["source_start"])
    end = int(row["source_end"])

    if start != end:
        raise SystemExit(
            f"ERROR: non-point locus {row['source_locus_id']}"
        )

    pos = start

    rec = {
        "row": row,
        "source_locus_id": row["source_locus_id"],
        "source_seqname": row["source_seqname"],
        "source_position": pos,
        "search_start": max(1, pos - SEARCH_FLANK),
        "search_end_requested": pos + SEARCH_FLANK,
        "parts": [],
        "observed_end": 0,
    }

    requests.append(rec)
    by_seq[row["source_seqname"]].append(rec)

wanted = set(by_seq)
seen = set()

current = None
position = 1

with fasta.open() as fh:

    for line in fh:

        if line.startswith(">"):
            current = line[1:].split()[0]
            position = 1

            if current in wanted:
                seen.add(current)

            continue

        seq = line.strip().upper()

        if not seq:
            continue

        chunk_start = position
        chunk_end = position + len(seq) - 1

        if current in by_seq:

            for rec in by_seq[current]:

                left = max(
                    chunk_start,
                    rec["search_start"],
                )

                right = min(
                    chunk_end,
                    rec["search_end_requested"],
                )

                if left <= right:

                    a = left - chunk_start
                    b = right - chunk_start + 1

                    rec["parts"].append(seq[a:b])
                    rec["observed_end"] = max(
                        rec["observed_end"],
                        right,
                    )

        position = chunk_end + 1

missing = sorted(wanted - seen)

if missing:
    raise SystemExit(
        "ERROR: missing source sequences: "
        + ",".join(missing)
    )

fasta_out = prefix.with_suffix(".fa")
meta_out = prefix.with_suffix(".metadata.tsv")

fields = [
    "query_id",
    "source_locus_id",
    "source_seqname",
    "source_position",
    "extract_start",
    "extract_end",
    "query_length",
    "boundary_offset0",
    "n_fraction",
    "standard_window_n_fraction",
    "query_strategy",
    "left_context_bp",
    "right_context_bp",
]

qc_lines = []

with fasta_out.open("w") as fout, \
     meta_out.open("w", newline="") as mout:

    writer = csv.DictWriter(
        mout,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()

    for rec in requests:

        search_seq = "".join(rec["parts"])

        if not search_seq:
            raise SystemExit(
                f"ERROR: empty sequence for "
                f"{rec['source_locus_id']}"
            )

        pos = rec["source_position"]

        # Published coordinate inside extracted ±10 kb region.
        offset_search = pos - rec["search_start"]

        if not (0 <= offset_search < len(search_seq)):
            raise SystemExit(
                f"ERROR: published coordinate absent from "
                f"{rec['source_locus_id']}"
            )

        site_base = search_seq[offset_search]

        if site_base not in {"A", "C", "G", "T"}:
            raise SystemExit(
                f"ERROR: published coordinate for "
                f"{rec['source_locus_id']} is '{site_base}' "
                "in CriGri_1.0; exact projection impossible"
            )

        # ----------------------------------------------------
        # Standard ±2 kb query.
        # ----------------------------------------------------

        standard_start = max(
            1,
            pos - flank,
        )

        standard_end = min(
            rec["observed_end"],
            pos + flank,
        )

        a = standard_start - rec["search_start"]
        b = standard_end - rec["search_start"] + 1

        standard_seq = search_seq[a:b]

        standard_n = (
            standard_seq.count("N")
            / len(standard_seq)
        )

        if standard_n <= MAX_N_FRACTION:

            query_seq = standard_seq
            extract_start = standard_start
            extract_end = standard_end
            strategy = "standard_2000bp_flank"

        else:

            # ------------------------------------------------
            # Find maximal continuous A/C/G/T block containing
            # the exact published coordinate.
            # ------------------------------------------------

            left = offset_search

            while (
                left > 0
                and search_seq[left - 1]
                in {"A", "C", "G", "T"}
            ):
                left -= 1

            right = offset_search

            while (
                right + 1 < len(search_seq)
                and search_seq[right + 1]
                in {"A", "C", "G", "T"}
            ):
                right += 1

            # Do not use >2 kb on either side even if the
            # contiguous block is much larger.
            desired_left = max(
                left,
                offset_search - flank,
            )

            desired_right = min(
                right,
                offset_search + flank,
            )

            query_seq = search_seq[
                desired_left:desired_right + 1
            ]

            extract_start = (
                rec["search_start"]
                + desired_left
            )

            extract_end = (
                rec["search_start"]
                + desired_right
            )

            strategy = (
                "adaptive_gap_free_block_"
                "containing_published_site"
            )

        query_offset = (
            pos - extract_start
        )

        left_context = query_offset
        right_context = (
            len(query_seq)
            - query_offset
            - 1
        )

        n_fraction = (
            query_seq.count("N")
            / len(query_seq)
        )

        # ----------------------------------------------------
        # Fail-closed QC
        # ----------------------------------------------------

        if query_seq[query_offset] not in {
            "A", "C", "G", "T"
        }:
            raise SystemExit(
                f"ERROR: exact site not resolved for "
                f"{rec['source_locus_id']}"
            )

        if len(query_seq) < MIN_QUERY_LENGTH:
            raise SystemExit(
                f"ERROR: usable gap-free query too short for "
                f"{rec['source_locus_id']}: "
                f"{len(query_seq)} bp"
            )

        if left_context < MIN_SIDE:
            raise SystemExit(
                f"ERROR: insufficient left context for "
                f"{rec['source_locus_id']}: "
                f"{left_context} bp"
            )

        if right_context < MIN_SIDE:
            raise SystemExit(
                f"ERROR: insufficient right context for "
                f"{rec['source_locus_id']}: "
                f"{right_context} bp"
            )

        if n_fraction > MAX_N_FRACTION:
            raise SystemExit(
                f"ERROR: final query still has too many N "
                f"for {rec['source_locus_id']}: "
                f"{n_fraction:.4f}"
            )

        fout.write(
            f">{rec['source_locus_id']}\n"
        )

        for i in range(0, len(query_seq), 80):
            fout.write(
                query_seq[i:i + 80] + "\n"
            )

        writer.writerow(
            {
                "query_id": rec["source_locus_id"],
                "source_locus_id": rec["source_locus_id"],
                "source_seqname": rec["source_seqname"],
                "source_position": pos,
                "extract_start": extract_start,
                "extract_end": extract_end,
                "query_length": len(query_seq),
                "boundary_offset0": query_offset,
                "n_fraction": f"{n_fraction:.6f}",
                "standard_window_n_fraction": (
                    f"{standard_n:.6f}"
                ),
                "query_strategy": strategy,
                "left_context_bp": left_context,
                "right_context_bp": right_context,
            }
        )

        qc_lines.append(
            (
                rec["source_locus_id"],
                standard_n,
                n_fraction,
                len(query_seq),
                left_context,
                right_context,
                strategy,
            )
        )

print("GAIDUKOV QUERY EXTRACTION: PASS")
print(f"Source loci:       {len(rows)}")
print(f"Sequence queries:  {len(requests)}")
print()

print(
    "locus\tstandard_N\tfinal_N\tquery_bp\t"
    "left_bp\tright_bp\tstrategy"
)

for x in qc_lines:
    print(
        f"{x[0]}\t"
        f"{x[1]:.4f}\t"
        f"{x[2]:.4f}\t"
        f"{x[3]}\t"
        f"{x[4]}\t"
        f"{x[5]}\t"
        f"{x[6]}"
    )

print()
print(f"FASTA:    {fasta_out}")
print(f"Metadata: {meta_out}")
