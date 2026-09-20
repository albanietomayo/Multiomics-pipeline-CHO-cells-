#!/usr/bin/env python3

import csv
import sys
from collections import Counter

guides_file, source_fa, canon_fa, output = sys.argv[1:5]


def rc(seq):
    return seq.translate(
        str.maketrans("ACGT", "TGCA")
    )[::-1]


with open(guides_file, encoding="utf-8") as fh:
    guides = list(csv.DictReader(fh, delimiter="\t"))


for g in guides:
    seq = g["reported_sequence_5to3"]

    if len(seq) != 23:
        raise SystemExit(
            f"ERROR: {g['guide_id']} length={len(seq)}"
        )

    if g["pam_side"] == "3prime":
        if seq[-2:] != "GG":
            raise SystemExit(
                f"ERROR: {g['guide_id']} invalid NGG"
            )

    elif g["pam_side"] == "5prime":
        if seq[:2] != "CC":
            raise SystemExit(
                f"ERROR: {g['guide_id']} invalid CCN"
            )

    else:
        raise SystemExit("ERROR: invalid PAM side")


def scan(path, assembly):

    hits = []

    name = None
    chunks = []

    def process(name, chunks):

        if name is None:
            return

        genome = "".join(chunks).upper()

        for g in guides:

            reported = g["reported_sequence_5to3"]

            for orientation, motif in (
                ("+", reported),
                ("-", rc(reported)),
            ):

                pos = 0

                while True:

                    i = genome.find(motif, pos)

                    if i < 0:
                        break

                    full_start = i + 1
                    full_end = i + 23

                    pam_side = g["pam_side"]

                    if orientation == "-":
                        pam_side = (
                            "5prime"
                            if pam_side == "3prime"
                            else "3prime"
                        )

                    if pam_side == "3prime":

                        prot_start = full_start
                        prot_end = full_end - 3

                        pam_start = full_end - 2
                        pam_end = full_end

                        cut_left = prot_end - 3
                        cut_right = prot_end - 2

                    else:

                        pam_start = full_start
                        pam_end = full_start + 2

                        prot_start = full_start + 3
                        prot_end = full_end

                        cut_left = prot_start + 2
                        cut_right = prot_start + 3

                    pam = genome[
                        pam_start - 1:pam_end
                    ]

                    pam_valid = (
                        pam[1:] == "GG"
                        if pam_side == "3prime"
                        else pam[:2] == "CC"
                    )

                    hits.append({
                        "assembly": assembly,
                        "site": g["site"],
                        "guide_id": g["guide_id"],
                        "reported_sequence_5to3":
                            reported,
                        "match_orientation":
                            orientation,
                        "seqname": name,
                        "full_start_1based":
                            full_start,
                        "full_end_1based":
                            full_end,
                        "genomic_pam_side":
                            pam_side,
                        "protospacer_start_1based":
                            prot_start,
                        "protospacer_end_1based":
                            prot_end,
                        "pam_start_1based":
                            pam_start,
                        "pam_end_1based":
                            pam_end,
                        "pam_sequence_reference":
                            pam,
                        "pam_valid":
                            str(pam_valid).upper(),
                        "cut_boundary_1based":
                            f"{cut_left}|{cut_right}",
                    })

                    pos = i + 1

    with open(path, encoding="ascii") as fh:

        for line in fh:

            if line.startswith(">"):

                process(name, chunks)

                name = line[1:].split()[0]
                chunks = []

            else:
                chunks.append(line.strip())

        process(name, chunks)

    return hits


rows = []

rows += scan(
    source_fa,
    "GCA_003668045.1",
)

rows += scan(
    canon_fa,
    "GCF_003668045.3",
)

fields = [
    "assembly",
    "site",
    "guide_id",
    "reported_sequence_5to3",
    "match_orientation",
    "seqname",
    "full_start_1based",
    "full_end_1based",
    "genomic_pam_side",
    "protospacer_start_1based",
    "protospacer_end_1based",
    "pam_start_1based",
    "pam_end_1based",
    "pam_sequence_reference",
    "pam_valid",
    "cut_boundary_1based",
]

with open(
    output,
    "w",
    newline="",
    encoding="utf-8",
) as out:

    w = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()
    w.writerows(rows)


counts = Counter(
    (x["assembly"], x["guide_id"])
    for x in rows
)

for assembly in (
    "GCA_003668045.1",
    "GCF_003668045.3",
):

    for g in guides:

        n = counts[
            assembly,
            g["guide_id"],
        ]

        print(
            assembly,
            g["guide_id"],
            f"exact_hits={n}",
            sep="\t",
        )
