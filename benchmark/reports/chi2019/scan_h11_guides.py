#!/usr/bin/env python3
import csv
import sys

queries = {
    "cH11g1": "GTATACACTTGAGCCAGTAGTGG",
    "cH11g2": "TATACACTTGAGCCAGTAGTGGG",
}

def rc(s):
    return s.translate(
        str.maketrans("ACGT", "TGCA")
    )[::-1]

def scan_fasta(path, assembly):
    hits = []

    name = None
    seq = []

    def process(name, seq):
        if name is None:
            return

        s = "".join(seq).upper()

        for guide_name, guide23 in queries.items():

            for strand, motif in (
                ("+", guide23),
                ("-", rc(guide23)),
            ):
                start = 0

                while True:
                    i = s.find(motif, start)
                    if i < 0:
                        break

                    full_start = i + 1
                    full_end = i + 23

                    if strand == "+":
                        ps = full_start
                        pe = full_start + 19
                        pam_s = full_start + 20
                        pam_e = full_start + 22
                        cut_left = pe - 3
                        cut_right = pe - 2
                    else:
                        pam_s = full_start
                        pam_e = full_start + 2
                        ps = full_start + 3
                        pe = full_start + 22
                        cut_left = ps + 2
                        cut_right = ps + 3

                    hits.append({
                        "assembly": assembly,
                        "guide": guide_name,
                        "seqname": name,
                        "full_start_1based": full_start,
                        "full_end_1based": full_end,
                        "strand": strand,
                        "protospacer_start_1based": ps,
                        "protospacer_end_1based": pe,
                        "pam_start_1based": pam_s,
                        "pam_end_1based": pam_e,
                        "pam_guide_orientation": guide23[20:23],
                        "cut_boundary_1based":
                            f"{cut_left}|{cut_right}",
                    })

                    start = i + 1

    with open(path, encoding="ascii") as f:
        for line in f:
            if line.startswith(">"):
                process(name, seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())

        process(name, seq)

    return hits


canon, old, output = sys.argv[1:4]

hits = []
hits += scan_fasta(canon, "GCF_003668045.3")
hits += scan_fasta(old, "GCF_000223135.1")

fields = [
    "assembly",
    "guide",
    "seqname",
    "full_start_1based",
    "full_end_1based",
    "strand",
    "protospacer_start_1based",
    "protospacer_end_1based",
    "pam_start_1based",
    "pam_end_1based",
    "pam_guide_orientation",
    "cut_boundary_1based",
]

with open(output, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    w.writeheader()
    w.writerows(hits)

for assembly in (
    "GCF_003668045.3",
    "GCF_000223135.1",
):
    for guide in ("cH11g1", "cH11g2"):
        n = sum(
            h["assembly"] == assembly
            and h["guide"] == guide
            for h in hits
        )
        print(f"{assembly}\t{guide}\texact_full23_hits={n}")
