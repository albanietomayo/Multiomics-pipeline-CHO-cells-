#!/usr/bin/env python3

import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

if len(sys.argv) != 5:
    raise SystemExit(
        "Usage: resolve_gaidukov_liftover.py "
        "loci.tsv metadata.tsv alignments.paf output_dir"
    )

loci_path = Path(sys.argv[1])
metadata_path = Path(sys.argv[2])
paf_path = Path(sys.argv[3])
outdir = Path(sys.argv[4])

outdir.mkdir(
    parents=True,
    exist_ok=True,
)

MIN_QUERY_COVERAGE = 0.95
MIN_IDENTITY = 0.99
MIN_MAPQ = 30

TARGET_ASSEMBLY = "GCF_003668045.3"

# ------------------------------------------------------------
# Metadata
# ------------------------------------------------------------

with metadata_path.open(newline="") as fh:
    metadata = {
        row["query_id"]: row
        for row in csv.DictReader(
            fh,
            delimiter="\t",
        )
    }

if len(metadata) != 21:
    raise SystemExit(
        f"ERROR: expected 21 query metadata rows, "
        f"found {len(metadata)}"
    )

# ------------------------------------------------------------
# PAF
# ------------------------------------------------------------

alignments = defaultdict(list)

with paf_path.open() as fh:

    for line in fh:

        if not line.strip():
            continue

        f = line.rstrip().split("\t")

        if len(f) < 12:
            raise SystemExit(
                "ERROR: malformed PAF"
            )

        tags = {}

        for field in f[12:]:

            p = field.split(
                ":",
                2,
            )

            if len(p) == 3:
                tags[p[0]] = p[2]

        qlen = int(f[1])
        qstart = int(f[2])
        qend = int(f[3])

        nmatch = int(f[9])
        blocklen = int(f[10])

        aln = {
            "qname": f[0],
            "qlen": qlen,
            "qstart": qstart,
            "qend": qend,
            "strand": f[4],
            "tname": f[5],
            "tlen": int(f[6]),
            "tstart": int(f[7]),
            "tend": int(f[8]),
            "nmatch": nmatch,
            "blocklen": blocklen,
            "mapq": int(f[11]),
            "qcov": (
                (qend - qstart) / qlen
                if qlen
                else 0
            ),
            "identity": (
                nmatch / blocklen
                if blocklen
                else 0
            ),
            "tp": tags.get("tp", ""),
            "cigar": tags.get("cg", ""),
        }

        alignments[
            aln["qname"]
        ].append(aln)

# ------------------------------------------------------------
# Exact source-point projection through minimap2 CIGAR
# ------------------------------------------------------------

def project_query_base(aln, query_pos0):

    cigar = aln["cigar"]

    if not cigar:
        return None

    ops = re.findall(
        r"(\d+)([MIDNSHP=X])",
        cigar,
    )

    if not ops:
        return None

    tcur = aln["tstart"]

    if aln["strand"] == "+":

        qcur = aln["qstart"]

        for count_text, op in ops:

            n = int(count_text)

            if op in {"M", "=", "X"}:

                if (
                    qcur
                    <= query_pos0
                    < qcur + n
                ):
                    return (
                        tcur
                        + query_pos0
                        - qcur
                    )

                qcur += n
                tcur += n

            elif op in {"I", "S"}:

                if (
                    qcur
                    <= query_pos0
                    < qcur + n
                ):
                    return None

                qcur += n

            elif op in {"D", "N"}:
                tcur += n

    else:

        qcur = aln["qend"] - 1

        for count_text, op in ops:

            n = int(count_text)

            if op in {"M", "=", "X"}:

                low = qcur - n + 1

                if (
                    low
                    <= query_pos0
                    <= qcur
                ):
                    return (
                        tcur
                        + qcur
                        - query_pos0
                    )

                qcur -= n
                tcur += n

            elif op in {"I", "S"}:

                low = qcur - n + 1

                if (
                    low
                    <= query_pos0
                    <= qcur
                ):
                    return None

                qcur -= n

            elif op in {"D", "N"}:
                tcur += n

    return None

# ------------------------------------------------------------
# Boundary/query QC
# ------------------------------------------------------------

qc_rows = []
resolved = {}

for qname in sorted(metadata):

    meta = metadata[qname]

    offset = int(
        meta["boundary_offset0"]
    )

    all_alignments = alignments.get(
        qname,
        [],
    )

    strong = []

    for aln in all_alignments:

        covers_point = (
            aln["qstart"]
            <= offset
            < aln["qend"]
        )

        if (
            aln["qcov"]
            >= MIN_QUERY_COVERAGE
            and aln["identity"]
            >= MIN_IDENTITY
            and covers_point
        ):
            strong.append(aln)

    status = "unresolved"
    chosen = None
    target_position = None

    if len(strong) > 1:

        status = "mapped_ambiguous"

    elif len(strong) == 1:

        candidate = strong[0]

        if (
            candidate["mapq"]
            >= MIN_MAPQ
            and candidate["tp"]
            in {"", "P"}
        ):

            projected = project_query_base(
                candidate,
                offset,
            )

            if projected is not None:

                chosen = candidate
                target_position = (
                    projected + 1
                )

                status = "mapped_unique"

    qc = {
        **meta,
        "n_total_alignments": len(
            all_alignments
        ),
        "n_strong_alignments": len(
            strong
        ),
        "mapping_status": status,
        "target_seqname": (
            chosen["tname"]
            if chosen
            else ""
        ),
        "target_position": (
            target_position
            if target_position is not None
            else ""
        ),
        "orientation": (
            chosen["strand"]
            if chosen
            else ""
        ),
        "query_coverage": (
            f"{chosen['qcov']:.6f}"
            if chosen
            else ""
        ),
        "identity": (
            f"{chosen['identity']:.6f}"
            if chosen
            else ""
        ),
        "mapq": (
            chosen["mapq"]
            if chosen
            else ""
        ),
        "cigar": (
            chosen["cigar"]
            if chosen
            else ""
        ),
    }

    qc_rows.append(qc)

    resolved[
        meta["source_locus_id"]
    ] = qc

# ------------------------------------------------------------
# Combine with biological/source table
# ------------------------------------------------------------

with loci_path.open(newline="") as fh:
    loci = list(
        csv.DictReader(
            fh,
            delimiter="\t",
        )
    )

if len(loci) != 21:
    raise SystemExit(
        "ERROR: expected 21 Gaidukov loci"
    )

mapped_rows = []

for locus in loci:

    locus_id = locus[
        "source_locus_id"
    ]

    if locus_id not in resolved:
        raise SystemExit(
            f"ERROR: missing mapping result "
            f"for {locus_id}"
        )

    m = resolved[locus_id]

    target_position = (
        m["target_position"]
        if m["mapping_status"]
        == "mapped_unique"
        else ""
    )

    mapped_rows.append(
        {
            **locus,
            "target_assembly_accession": TARGET_ASSEMBLY,
            "target_seqname": (
                m["target_seqname"]
            ),
            "target_start": target_position,
            "target_end": target_position,
            "mapping_orientation": (
                m["orientation"]
            ),
            "query_coverage": (
                m["query_coverage"]
            ),
            "identity": (
                m["identity"]
            ),
            "mapq": (
                m["mapq"]
            ),
            "mapping_status": (
                m["mapping_status"]
            ),
            "mapping_method": (
                "2000bp_flank_minimap2_asm5_"
                "exact_CIGAR_projection"
            ),
        }
    )

# ------------------------------------------------------------
# Write files
# ------------------------------------------------------------

qc_file = (
    outdir
    / "gaidukov_mapping_qc.tsv"
)

mapped_file = (
    outdir
    / "gaidukov_loci_mapped.tsv"
)

with qc_file.open(
    "w",
    newline="",
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            qc_rows[0].keys()
        ),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(qc_rows)


with mapped_file.open(
    "w",
    newline="",
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            mapped_rows[0].keys()
        ),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(mapped_rows)

counts = Counter(
    x["mapping_status"]
    for x in mapped_rows
)

print("GAIDUKOV LIFTOVER RESOLUTION: PASS")
print(f"Source loci:        {len(mapped_rows)}")

for status in sorted(counts):
    print(
        f"{status:22s} "
        f"{counts[status]}"
    )

print()
print(f"Mapping QC: {qc_file}")
print(f"Mapped loci: {mapped_file}")
