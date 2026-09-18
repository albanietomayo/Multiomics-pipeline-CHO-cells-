#!/usr/bin/env python3

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

if len(sys.argv) != 5:
    raise SystemExit(
        "Usage: resolve_dhiman_liftover.py "
        "loci.tsv metadata.tsv alignments.paf output_dir"
    )

loci_path = Path(sys.argv[1])
metadata_path = Path(sys.argv[2])
paf_path = Path(sys.argv[3])
outdir = Path(sys.argv[4])

outdir.mkdir(parents=True, exist_ok=True)

MIN_QUERY_COVERAGE = 0.95
MIN_IDENTITY = 0.99
MIN_MAPQ = 30

# ------------------------------------------------------------
# Read query metadata
# ------------------------------------------------------------

with metadata_path.open(newline="") as fh:
    metadata = {
        row["query_id"]: row
        for row in csv.DictReader(fh, delimiter="\t")
    }

# ------------------------------------------------------------
# Parse minimap2 PAF
# ------------------------------------------------------------

alignments = defaultdict(list)

with paf_path.open() as fh:
    for line in fh:
        if not line.strip():
            continue

        f = line.rstrip().split("\t")

        if len(f) < 12:
            raise SystemExit("ERROR: malformed PAF line")

        tags = {}

        for field in f[12:]:
            pieces = field.split(":", 2)

            if len(pieces) == 3:
                tags[pieces[0]] = pieces[2]

        qname = f[0]
        qlen = int(f[1])
        qstart = int(f[2])
        qend = int(f[3])
        strand = f[4]
        tname = f[5]
        tlen = int(f[6])
        tstart = int(f[7])
        tend = int(f[8])
        nmatch = int(f[9])
        blocklen = int(f[10])
        mapq = int(f[11])

        qcov = (
            (qend - qstart) / qlen
            if qlen
            else 0.0
        )

        identity = (
            nmatch / blocklen
            if blocklen
            else 0.0
        )

        alignments[qname].append(
            {
                "qname": qname,
                "qlen": qlen,
                "qstart": qstart,
                "qend": qend,
                "strand": strand,
                "tname": tname,
                "tlen": tlen,
                "tstart": tstart,
                "tend": tend,
                "nmatch": nmatch,
                "blocklen": blocklen,
                "mapq": mapq,
                "qcov": qcov,
                "identity": identity,
                "tp": tags.get("tp", ""),
                "cigar": tags.get("cg", ""),
            }
        )

# ------------------------------------------------------------
# Exact projection of the source boundary through the CIGAR.
# Return 0-based target base coordinate.
# ------------------------------------------------------------

def project_query_base(aln, query_pos0):
    cigar = aln["cigar"]

    if not cigar:
        return None

    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)

    if not ops:
        return None

    tcur = aln["tstart"]

    if aln["strand"] == "+":
        qcur = aln["qstart"]

        for count_text, op in ops:
            n = int(count_text)

            if op in {"M", "=", "X"}:
                if qcur <= query_pos0 < qcur + n:
                    return tcur + (query_pos0 - qcur)

                qcur += n
                tcur += n

            elif op in {"I", "S"}:
                if qcur <= query_pos0 < qcur + n:
                    return None
                qcur += n

            elif op in {"D", "N"}:
                tcur += n

            elif op in {"H", "P"}:
                pass

    else:
        qcur = aln["qend"] - 1

        for count_text, op in ops:
            n = int(count_text)

            if op in {"M", "=", "X"}:
                low = qcur - n + 1
                high = qcur

                if low <= query_pos0 <= high:
                    return tcur + (qcur - query_pos0)

                qcur -= n
                tcur += n

            elif op in {"I", "S"}:
                low = qcur - n + 1
                high = qcur

                if low <= query_pos0 <= high:
                    return None

                qcur -= n

            elif op in {"D", "N"}:
                tcur += n

            elif op in {"H", "P"}:
                pass

    return None

# ------------------------------------------------------------
# Resolve each boundary fail-closed.
# ------------------------------------------------------------

boundary_results = []

for qname, meta in sorted(metadata.items()):
    offset = int(meta["boundary_offset0"])

    candidates = []

    for aln in alignments.get(qname, []):
        boundary_covered = (
            aln["qstart"]
            <= offset
            < aln["qend"]
        )

        strong = (
            aln["qcov"] >= MIN_QUERY_COVERAGE
            and aln["identity"] >= MIN_IDENTITY
            and boundary_covered
        )

        if strong:
            candidates.append(aln)

    status = "unresolved"
    chosen = None
    projected = None

    if len(candidates) == 1:
        candidate = candidates[0]

        if (
            candidate["mapq"] >= MIN_MAPQ
            and candidate["tp"] in {"P", ""}
        ):
            projected = project_query_base(
                candidate,
                offset,
            )

            if projected is not None:
                chosen = candidate
                status = "mapped_unique"

    elif len(candidates) > 1:
        status = "mapped_ambiguous"

    row = {
        **meta,
        "n_total_alignments": len(
            alignments.get(qname, [])
        ),
        "n_strong_alignments": len(candidates),
        "mapping_status": status,
        "target_seqname": (
            chosen["tname"] if chosen else ""
        ),
        "target_position": (
            projected + 1
            if projected is not None
            else ""
        ),
        "orientation": (
            chosen["strand"] if chosen else ""
        ),
        "query_coverage": (
            f"{chosen['qcov']:.6f}"
            if chosen else ""
        ),
        "identity": (
            f"{chosen['identity']:.6f}"
            if chosen else ""
        ),
        "mapq": (
            chosen["mapq"] if chosen else ""
        ),
        "cigar": (
            chosen["cigar"] if chosen else ""
        ),
    }

    boundary_results.append(row)

boundary_by_locus = defaultdict(dict)

for row in boundary_results:
    boundary_by_locus[row["source_locus_id"]][
        row["boundary"]
    ] = row

# ------------------------------------------------------------
# Resolve complete loci
# ------------------------------------------------------------

with loci_path.open(newline="") as fh:
    loci = list(csv.DictReader(fh, delimiter="\t"))

mapped_loci = []

for source in loci:
    locus = source["source_locus_id"]
    start = int(source["source_start"])
    end = int(source["source_end"])

    required = ["START"]

    if end != start:
        required.append("END")

    boundaries = boundary_by_locus[locus]

    missing = [
        x for x in required
        if x not in boundaries
    ]

    if missing:
        raise SystemExit(
            f"ERROR: missing boundary result for "
            f"{locus}: {missing}"
        )

    statuses = {
        boundaries[x]["mapping_status"]
        for x in required
    }

    final_status = "unresolved"
    target_seq = ""
    target_start = ""
    target_end = ""
    orientation = ""
    target_span = ""
    source_span = end - start + 1
    span_ratio = ""

    if statuses == {"mapped_unique"}:
        seqs = {
            boundaries[x]["target_seqname"]
            for x in required
        }

        strands = {
            boundaries[x]["orientation"]
            for x in required
        }

        if len(seqs) == 1 and len(strands) == 1:
            positions = [
                int(boundaries[x]["target_position"])
                for x in required
            ]

            target_seq = next(iter(seqs))
            orientation = next(iter(strands))
            target_start = min(positions)
            target_end = max(positions)

            if len(required) == 1:
                target_end = target_start

            target_span = (
                target_end - target_start + 1
            )

            span_ratio = (
                target_span / source_span
                if source_span
                else ""
            )

            # A huge span discrepancy is retained but not
            # automatically accepted as unique.
            if (
                len(required) == 1
                or 0.50 <= span_ratio <= 2.00
            ):
                final_status = "mapped_unique"
            else:
                final_status = "span_inconsistent"

        else:
            final_status = "boundary_inconsistent"

    elif "mapped_ambiguous" in statuses:
        final_status = "mapped_ambiguous"

    mapped_loci.append(
        {
            **source,
            "target_assembly_accession": "GCF_003668045.3",
            "target_seqname": target_seq,
            "target_start": target_start,
            "target_end": target_end,
            "orientation": orientation,
            "source_span_bp": source_span,
            "target_span_bp": target_span,
            "span_ratio": (
                f"{span_ratio:.6f}"
                if isinstance(span_ratio, float)
                else ""
            ),
            "mapping_status": final_status,
            "mapping_method": (
                "2000bp_boundary_flanks_minimap2_asm5"
            ),
        }
    )

# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------

boundary_file = outdir / "dhiman_boundary_mapping_qc.tsv"

boundary_fields = list(boundary_results[0].keys())

with boundary_file.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=boundary_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(boundary_results)

loci_file = outdir / "dhiman_loci_mapped.tsv"

locus_fields = list(mapped_loci[0].keys())

with loci_file.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=locus_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(mapped_loci)

counts = defaultdict(int)

for row in mapped_loci:
    counts[row["mapping_status"]] += 1

print("DHIMAN LIFTOVER RESOLUTION: PASS")
print(f"Boundary queries:     {len(boundary_results)}")
print(f"Source loci:          {len(mapped_loci)}")

for key in sorted(counts):
    print(f"{key:22s} {counts[key]}")

print()
print(f"Boundary QC: {boundary_file}")
print(f"Mapped loci: {loci_file}")
