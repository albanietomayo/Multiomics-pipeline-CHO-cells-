#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --job-name=gaidukov_val
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=benchmark/reports/gaidukov2018/slurm_validate_%j.out
#SBATCH --error=benchmark/reports/gaidukov2018/slurm_validate_%j.err

set -euo pipefail

WT="/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_benchmark"

TOOL="$WT/benchmark/raw/tools/minimap2-2.31_x64-linux"

MINIMAP2="$TOOL/minimap2"
K8="$TOOL/k8"
PAFTOOLS="$TOOL/paftools.js"

TARGET="$WT/resources/reference/CriGri-PICRH-1.0/genome.fa"

LIFT="$WT/benchmark/intermediate/gaidukov2018/liftover"
REPORT="$WT/benchmark/reports/gaidukov2018"

QUERIES="$LIFT/gaidukov_queries.fa"
META="$LIFT/gaidukov_queries.metadata.tsv"
CUSTOM="$LIFT/gaidukov_loci_mapped.tsv"
ORIGINAL_PAF="$LIFT/gaidukov_queries.paf"

POINTS="$LIFT/gaidukov_query_points.bed"

PAFTOOLS_PAF="$LIFT/gaidukov_paftools_asm5.paf"
PAFTOOLS_BED="$LIFT/gaidukov_paftools_asm5_liftover.bed"

UNRES_FASTA="$LIFT/gaidukov_unresolved_queries.fa"
UNRES_BED="$LIFT/gaidukov_unresolved_points.bed"

ASM10_PAF="$LIFT/gaidukov_unresolved_asm10.paf"
ASM10_BED="$LIFT/gaidukov_unresolved_asm10_liftover.bed"

cd "$WT"

echo "============================================================"
echo "ENVIRONMENT"
echo "============================================================"

echo "Host:       $(hostname)"
echo "Date UTC:   $(date -u +%FT%TZ)"
echo "Git commit: $(git rev-parse HEAD)"
echo "minimap2:   $("$MINIMAP2" --version)"

echo
echo "============================================================"
echo "A. OFFICIAL PAFTOOLS-COMPATIBLE asm5 ALIGNMENT"
echo "============================================================"

# Important:
# Do NOT use --eqx here because paftools liftover parses M/I/D
# CIGAR operations. This is an independent projection run.
"$MINIMAP2" \
    -t "${SLURM_CPUS_PER_TASK:-4}" \
    -x asm5 \
    -c \
    --cs \
    --secondary=yes \
    -N 20 \
    "$TARGET" \
    "$QUERIES" \
    > "$PAFTOOLS_PAF"

[[ -s "$PAFTOOLS_PAF" ]] || {
    echo "ERROR: empty asm5 validation PAF"
    exit 20
}

echo "PAF records: $(wc -l < "$PAFTOOLS_PAF")"

echo
echo "============================================================"
echo "B. OFFICIAL PAFTOOLS LIFTOVER"
echo "============================================================"

"$K8" "$PAFTOOLS" liftover \
    -q 30 \
    -l 0 \
    -d 1 \
    "$PAFTOOLS_PAF" \
    "$POINTS" \
    > "$PAFTOOLS_BED"

echo "paftools lifted records: $(wc -l < "$PAFTOOLS_BED")"

cat "$PAFTOOLS_BED"

echo
echo "============================================================"
echo "C. COMPARE PAFTOOLS VS CUSTOM STRICT RESOLVER"
echo "============================================================"

python3 - \
    "$CUSTOM" \
    "$PAFTOOLS_BED" \
    "$REPORT/gaidukov_paftools_crosscheck.tsv" <<'PY'
import csv
import sys
from collections import defaultdict
from pathlib import Path

custom_file = Path(sys.argv[1])
paftools_file = Path(sys.argv[2])
report_file = Path(sys.argv[3])

with custom_file.open(newline="") as fh:
    custom = list(
        csv.DictReader(
            fh,
            delimiter="\t",
        )
    )

pt = defaultdict(list)

with paftools_file.open() as fh:

    for line in fh:

        if not line.strip():
            continue

        f = line.rstrip().split("\t")

        if len(f) < 6:
            raise SystemExit(
                "ERROR: malformed paftools BED"
            )

        # paftools constructs:
        # queryid_start_end
        name = f[3]

        try:
            qid, qstart, qend = name.rsplit("_", 2)
        except ValueError:
            raise SystemExit(
                f"ERROR: cannot parse paftools name {name}"
            )

        pt[qid].append(
            {
                "target_seqname": f[0],
                "target_start0": int(f[1]),
                "target_end0": int(f[2]),
                "strand": f[5],
            }
        )

rows = []
errors = []

mapped = [
    x for x in custom
    if x["mapping_status"] == "mapped_unique"
]

if len(mapped) != 17:
    raise SystemExit(
        f"ERROR: expected 17 mapped_unique, "
        f"found {len(mapped)}"
    )

for row in mapped:

    qid = row["source_locus_id"]

    expected = {
        "target_seqname": row["target_seqname"],
        "target_start0": int(row["target_start"]) - 1,
        "target_end0": int(row["target_end"]),
        "strand": row["mapping_orientation"],
    }

    candidates = pt.get(qid, [])

    status = "PASS"

    if len(candidates) != 1:
        status = "FAIL"
        errors.append(
            f"{qid}: paftools candidates={len(candidates)}"
        )

        observed = {
            "target_seqname": "",
            "target_start0": "",
            "target_end0": "",
            "strand": "",
        }

    else:
        observed = candidates[0]

        if observed != expected:
            status = "FAIL"

            errors.append(
                f"{qid}: expected={expected}, "
                f"observed={observed}"
            )

    rows.append(
        {
            "source_locus_id": qid,
            "custom_target": (
                f"{expected['target_seqname']}:"
                f"{expected['target_start0']}-"
                f"{expected['target_end0']}"
            ),
            "custom_strand": expected["strand"],
            "paftools_target": (
                (
                    f"{observed['target_seqname']}:"
                    f"{observed['target_start0']}-"
                    f"{observed['target_end0']}"
                )
                if observed["target_seqname"]
                else ""
            ),
            "paftools_strand": observed["strand"],
            "crosscheck": status,
        }
    )

with report_file.open("w", newline="") as fh:

    fields = list(rows[0].keys())

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(rows)

print("PAFTOOLS CROSS-CHECK")
print(f"Custom mapped_unique: {len(mapped)}")
print(
    "Exact paftools matches:",
    sum(x["crosscheck"] == "PASS" for x in rows),
)

if errors:

    print()
    print("DISCREPANCIES:")

    for error in errors:
        print(" ", error)

    raise SystemExit(
        "ERROR: custom/paftools coordinate discrepancy"
    )

print("INDEPENDENT PROJECTION VALIDATION: PASS")
PY

echo
echo "Cross-check table:"
column -t -s $'\t' \
    "$REPORT/gaidukov_paftools_crosscheck.tsv"

echo
echo "============================================================"
echo "D. DIAGNOSE THE FOUR UNRESOLVED asm5 LOCI"
echo "============================================================"

python3 - \
    "$META" \
    "$CUSTOM" \
    "$ORIGINAL_PAF" \
    "$REPORT/gaidukov_unresolved_asm5_diagnostics.tsv" <<'PY'
import csv
import sys
from collections import defaultdict
from pathlib import Path

meta_file = Path(sys.argv[1])
custom_file = Path(sys.argv[2])
paf_file = Path(sys.argv[3])
report_file = Path(sys.argv[4])

with meta_file.open(newline="") as fh:
    meta = {
        row["query_id"]: row
        for row in csv.DictReader(
            fh,
            delimiter="\t",
        )
    }

with custom_file.open(newline="") as fh:
    custom = list(
        csv.DictReader(
            fh,
            delimiter="\t",
        )
    )

unresolved = {
    row["source_locus_id"]
    for row in custom
    if row["mapping_status"] == "unresolved"
}

alns = defaultdict(list)

with paf_file.open() as fh:

    for line in fh:

        if not line.strip():
            continue

        f = line.rstrip().split("\t")

        qid = f[0]

        if qid not in unresolved:
            continue

        tags = {}

        for item in f[12:]:

            p = item.split(":", 2)

            if len(p) == 3:
                tags[p[0]] = p[2]

        qlen = int(f[1])
        qstart = int(f[2])
        qend = int(f[3])

        block = int(f[10])
        matches = int(f[9])

        offset = int(
            meta[qid]["boundary_offset0"]
        )

        qcov = (
            (qend - qstart) / qlen
            if qlen
            else 0.0
        )

        identity = (
            matches / block
            if block
            else 0.0
        )

        covers = (
            qstart <= offset < qend
        )

        reasons = []

        if qcov < 0.95:
            reasons.append("qcov<0.95")

        if identity < 0.99:
            reasons.append("identity<0.99")

        if not covers:
            reasons.append("site_not_covered")

        if int(f[11]) < 30:
            reasons.append("MAPQ<30")

        if tags.get("tp", "") not in {"P", ""}:
            reasons.append(
                "non_primary:" + tags.get("tp", "")
            )

        alns[qid].append(
            {
                "source_locus_id": qid,
                "target_seqname": f[5],
                "strand": f[4],
                "qstart": qstart,
                "qend": qend,
                "query_length": qlen,
                "query_coverage": f"{qcov:.6f}",
                "identity": f"{identity:.6f}",
                "mapq": f[11],
                "tp": tags.get("tp", ""),
                "site_covered": str(covers).upper(),
                "failure_reasons": (
                    ";".join(reasons)
                    if reasons
                    else "none"
                ),
            }
        )

fields = [
    "source_locus_id",
    "target_seqname",
    "strand",
    "qstart",
    "qend",
    "query_length",
    "query_coverage",
    "identity",
    "mapq",
    "tp",
    "site_covered",
    "failure_reasons",
]

rows = []

for qid in sorted(unresolved):

    if not alns[qid]:

        rows.append(
            {
                "source_locus_id": qid,
                "target_seqname": "",
                "strand": "",
                "qstart": "",
                "qend": "",
                "query_length": meta[qid]["query_length"],
                "query_coverage": "",
                "identity": "",
                "mapq": "",
                "tp": "",
                "site_covered": "",
                "failure_reasons": "no_alignment",
            }
        )

    else:
        rows.extend(alns[qid])

with report_file.open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(rows)

print("UNRESOLVED asm5 DIAGNOSTICS")
print()

for row in rows:

    print(
        "\t".join(
            str(row[x])
            for x in fields
        )
    )
PY

echo
echo "============================================================"
echo "E. EXPLORATORY asm10 FOR THE FOUR UNRESOLVED LOCI"
echo "============================================================"

"$MINIMAP2" \
    -t "${SLURM_CPUS_PER_TASK:-4}" \
    -x asm10 \
    -c \
    --cs \
    --secondary=yes \
    -N 20 \
    "$TARGET" \
    "$UNRES_FASTA" \
    > "$ASM10_PAF"

echo "asm10 PAF records: $(wc -l < "$ASM10_PAF")"

"$K8" "$PAFTOOLS" liftover \
    -q 0 \
    -l 0 \
    -d 1 \
    "$ASM10_PAF" \
    "$UNRES_BED" \
    > "$ASM10_BED"

echo
echo "asm10 paftools candidate projections:"
cat "$ASM10_BED" || true

echo
echo "============================================================"
echo "F. EVALUATE asm10 AGAINST ORIGINAL STRICT CRITERIA"
echo "============================================================"

python3 - \
    "$META" \
    "$CUSTOM" \
    "$ASM10_PAF" \
    "$REPORT/gaidukov_unresolved_asm10_diagnostics.tsv" <<'PY'
import csv
import sys
from collections import defaultdict
from pathlib import Path

meta_file = Path(sys.argv[1])
custom_file = Path(sys.argv[2])
paf_file = Path(sys.argv[3])
report_file = Path(sys.argv[4])

with meta_file.open(newline="") as fh:
    meta = {
        row["query_id"]: row
        for row in csv.DictReader(
            fh,
            delimiter="\t",
        )
    }

with custom_file.open(newline="") as fh:
    custom = list(
        csv.DictReader(
            fh,
            delimiter="\t",
        )
    )

unresolved = {
    row["source_locus_id"]
    for row in custom
    if row["mapping_status"] == "unresolved"
}

alignments = defaultdict(list)

with paf_file.open() as fh:

    for line in fh:

        if not line.strip():
            continue

        f = line.rstrip().split("\t")

        qid = f[0]

        tags = {}

        for item in f[12:]:

            p = item.split(":", 2)

            if len(p) == 3:
                tags[p[0]] = p[2]

        qlen = int(f[1])
        qs = int(f[2])
        qe = int(f[3])

        matches = int(f[9])
        block = int(f[10])
        mapq = int(f[11])

        offset = int(
            meta[qid]["boundary_offset0"]
        )

        qcov = (qe - qs) / qlen
        identity = matches / block if block else 0.0
        covers = qs <= offset < qe

        strong = (
            qcov >= 0.95
            and identity >= 0.99
            and covers
        )

        alignments[qid].append(
            {
                "source_locus_id": qid,
                "target_seqname": f[5],
                "target_start0": f[7],
                "target_end0": f[8],
                "strand": f[4],
                "query_coverage": f"{qcov:.6f}",
                "identity": f"{identity:.6f}",
                "mapq": mapq,
                "tp": tags.get("tp", ""),
                "site_covered": str(covers).upper(),
                "strong_pre_mapq": str(strong).upper(),
            }
        )

fields = [
    "source_locus_id",
    "target_seqname",
    "target_start0",
    "target_end0",
    "strand",
    "query_coverage",
    "identity",
    "mapq",
    "tp",
    "site_covered",
    "strong_pre_mapq",
]

rows = []

print(
    "locus\talignments\tstrong\t"
    "strict_fallback_eligible"
)

for qid in sorted(unresolved):

    a = alignments.get(qid, [])

    strong = [
        x for x in a
        if x["strong_pre_mapq"] == "TRUE"
    ]

    eligible = (
        len(strong) == 1
        and int(strong[0]["mapq"]) >= 30
        and strong[0]["tp"] in {"", "P"}
    )

    print(
        f"{qid}\t"
        f"{len(a)}\t"
        f"{len(strong)}\t"
        f"{str(eligible).upper()}"
    )

    if a:
        rows.extend(a)
    else:
        rows.append(
            {
                "source_locus_id": qid,
                "target_seqname": "",
                "target_start0": "",
                "target_end0": "",
                "strand": "",
                "query_coverage": "",
                "identity": "",
                "mapq": "",
                "tp": "",
                "site_covered": "",
                "strong_pre_mapq": "FALSE",
            }
        )

with report_file.open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(rows)

print()
print("Detailed asm10 alignments:")

for row in rows:
    print(
        "\t".join(
            str(row[x])
            for x in fields
        )
    )
PY

echo
echo "============================================================"
echo "G. CHECKSUMS"
echo "============================================================"

sha256sum \
    "$PAFTOOLS_PAF" \
    "$PAFTOOLS_BED" \
    "$REPORT/gaidukov_paftools_crosscheck.tsv" \
    "$REPORT/gaidukov_unresolved_asm5_diagnostics.tsv" \
    "$ASM10_PAF" \
    "$REPORT/gaidukov_unresolved_asm10_diagnostics.tsv"

echo
echo "GAIDUKOV INDEPENDENT VALIDATION JOB: PASS"
