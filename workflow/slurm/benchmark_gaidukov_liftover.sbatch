#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --job-name=gaidukov_lift
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=benchmark/reports/gaidukov2018/slurm_%j.out
#SBATCH --error=benchmark/reports/gaidukov2018/slurm_%j.err

set -euo pipefail

WT="/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_benchmark"

MINIMAP2="$WT/benchmark/raw/tools/minimap2-2.31_x64-linux/minimap2"

SOURCE_FASTA="$WT/benchmark/raw/gaidukov2018/reference/GCF_000223135.1/genome.fa"
TARGET_FASTA="$WT/resources/reference/CriGri-PICRH-1.0/genome.fa"

LOCI="$WT/benchmark/intermediate/gaidukov2018/gaidukov_source_loci.tsv"

LIFT="$WT/benchmark/intermediate/gaidukov2018/liftover"
PREFIX="$LIFT/gaidukov_queries"

cd "$WT"

mkdir -p "$LIFT"

echo "============================================================"
echo "ENVIRONMENT"
echo "============================================================"

echo "Host:       $(hostname)"
echo "Date UTC:   $(date -u +%FT%TZ)"
echo "Git commit: $(git rev-parse HEAD)"
echo "minimap2:   $MINIMAP2"
echo "version:    $("$MINIMAP2" --version)"

echo
echo "============================================================"
echo "BUILD GAIDUKOV ±2 kb QUERIES"
echo "============================================================"

python3 \
    benchmark/scripts/build_gaidukov_liftover_queries.py \
    "$LOCI" \
    "$SOURCE_FASTA" \
    2000 \
    "$PREFIX"

echo
echo "============================================================"
echo "ALIGN TO GCF_003668045.3"
echo "============================================================"

"$MINIMAP2" \
    -t "${SLURM_CPUS_PER_TASK:-4}" \
    -x asm5 \
    -c \
    --eqx \
    --secondary=yes \
    -N 20 \
    "$TARGET_FASTA" \
    "${PREFIX}.fa" \
    > "${PREFIX}.paf"

[[ -s "${PREFIX}.paf" ]] || {
    echo "ERROR: minimap2 returned empty PAF"
    exit 21
}

echo "PAF records: $(wc -l < "${PREFIX}.paf")"

echo
echo "============================================================"
echo "STRICT RESOLUTION"
echo "============================================================"

python3 \
    benchmark/scripts/resolve_gaidukov_liftover.py \
    "$LOCI" \
    "${PREFIX}.metadata.tsv" \
    "${PREFIX}.paf" \
    "$LIFT"

echo
echo "============================================================"
echo "MAPPING SUMMARY"
echo "============================================================"

python3 - "$LIFT/gaidukov_loci_mapped.tsv" <<'PY'
import csv
import sys
from collections import Counter

with open(sys.argv[1], newline="") as fh:
    rows = list(
        csv.DictReader(
            fh,
            delimiter="\t",
        )
    )

counts = Counter(
    x["mapping_status"]
    for x in rows
)

print("Total loci:", len(rows))

for status in sorted(counts):
    print(
        f"{status:22s}",
        counts[status]
    )

print()
print(
    "locus\tclone\ttarget\tposition\t"
    "orientation\tidentity\tMAPQ\tstatus"
)

for row in rows:
    print(
        "\t".join(
            [
                row["source_locus_id"],
                row["lp_clone"],
                row["target_seqname"],
                row["target_start"],
                row["mapping_orientation"],
                row["identity"],
                row["mapq"],
                row["mapping_status"],
            ]
        )
    )
PY

echo
echo "============================================================"
echo "CHECKSUMS"
echo "============================================================"

sha256sum \
    "${PREFIX}.fa" \
    "${PREFIX}.metadata.tsv" \
    "${PREFIX}.paf" \
    "$LIFT/gaidukov_mapping_qc.tsv" \
    "$LIFT/gaidukov_loci_mapped.tsv"

echo
echo "GAIDUKOV LIFTOVER JOB: PASS"
