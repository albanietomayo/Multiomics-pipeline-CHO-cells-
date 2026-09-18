#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --job-name=cho_bm_lift
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=benchmark/reports/dhiman2020/slurm_%j.out
#SBATCH --error=benchmark/reports/dhiman2020/slurm_%j.err

set -euo pipefail

WT="/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_benchmark"

MINIMAP2="$WT/benchmark/raw/tools/minimap2-2.31_x64-linux/minimap2"

SOURCE_FASTA="$WT/benchmark/raw/reference/GCF_003668045.1/genome.fa"
TARGET_FASTA="$WT/resources/reference/CriGri-PICRH-1.0/genome.fa"

LOCI="$WT/benchmark/intermediate/dhiman2020/dhiman_unique_loci.tsv"

LIFT="$WT/benchmark/intermediate/dhiman2020/liftover"
PREFIX="$LIFT/dhiman_boundaries"

cd "$WT"

mkdir -p "$LIFT"

echo "============================================================"
echo "ENVIRONMENT"
echo "============================================================"

echo "Host:       $(hostname)"
echo "Date UTC:   $(date -u +%FT%TZ)"
echo "Git commit: $(git rev-parse HEAD)"

[[ -x "$MINIMAP2" ]] || {
    echo "ERROR: pinned minimap2 unavailable"
    exit 20
}

echo "minimap2:   $MINIMAP2"
echo "version:    $("$MINIMAP2" --version)"

echo
echo "============================================================"
echo "BUILD SOURCE FLANKS"
echo "============================================================"

python3 benchmark/scripts/build_dhiman_liftover_queries.py \
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

python3 benchmark/scripts/resolve_dhiman_liftover.py \
    "$LOCI" \
    "${PREFIX}.metadata.tsv" \
    "${PREFIX}.paf" \
    "$LIFT"

echo
echo "============================================================"
echo "MAPPED LOCI"
echo "============================================================"

column -t -s $'\t' \
    "$LIFT/dhiman_loci_mapped.tsv"

echo
echo "============================================================"
echo "CHECKSUMS"
echo "============================================================"

sha256sum \
    "${PREFIX}.fa" \
    "${PREFIX}.metadata.tsv" \
    "${PREFIX}.paf" \
    "$LIFT/dhiman_boundary_mapping_qc.tsv" \
    "$LIFT/dhiman_loci_mapped.tsv"

echo
echo "LIFTOVER JOB: PASS"
