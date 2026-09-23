#!/usr/bin/env bash
#SBATCH --job-name=rnaseq_study_aggregation
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=rnaseq_study_aggregation_%j.out

set -euo pipefail

REPO="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"

cd "$REPO"

LOGCPM="benchmarks/rnaseq/2026-09-23/rna_normalization_v1/rnaseq_tmm_logcpm_by_experiment_1473.tsv.gz"
TPM="benchmarks/rnaseq/2026-09-23/rna_normalization_v1/rnaseq_tpm_by_experiment_1473.tsv.gz"
META="benchmarks/rnaseq/2026-09-22/rna_experiment_matrix_v1/experiment_metadata.tsv"
OUT="benchmarks/rnaseq/2026-09-23/rna_study_aggregation_v1"

echo "RNA study-level aggregation"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Host: $(hostname)"
echo "Repository: $REPO"
echo "Git commit: $(git rev-parse HEAD)"
echo

for f in "$LOGCPM" "$TPM" "$META"; do
    [[ -s "$f" ]] || {
        echo "ERROR: missing or empty input: $f" >&2
        exit 1
    }
done

[[ ! -e "$OUT" ]] || {
    echo "ERROR: output already exists: $OUT" >&2
    exit 1
}

python3 workflow/scripts/aggregate_rnaseq_by_study.py \
    --logcpm "$LOGCPM" \
    --tpm "$TPM" \
    --metadata "$META" \
    --outdir "$OUT"

echo
echo "Production aggregation completed"
