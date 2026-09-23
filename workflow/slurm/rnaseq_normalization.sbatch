#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --job-name=rnaseq_normalization
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=rnaseq_normalization_%j.out

set -euo pipefail

PROJECT=/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_git
SOURCE=/cephyr/users/mayoa/Vera/hpc_rnaseq_production_results/2026-09-10/final_matrix/rnaseq_gene_counts_v1
EXPERIMENTS="$PROJECT/benchmarks/rnaseq/2026-09-22/rna_experiment_matrix_v1"
OUTDIR="$PROJECT/benchmarks/rnaseq/2026-09-23/rna_normalization_v1"
ENV_PREFIX="${RNASEQ_NORM_ENV_PREFIX:-$HOME/.conda/envs/tfm_rnaseq_normalization}"

module purge
module load Miniforge3/24.1.2-0

if [[ ! -x "$ENV_PREFIX/bin/Rscript" ]]; then
    echo "ERROR: $ENV_PREFIX/bin/Rscript is missing or not executable. Create the RNA normalization Conda environment before sbatch submission." >&2
    exit 2
fi

cd "$PROJECT"
"$ENV_PREFIX/bin/Rscript" workflow/scripts/normalize_rnaseq_experiments.R \
    --counts "$EXPERIMENTS/rnaseq_gene_counts_by_experiment_1473.tsv.gz" \
    --metadata "$EXPERIMENTS/experiment_metadata.tsv" \
    --qc "$EXPERIMENTS/experiment_raw_count_qc.tsv" \
    --annotation "$SOURCE/gene_annotation.tsv.gz" \
    --outdir "$OUTDIR"
