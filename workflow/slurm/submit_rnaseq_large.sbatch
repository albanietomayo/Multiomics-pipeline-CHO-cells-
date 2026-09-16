#!/bin/bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --job-name=rnaseq_LARGE
#SBATCH --array=1-9%4
#SBATCH --cpus-per-task=16
#SBATCH --mem=94800M
#SBATCH --time=12:00:00
#SBATCH --output=/cephyr/users/mayoa/Vera/hpc_rnaseq_production_results/2026-09-10/slurm/rnaseq_LARGE_%A_%a.out

set -euo pipefail

BATCH_ID="$(printf 'LARGE_%03d' "$SLURM_ARRAY_TASK_ID")"
export BATCH_ID

export RNASEQ_SAVE_ROOT="$HOME/hpc_rnaseq_production_results/2026-09-10"
export RNASEQ_PLAN_REL="benchmarks/rnaseq/2026-09-10/rnaseq_batches.tsv"
export MEM_MB=40000

echo "============================================================"
echo "Starting $BATCH_ID"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"
echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
echo "CPUS=$SLURM_CPUS_PER_TASK"
echo "TMPDIR=$TMPDIR"
echo "============================================================"

exec bash "$HOME/TFM_multiomics_pipeline_git/workflow/slurm/rnaseq_batch_worker.sbatch"
