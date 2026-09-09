#!/usr/bin/env bash
#SBATCH -A C3SE2026-1-42
#SBATCH -p vera
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --job-name=rnaseq_smoke
#SBATCH --output=/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline/hpc_logs/rnaseq_smoke_%j.out

set -uo pipefail

PROJECT_HOME="$HOME/TFM_multiomics_pipeline"
WORK="$TMPDIR/TFM_multiomics_pipeline"
SAVE="$PROJECT_HOME/hpc_pilot_results/$SLURM_JOB_ID"
MAIN_ENV="$TMPDIR/cho-multiomics"

echo "======================================"
echo " RNA-seq HPC SMOKE TEST"
echo "======================================"
echo "Run: SRR10493862"
echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "TMPDIR: $TMPDIR"
echo "Start: $(date)"
echo

mkdir -p "$WORK" "$SAVE"

echo "=== COPYING PROJECT TO LOCAL SCRATCH ==="
cp -a "$PROJECT_HOME"/. "$WORK"/
cd "$WORK"

echo
echo "=== INITIAL TMPDIR SPACE ==="
df -h "$TMPDIR"

echo
echo "=== SOFTWARE SETUP ==="
module purge
module load Miniforge3/24.1.2-0

export CONDA_PKGS_DIRS="$TMPDIR/conda-pkgs"
export XDG_CACHE_HOME="$TMPDIR/.cache"
mkdir -p "$CONDA_PKGS_DIRS" "$XDG_CACHE_HOME"

echo "Creating main Snakemake environment..."

conda env create \
    --prefix "$MAIN_ENV" \
    --file environment.yml

if [ $? -ne 0 ]; then
    echo "ERROR: main Conda environment could not be created."
    exit 1
fi

SNAKEMAKE="$MAIN_ENV/bin/snakemake"

echo
echo "Snakemake version:"
"$SNAKEMAKE" --version

echo
echo "=== STARTING WORKFLOW ==="
echo "Target: SRR10493862 gene counts"
echo

"$SNAKEMAKE" \
    results/rnaseq/counts_by_run/SRR10493862/gene_counts.tsv \
    --use-conda \
    --conda-prefix "$TMPDIR/snakemake-conda" \
    --cores 8 \
    --resources mem_mb=80000 \
    --rerun-incomplete \
    --notemp \
    -p

STATUS=$?

echo
echo "=== WORKFLOW EXIT STATUS: $STATUS ==="

{
    echo "RNA-seq HPC smoke test"
    echo "Run: SRR10493862"
    echo "Job ID: $SLURM_JOB_ID"
    echo "Host: $(hostname)"
    echo "CPUs: $SLURM_CPUS_PER_TASK"
    echo "Workflow exit status: $STATUS"
    echo "End: $(date)"
    echo
    echo "=== TMPDIR ==="
    df -h "$TMPDIR"
    echo
    echo "=== PROJECT SIZE ==="
    du -sh "$WORK" 2>/dev/null || true
    echo
    echo "=== REFERENCE SIZE ==="
    du -sh resources/reference 2>/dev/null || true
    echo
    echo "=== RNA-seq RESULTS SIZE ==="
    du -sh results/rnaseq 2>/dev/null || true
    echo
    echo "=== BAM FILES ==="
    find results/rnaseq -type f -name '*.bam' -printf '%p\t%s bytes\n' 2>/dev/null || true
    echo
    echo "=== GENE COUNTS ==="
    ls -lh results/rnaseq/counts_by_run/SRR10493862/gene_counts.tsv 2>/dev/null || true
} > "$SAVE/resource_summary.txt"

tar \
    --exclude='*.bam' \
    --exclude='*.bai' \
    -czf "$SAVE/rnaseq_small_outputs.tar.gz" \
    results/rnaseq \
    2>/dev/null || true

cp -a .snakemake/log "$SAVE/snakemake_logs" 2>/dev/null || true

echo
echo "Small outputs saved in:"
echo "$SAVE"

echo
echo "End: $(date)"

exit "$STATUS"
