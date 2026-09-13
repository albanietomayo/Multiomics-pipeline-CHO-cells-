#!/bin/bash
#SBATCH --job-name=chipseq_metadata
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
: "${TMPDIR:?The job requires node-local TMPDIR}"

PROJECT="$PWD"
SAVE="$(dirname "$PROJECT")/job_${SLURM_JOB_ID}"
mkdir -p "$SAVE"
trap 'CHIP_EXIT=$?; printf "exit_status\t%s\n" "$CHIP_EXIT" > "$SAVE/job_status.tsv"' EXIT

exec 9>results/metadata/chipseq/audit.lock
flock -n 9 || { echo "Another ChIP metadata job is active."; exit 1; }
sha256sum --check ../source.sha256 > "$SAVE/source_validation.log"

module purge
module load Miniforge3/24.1.2-0
export CONDA_PKGS_DIRS="$TMPDIR/chipseq-conda-pkgs"
export XDG_CACHE_HOME="$TMPDIR/chipseq-cache"
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$CONDA_PKGS_DIRS" "$XDG_CACHE_HOME"
CHIP_ENV="$TMPDIR/cho-chipseq-metadata"

conda env create --prefix "$CHIP_ENV" --file environment.yml \
    > "$SAVE/environment_creation.log" 2>&1
export PATH="$CHIP_ENV/bin:$PATH"
snakemake --version > "$SAVE/snakemake_version.txt"
conda list --prefix "$CHIP_ENV" --explicit > "$SAVE/environment_explicit.txt"

CHIP_ARGS=(
    --snakefile workflow/rules/chipseq_metadata.smk
    --cores 1
    --rerun-incomplete
    --forcerun chipseq_xml_inventory
)

snakemake chipseq_metadata_all "${CHIP_ARGS[@]}" -n \
    > "$SAVE/dry_run.log" 2>&1

/usr/bin/time -v -o "$SAVE/workflow_resources.txt" \
    snakemake chipseq_metadata_all "${CHIP_ARGS[@]}" --keep-going -p \
    > "$SAVE/snakemake_execution.log" 2>&1

echo "[OK] Workflow completed. Logs: $SAVE"
