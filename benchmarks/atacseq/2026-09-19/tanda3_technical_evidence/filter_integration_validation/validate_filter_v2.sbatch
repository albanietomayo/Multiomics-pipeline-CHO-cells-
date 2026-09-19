#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=atac_filter_check
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO="/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_git"
SAVE="$HOME/atac_filter_integration_validation/job_${SLURM_JOB_ID}"
SCRATCH="${TMPDIR:?TMPDIR is required}/atac_filter_validation_${SLURM_JOB_ID}"

mkdir -p "$SAVE" "$SCRATCH"

cd "$REPO"

echo "============================================================"
echo "REPOSITORY STATE"
echo "============================================================"

BRANCH="$(git branch --show-current)"

[[ "$BRANCH" == "atacseq-downstream-integration" ]] || {
    echo "ERROR: unexpected branch: $BRANCH" >&2
    exit 10
}

git status --short
git diff --check

git diff -- \
  config/config.yaml \
  workflow/rules/atacseq.smk \
  > "$SAVE/integration.diff"

echo
echo "============================================================"
echo "CREATE SNAKEMAKE ENVIRONMENT IN SCRATCH"
echo "============================================================"

module purge
module load Miniforge3/24.1.2-0

export CONDA_PKGS_DIRS="$SCRATCH/conda-pkgs"
export XDG_CACHE_HOME="$SCRATCH/cache"

mkdir -p "$CONDA_PKGS_DIRS" "$XDG_CACHE_HOME"

MAIN_ENV="$SCRATCH/main-env"

conda env create \
  --prefix "$MAIN_ENV" \
  --file environment.yml \
  > "$SAVE/environment_creation.log" 2>&1

SNAKEMAKE="$MAIN_ENV/bin/snakemake"

"$SNAKEMAKE" --version \
  > "$SAVE/snakemake_version.txt"

echo "Snakemake: $("$SNAKEMAKE" --version)"

echo
echo "============================================================"
echo "PARSE / RULE CHECK"
echo "============================================================"

grep -n '^rule atacseq_filter_bam:' \
  workflow/rules/atacseq.smk

if grep -q '^rule atacseq_filter_single_bam:' \
    workflow/rules/atacseq.smk
then
    echo "ERROR: old SINGLE-only rule still exists" >&2
    exit 20
fi

echo "Unified filter rule: present"
echo "Old SINGLE-only rule: absent"

echo
echo "============================================================"
echo "DRY-RUN FILTERED BAM"
echo "============================================================"

for RUN in SRR12774931 SRR12774932; do

    TARGET="results/atacseq/filtered_bam/${RUN}/filtered.bam"

    echo
    echo "===== $RUN : FILTER ====="

    "$SNAKEMAKE" \
      --snakefile Snakefile \
      --cores 1 \
      --dry-run \
      --printshellcmds \
      --forcerun atacseq_filter_bam \
      "$TARGET" \
      > "$SAVE/${RUN}.filter.dryrun.log" 2>&1

    echo "FILTER DRY-RUN: OK"

    grep -E \
      'rule atacseq_filter_bam|run_accession:|structure|SRR1277493[12]' \
      "$SAVE/${RUN}.filter.dryrun.log" \
      | head -n 80 || true

done

echo
echo "============================================================"
echo "DRY-RUN THROUGH TSS"
echo "============================================================"

for RUN in SRR12774931 SRR12774932; do

    TARGET="results/atacseq/qc/${RUN}/tss_enrichment_summary.tsv"

    echo
    echo "===== $RUN : TSS ====="

    "$SNAKEMAKE" \
      --snakefile Snakefile \
      --cores 1 \
      --dry-run \
      --printshellcmds \
      "$TARGET" \
      > "$SAVE/${RUN}.tss.dryrun.log" 2>&1

    echo "TSS DRY-RUN: OK"

    grep -E \
      'atacseq_filter_bam|atacseq_tss_enrichment|SRR1277493[12]' \
      "$SAVE/${RUN}.tss.dryrun.log" \
      | head -n 100 || true

done

echo
echo "============================================================"
echo "FINAL VALIDATION"
echo "============================================================"

printf \
  "Unified ATAC filtering rule parsed successfully and both eligible SINGLE CHO runs resolved through filtered BAM and TSS dry-runs.\n" \
  > "$SAVE/complete.ok"

cat "$SAVE/complete.ok"

echo
echo "RESULTS:"
echo "$SAVE"
