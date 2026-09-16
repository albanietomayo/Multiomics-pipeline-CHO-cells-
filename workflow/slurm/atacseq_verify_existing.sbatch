#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=atac_verify_existing
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
REPO="${ATAC_PROJECT_HOME:?set ATAC_PROJECT_HOME}"
EXISTING="${ATAC_EXISTING_JOB:?set ATAC_EXISTING_JOB}"
WORK="${TMPDIR:?TMPDIR required}/atac_verify_${SLURM_JOB_ID:?}"
SAVE="${ATAC_VERIFY_SAVE:?set ATAC_VERIFY_SAVE}/job_${SLURM_JOB_ID}"
[[ ! -e "$WORK" && ! -e "$SAVE" ]] || exit 2
mkdir -p "$WORK" "$SAVE"
trap 'rc=$?; printf "exit_status\t%s\n" "$rc" > "$SAVE/job_status.tsv"; exit "$rc"' EXIT
for ITEM in config workflow environment.yml Snakefile.atac; do cp -a "$REPO/$ITEM" "$WORK/"; done
git -C "$REPO" rev-parse HEAD > "$SAVE/source_commit.txt"
git -C "$REPO" diff HEAD > "$SAVE/source_changes.patch"
cd "$WORK"
find config workflow -type f ! -name '*.pyc' -print0 | sort -z | xargs -0 sha256sum > "$SAVE/source_SHA256SUMS.txt"
module purge
module load Miniforge3/24.1.2-0
export CONDA_PKGS_DIRS="$TMPDIR/atac-verify-pkgs"
conda env create --prefix "$TMPDIR/atac-verify-main" --file environment.yml > "$SAVE/environment.log" 2>&1
export PATH="$TMPDIR/atac-verify-main/bin:$PATH"
snakemake --snakefile Snakefile.atac atacseq_existing_all --use-conda \
  --conda-prefix "$TMPDIR/atac-verify-envs" --cores 2 \
  --default-resources "tmpdir=$TMPDIR" \
  --config "atacseq_existing_job=$EXISTING" -p > "$SAVE/snakemake_execution.log" 2>&1
cp -a results/atacseq/qc/verification_existing "$SAVE/"
(cd "$SAVE" && find verification_existing -type f -print0 | sort -z | xargs -0 sha256sum) > "$SAVE/SHA256SUMS.txt"
printf 'Technical checks passed; biological acceptance not assessed.\n' > "$SAVE/complete.ok"
