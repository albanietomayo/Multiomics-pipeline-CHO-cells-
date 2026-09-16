#!/bin/bash
#SBATCH --job-name=chipseq_filtering
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
set -euo pipefail
: "${TMPDIR:?Node-local TMPDIR is required}"
CHIP_SOURCE="$PWD"
CHIP_SAVE="$(dirname "$CHIP_SOURCE")/job_${SLURM_JOB_ID}"
CHIP_WORK="$TMPDIR/chipseq-filtering"
mkdir -p "$CHIP_SAVE" "$CHIP_WORK"
CHIP_STAGE=source_validation
finish() {
    CHIP_EXIT=$?
    trap - EXIT
    set +e
    if [[ -d "$CHIP_WORK/results/chipseq/filtering" && "$CHIP_EXIT" != 0 ]]; then
        rsync -a --exclude='*.bam' --exclude='*.csi' --exclude='*.part' \
            "$CHIP_WORK/results/chipseq/filtering/" "$CHIP_SAVE/failure_reports/" \
            > "$CHIP_SAVE/failure_copy.log" 2>&1
    fi
    if [[ -d "$CHIP_WORK/.snakemake/log" ]]; then
        rsync -a "$CHIP_WORK/.snakemake/log/" "$CHIP_SAVE/snakemake_logs/"
    fi
    printf 'stage\t%s\nexit_status\t%s\n' "$CHIP_STAGE" "$CHIP_EXIT" > "$CHIP_SAVE/job_status.tsv"
    echo "Filtering exit=$CHIP_EXIT; stage=$CHIP_STAGE; logs=$CHIP_SAVE"
    exit "$CHIP_EXIT"
}
trap finish EXIT
sha256sum --check ../source.sha256 > "$CHIP_SAVE/source_validation.log"
rsync -a "$CHIP_SOURCE/" "$CHIP_WORK/"
cd "$CHIP_WORK"
module purge
module load Miniforge3/24.1.2-0
export CONDA_PKGS_DIRS="$TMPDIR/chipseq-filtering-conda-pkgs"
export XDG_CACHE_HOME="$TMPDIR/chipseq-filtering-cache"
export PYTHONDONTWRITEBYTECODE=1
CHIP_ENV="$TMPDIR/chipseq-filtering-env"
CHIP_RULE_ENVS="$TMPDIR/chipseq-filtering-rule-envs"
CHIP_STAGE=environment_creation
conda env create --prefix "$CHIP_ENV" --file environment.yml > "$CHIP_SAVE/environment_creation.log" 2>&1
export PATH="$CHIP_ENV/bin:$PATH"
snakemake --version > "$CHIP_SAVE/snakemake_version.txt"
conda list --prefix "$CHIP_ENV" --explicit > "$CHIP_SAVE/environment_explicit.txt"
CHIP_ARGS=(--snakefile workflow/rules/chipseq_filtering.smk --cores "$SLURM_CPUS_PER_TASK"
    --resources mem_mb=20000 --use-conda --conda-prefix "$CHIP_RULE_ENVS" --rerun-incomplete)
CHIP_STAGE=eligibility
snakemake chipseq_experimental_eligibility "${CHIP_ARGS[@]}" -p > "$CHIP_SAVE/eligibility.log" 2>&1
CHIP_STAGE=dry_run
snakemake chipseq_filtering_all "${CHIP_ARGS[@]}" -n > "$CHIP_SAVE/dry_run.log" 2>&1
CHIP_STAGE=filtering
/usr/bin/time -v -o "$CHIP_SAVE/workflow_resources.txt" \
    snakemake chipseq_filtering_all "${CHIP_ARGS[@]}" -p > "$CHIP_SAVE/snakemake_execution.log" 2>&1
CHIP_STAGE=environment_provenance
mkdir -p "$CHIP_SAVE/rule_environments"
for CHIP_RULE_ENV in "$CHIP_RULE_ENVS"/*; do
    if [[ -d "$CHIP_RULE_ENV/conda-meta" ]]; then
        conda list --prefix "$CHIP_RULE_ENV" --explicit \
            > "$CHIP_SAVE/rule_environments/$(basename "$CHIP_RULE_ENV").txt"
    fi
done
CHIP_STAGE=output_copy
python workflow/scripts/chipseq_filtering_support.py publish --save "$CHIP_SAVE" \
    > "$CHIP_SAVE/output_copy_validation.log" 2>&1
CHIP_STAGE=completed
echo "[OK] Duplicate marking, nuclear filtering and QC completed; outputs copied and verified."
