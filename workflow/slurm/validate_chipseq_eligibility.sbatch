#!/bin/bash
#SBATCH --job-name=chipseq_eligibility_check
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
: "${TMPDIR:?Node-local TMPDIR is required}"
CHIP_SOURCE="$PWD"
CHIP_SAVE="$(dirname "$CHIP_SOURCE")/job_${SLURM_JOB_ID}"
CHIP_WORK="$TMPDIR/chipseq-eligibility-validation"
CHIP_STAGE=initialization
mkdir -p "$CHIP_SAVE" "$CHIP_WORK"

finish() {
    CHIP_EXIT=$?
    trap - EXIT
    set +e
    if [[ -d "$CHIP_WORK/results" ]]; then
        rsync -a "$CHIP_WORK/results/" "$CHIP_SAVE/results/" \
            > "$CHIP_SAVE/result_copy.log" 2>&1
        CHIP_COPY_EXIT=$?
        if [[ "$CHIP_EXIT" == 0 && "$CHIP_COPY_EXIT" != 0 ]]; then
            CHIP_EXIT=$CHIP_COPY_EXIT
            CHIP_STAGE=copy_results
        fi
    fi
    printf 'stage\t%s\nexit_status\t%s\n' \
        "$CHIP_STAGE" "$CHIP_EXIT" > "$CHIP_SAVE/validation_status.tsv"
    echo "Validation exit=$CHIP_EXIT; stage=$CHIP_STAGE; logs=$CHIP_SAVE"
    exit "$CHIP_EXIT"
}
trap finish EXIT

CHIP_STAGE=source_validation
sha256sum --check ../source.sha256 > "$CHIP_SAVE/source_validation.log"
rsync -a "$CHIP_SOURCE/" "$CHIP_WORK/"
cd "$CHIP_WORK"

module purge
module load Miniforge3/24.1.2-0
export CONDA_PKGS_DIRS="$TMPDIR/chipseq-validation-conda-pkgs"
export XDG_CACHE_HOME="$TMPDIR/chipseq-validation-cache"
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$CONDA_PKGS_DIRS" "$XDG_CACHE_HOME"

CHIP_STAGE=environment_creation
CHIP_ENV="$TMPDIR/chipseq-validation-env"
conda env create --prefix "$CHIP_ENV" --file environment.yml \
    > "$CHIP_SAVE/environment_creation.log" 2>&1
export PATH="$CHIP_ENV/bin:$PATH"
snakemake --version > "$CHIP_SAVE/snakemake_version.txt"
conda list --prefix "$CHIP_ENV" --explicit \
    > "$CHIP_SAVE/environment_explicit.txt"

CHIP_ARGS=(
    --cores "${SLURM_CPUS_PER_TASK:-2}"
    --use-conda
    --conda-prefix "$TMPDIR/chipseq-validation-rule-envs"
)
CHIP_PRE=(--snakefile workflow/rules/chipseq_preprocessing.smk)

CHIP_STAGE=pilot_eligibility
snakemake chipseq_experimental_eligibility \
    "${CHIP_PRE[@]}" "${CHIP_ARGS[@]}" -p \
    > "$CHIP_SAVE/eligibility_execution.log" 2>&1

CHIP_STAGE=pilot_manifest
snakemake results/metadata/validation_fastq_manifest.tsv \
    "${CHIP_PRE[@]}" "${CHIP_ARGS[@]}" -p \
    > "$CHIP_SAVE/manifest_setup.log" 2>&1
python workflow/scripts/prepare_chipseq_pilot_runs.py \
    --pilot config/chipseq_pilot.json --samples config/samples.tsv \
    --check-manifest results/metadata/validation_fastq_manifest.tsv \
    > "$CHIP_SAVE/manifest_validation.log" 2>&1

CHIP_STAGE=preprocessing_dry_run
snakemake chipseq_preprocessing_all \
    "${CHIP_PRE[@]}" "${CHIP_ARGS[@]}" -n \
    > "$CHIP_SAVE/preprocessing_dry_run.log" 2>&1

CHIP_STAGE=metadata_configuration
python - <<'PY'
from pathlib import Path
import yaml

config = yaml.safe_load(Path("config/chipseq_metadata.yaml").read_text())
config["chipseq_metadata"].update({
    "samples": "config/samples.tsv",
    "plan_dir": "results/chipseq/eligibility_validation/plan",
    "snapshot_dir": "results/chipseq/eligibility_validation/xml",
})
Path("config/chipseq_metadata_validation.yaml").write_text(
    yaml.safe_dump(config, sort_keys=False)
)
PY
cp config/chipseq_metadata_validation.yaml "$CHIP_SAVE/"

CHIP_META=(
    --snakefile workflow/rules/chipseq_metadata.smk
    --configfile config/chipseq_metadata_validation.yaml
)

CHIP_STAGE=metadata_plan
snakemake results/chipseq/eligibility_validation/plan/ena_request_plan.tsv \
    "${CHIP_META[@]}" "${CHIP_ARGS[@]}" -p \
    > "$CHIP_SAVE/metadata_plan.log" 2>&1

CHIP_STAGE=metadata_dry_run
snakemake chipseq_metadata_all \
    "${CHIP_META[@]}" "${CHIP_ARGS[@]}" -n \
    > "$CHIP_SAVE/metadata_dry_run.log" 2>&1

CHIP_STAGE=completed
echo "[OK] Eligibility executed; pilot manifest checked; both DAG dry-runs passed."
echo "[INFO] No FASTQ download, preprocessing or XML retrieval was executed."
