#!/bin/bash
# --check prepares and validates an isolated submission; --submit also sends it.
set -euo pipefail
CHIP_MODE="${1:---check}"
if [[ $# -gt 1 || ( "$CHIP_MODE" != --check && "$CHIP_MODE" != --submit ) ]]; then
    echo "Usage: bash workflow/slurm/submit_chipseq_preprocessing.sh [--check|--submit]" >&2
    exit 2
fi
CHIP_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$CHIP_ROOT"
[[ "$(git branch --show-current)" == chipseq-metadata-audit ]] || {
    echo "ERROR: expected branch chipseq-metadata-audit" >&2
    exit 1
}
CHIP_FILES=(
    environment.yml config/config.yaml config/samples.tsv config/chipseq_pilot.json
    config/chipseq_experimental_eligibility.json
    workflow/rules/common.smk workflow/rules/chipseq_preprocessing.smk
    workflow/rules/chipseq_eligibility.smk
    workflow/scripts/build_fastq_manifest.py workflow/scripts/build_validation_manifest.py
    workflow/scripts/download_fastq.py workflow/scripts/download_validation_subset.py
    workflow/scripts/summarize_preprocessing_qc.py workflow/scripts/prepare_chipseq_pilot_runs.py
    workflow/scripts/validate_chipseq_eligibility.py
    workflow/envs/qc.yaml workflow/envs/preprocessing.yaml workflow/envs/reporting.yaml
    workflow/slurm/chipseq_preprocessing.sbatch
    workflow/slurm/submit_chipseq_preprocessing.sh
    snapshots/chipseq/control_validation_001/chipseq_conditions.tsv
    snapshots/chipseq/control_validation_001/control_candidate_summary.json
    snapshots/chipseq/control_validation_001/request_plan_provenance.json
)
for CHIP_FILE in "${CHIP_FILES[@]}"; do
    [[ -f "$CHIP_FILE" ]] || { echo "ERROR: missing $CHIP_FILE" >&2; exit 1; }
done
bash -n workflow/slurm/chipseq_preprocessing.sbatch
CHIP_SLURM="$CHIP_ROOT/results/chipseq/preprocessing/slurm"
mkdir -p "$CHIP_SLURM"
CHIP_SUBMIT=$(mktemp -d "$CHIP_SLURM/submission_XXXXXX")
for CHIP_FILE in "${CHIP_FILES[@]}"; do
    mkdir -p "$CHIP_SUBMIT/project/$(dirname "$CHIP_FILE")"
    cp -p -- "$CHIP_FILE" "$CHIP_SUBMIT/project/$CHIP_FILE"
done
git rev-parse HEAD > "$CHIP_SUBMIT/base_commit.txt"
git status --short > "$CHIP_SUBMIT/git_status.txt"
git diff --binary -- "${CHIP_FILES[@]}" > "$CHIP_SUBMIT/source_changes.patch"
(
    cd "$CHIP_SUBMIT/project"
    sha256sum -- "${CHIP_FILES[@]}" > ../source.sha256
    python3 workflow/scripts/validate_chipseq_eligibility.py \
        --decisions config/chipseq_experimental_eligibility.json \
        --samples config/samples.tsv \
        --conditions snapshots/chipseq/control_validation_001/chipseq_conditions.tsv \
        --pilot config/chipseq_pilot.json --report ../eligibility_preflight.json
    python3 workflow/scripts/prepare_chipseq_pilot_runs.py \
        --pilot config/chipseq_pilot.json --samples config/samples.tsv \
        --output ../pilot_selection_check.tsv
    sha256sum --check ../source.sha256 > ../source_validation.log
)
if [[ "$CHIP_MODE" == --check ]]; then
    echo "[OK] Source copy, eligibility and pilot selection checked. No job submitted."
    echo "CHECK DIRECTORY: $CHIP_SUBMIT"
    echo "Snakemake DAG validation remains pending in the VERA environment."
    exit 0
fi
CHIP_SUBMISSION=$(sbatch --parsable --chdir="$CHIP_SUBMIT/project" \
    --output="$CHIP_SUBMIT/slurm_%j.out" \
    "$CHIP_SUBMIT/project/workflow/slurm/chipseq_preprocessing.sbatch")
CHIP_JOB="${CHIP_SUBMISSION%%;*}"
[[ "$CHIP_JOB" =~ ^[0-9]+$ ]] || { echo "ERROR: unexpected sbatch response: $CHIP_SUBMISSION" >&2; exit 1; }
printf '%s\n' "$CHIP_JOB" > "$CHIP_SUBMIT/job_id.txt"
printf '%s\n' "$CHIP_SUBMIT" > "$CHIP_SLURM/latest_submission.txt"
echo "JOB SUBMITTED: $CHIP_JOB"
echo "LOGS: $CHIP_SUBMIT"
squeue -j "$CHIP_JOB" -o "%.18i %.26j %.12T %.10M %.6C %R" || true
