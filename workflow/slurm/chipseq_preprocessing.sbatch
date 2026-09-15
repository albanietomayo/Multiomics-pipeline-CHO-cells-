#!/bin/bash
#SBATCH --job-name=chipseq_preprocessing
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail
: "${TMPDIR:?Node-local TMPDIR is required}"
CHIP_SOURCE="$PWD"
CHIP_SAVE="$(dirname "$CHIP_SOURCE")/job_${SLURM_JOB_ID}"
CHIP_WORK="$TMPDIR/chipseq-pilot"
mkdir -p "$CHIP_SAVE" "$CHIP_WORK"

finish() {
    CHIP_EXIT=$?
    trap - EXIT
    set +e
    CHIP_WORKFLOW_EXIT=$CHIP_EXIT
    if [[ "$CHIP_EXIT" == 0 ]]; then
        rsync -a "$CHIP_WORK/data/" "$CHIP_SOURCE/data/" || CHIP_EXIT=$?
        rsync -a "$CHIP_WORK/results/" "$CHIP_SOURCE/results/" || CHIP_EXIT=$?
        cp "$CHIP_WORK/config/chipseq_pilot_runs.tsv" "$CHIP_SOURCE/config/" || CHIP_EXIT=$?
        if [[ "$CHIP_EXIT" == 0 ]]; then
            (cd "$CHIP_SOURCE" && sha256sum --check "$CHIP_SAVE/output.sha256") \
                > "$CHIP_SAVE/output_copy_validation.log" 2>&1 || CHIP_EXIT=$?
        fi
    elif [[ -d "$CHIP_WORK/results" ]]; then
        rsync -a --exclude='*.fastq.gz' --exclude='*.part' \
            "$CHIP_WORK/results/" "$CHIP_SOURCE/results/" \
            > "$CHIP_SAVE/failure_copy.log" 2>&1
    fi
    printf 'workflow_exit_status\t%s\nexit_status\t%s\n' \
        "$CHIP_WORKFLOW_EXIT" "$CHIP_EXIT" > "$CHIP_SAVE/job_status.tsv"
    if [[ "$CHIP_EXIT" == 0 ]]; then
        echo "[OK] Pilot completed; outputs copied and verified: $CHIP_SOURCE"
    else
        echo "[ERROR] Pilot or output copy failed; see $CHIP_SAVE"
    fi
    exit "$CHIP_EXIT"
}
trap finish EXIT

sha256sum --check ../source.sha256 > "$CHIP_SAVE/source_validation.log"
rsync -a --exclude=data --exclude=results --exclude=.snakemake \
    "$CHIP_SOURCE/" "$CHIP_WORK/"
cd "$CHIP_WORK"

module purge
module load Miniforge3/24.1.2-0
export CONDA_PKGS_DIRS="$TMPDIR/chipseq-conda-pkgs"
export XDG_CACHE_HOME="$TMPDIR/chipseq-cache"
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$CONDA_PKGS_DIRS" "$XDG_CACHE_HOME"
CHIP_ENV="$TMPDIR/cho-chipseq-pilot"
conda env create --prefix "$CHIP_ENV" --file environment.yml \
    > "$CHIP_SAVE/environment_creation.log" 2>&1
export PATH="$CHIP_ENV/bin:$PATH"
snakemake --version > "$CHIP_SAVE/snakemake_version.txt"
conda list --prefix "$CHIP_ENV" --explicit > "$CHIP_SAVE/environment_explicit.txt"

CHIP_ARGS=(
    --snakefile workflow/rules/chipseq_preprocessing.smk
    --cores "${SLURM_CPUS_PER_TASK:-4}"
    --use-conda
    --conda-prefix "$TMPDIR/chipseq-rule-envs"
    --rerun-incomplete
)

# ChIP eligibility integration v1
snakemake chipseq_experimental_eligibility "${CHIP_ARGS[@]}" -p \
    > "$CHIP_SAVE/eligibility_validation.log" 2>&1

snakemake results/metadata/validation_fastq_manifest.tsv "${CHIP_ARGS[@]}" -p \
    > "$CHIP_SAVE/manifest_setup.log" 2>&1
python workflow/scripts/prepare_chipseq_pilot_runs.py \
    --pilot config/chipseq_pilot.json --samples config/samples.tsv \
    --check-manifest results/metadata/validation_fastq_manifest.tsv \
    > "$CHIP_SAVE/manifest_validation.log" 2>&1

snakemake chipseq_preprocessing_all "${CHIP_ARGS[@]}" -n \
    > "$CHIP_SAVE/dry_run.log" 2>&1
/usr/bin/time -v -o "$CHIP_SAVE/workflow_resources.txt" \
    snakemake chipseq_preprocessing_all "${CHIP_ARGS[@]}" -p \
    > "$CHIP_SAVE/snakemake_execution.log" 2>&1

python - <<'PY'
import csv
import json
from pathlib import Path
import yaml

pilot = json.loads(Path("config/chipseq_pilot.json").read_text())
expected = {r["run_accession"] for r in pilot["runs"]}
config = yaml.safe_load(Path("config/config.yaml").read_text())
for key in ("qc_by_fastq", "qc_by_job"):
    path = Path(config["preprocessing"][key])
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 2 or {r["run_accession"] for r in rows} != expected:
        raise ValueError(f"QC summary does not match the pilot: {path}")
print("[OK] Both QC tables contain exactly the two pilot runs.")
PY

find data results config/chipseq_pilot_runs.tsv -type f -print0 |
    sort -z | xargs -0 sha256sum > "$CHIP_SAVE/output.sha256"
