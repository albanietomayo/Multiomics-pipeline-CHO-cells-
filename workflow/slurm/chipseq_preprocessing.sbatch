#!/bin/bash
#SBATCH --job-name=chipseq_preprocessing
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00

set -euo pipefail

: "${TMPDIR:?Node-local TMPDIR is required}"

CHIP_SOURCE="$PWD"
CHIP_SAVE="$(dirname "$CHIP_SOURCE")/job_${SLURM_JOB_ID}"
CHIP_WORK="$TMPDIR/chipseq-incremental"

PROCESSING_RUNS="results/metadata/chipseq/ena_snapshot_001/analysis_plan/chipseq_processing_runs.tsv"
PROCESSING_SUMMARY="results/metadata/chipseq/ena_snapshot_001/analysis_plan/processing_run_summary.json"

mkdir -p \
    "$CHIP_SAVE" \
    "$CHIP_WORK"

finish() {
    CHIP_EXIT=$?

    trap - EXIT
    set +e

    CHIP_WORKFLOW_EXIT=$CHIP_EXIT

    if [[ "$CHIP_EXIT" == 0 ]]; then

        rsync -a \
            "$CHIP_WORK/data/" \
            "$CHIP_SOURCE/data/" \
            || CHIP_EXIT=$?

        rsync -a \
            "$CHIP_WORK/results/" \
            "$CHIP_SOURCE/results/" \
            || CHIP_EXIT=$?

        if [[ "$CHIP_EXIT" == 0 ]]; then
            (
                cd "$CHIP_SOURCE"
                sha256sum \
                    --check "$CHIP_SAVE/output.sha256"
            ) > \
                "$CHIP_SAVE/output_copy_validation.log" \
                2>&1 \
                || CHIP_EXIT=$?
        fi

    elif [[ -d "$CHIP_WORK/results" ]]; then

        rsync -a \
            --exclude='*.fastq.gz' \
            --exclude='*.part' \
            "$CHIP_WORK/results/" \
            "$CHIP_SOURCE/results/" \
            > "$CHIP_SAVE/failure_copy.log" \
            2>&1
    fi

    printf \
        'workflow_exit_status\t%s\nexit_status\t%s\n' \
        "$CHIP_WORKFLOW_EXIT" \
        "$CHIP_EXIT" \
        > "$CHIP_SAVE/job_status.tsv"

    if [[ "$CHIP_EXIT" == 0 ]]; then
        echo \
            "[OK] Incremental ChIP preprocessing completed; outputs copied and verified."
    else
        echo \
            "[ERROR] Incremental ChIP preprocessing failed; see $CHIP_SAVE"
    fi

    exit "$CHIP_EXIT"
}

trap finish EXIT

sha256sum \
    --check ../source.sha256 \
    > "$CHIP_SAVE/source_validation.log"

rsync -a \
    --exclude=data \
    --exclude=results \
    --exclude=.snakemake \
    "$CHIP_SOURCE/" \
    "$CHIP_WORK/"

cd "$CHIP_WORK"

module purge
module load Miniforge3/24.1.2-0

export CONDA_PKGS_DIRS="$TMPDIR/chipseq-conda-pkgs"
export XDG_CACHE_HOME="$TMPDIR/chipseq-cache"
export PYTHONDONTWRITEBYTECODE=1

mkdir -p \
    "$CONDA_PKGS_DIRS" \
    "$XDG_CACHE_HOME"

CHIP_ENV="$TMPDIR/cho-chipseq-incremental"

conda env create \
    --prefix "$CHIP_ENV" \
    --file environment.yml \
    > "$CHIP_SAVE/environment_creation.log" \
    2>&1

export PATH="$CHIP_ENV/bin:$PATH"

snakemake --version \
    > "$CHIP_SAVE/snakemake_version.txt"

python --version \
    > "$CHIP_SAVE/python_version.txt" \
    2>&1

python - <<'PY' \
    > "$CHIP_SAVE/pandas_version.txt"
import pandas
print(pandas.__version__)
PY

conda list \
    --prefix "$CHIP_ENV" \
    --explicit \
    > "$CHIP_SAVE/environment_explicit.txt"

META_ARGS=(
    --snakefile workflow/rules/chipseq_metadata.smk
    --cores "${SLURM_CPUS_PER_TASK:-4}"
    --rerun-incomplete
)

PREP_ARGS=(
    --snakefile workflow/rules/chipseq_preprocessing.smk
    --cores "${SLURM_CPUS_PER_TASK:-4}"
    --use-conda
    --conda-prefix "$TMPDIR/chipseq-rule-envs"
    --rerun-incomplete
)

echo "============================================================" \
    > "$CHIP_SAVE/runtime_chain.txt"

echo "1 metadata/planning" \
    >> "$CHIP_SAVE/runtime_chain.txt"

echo "2 FASTQ manifest + preprocessing/QC" \
    >> "$CHIP_SAVE/runtime_chain.txt"

echo "============================================================" \
    >> "$CHIP_SAVE/runtime_chain.txt"

# ------------------------------------------------------------
# Stage 1: active metadata -> incremental eligibility ->
#          analysis plan -> deduplicated processing-run set
# ------------------------------------------------------------

snakemake \
    chipseq_metadata_all \
    "${META_ARGS[@]}" \
    -n \
    > "$CHIP_SAVE/metadata_dry_run.log" \
    2>&1

/usr/bin/time \
    -v \
    -o "$CHIP_SAVE/metadata_resources.txt" \
    snakemake \
        chipseq_metadata_all \
        "${META_ARGS[@]}" \
        -p \
        > "$CHIP_SAVE/metadata_execution.log" \
        2>&1

[[ -s "$PROCESSING_RUNS" ]] || {
    echo "ERROR: processing run table was not generated" >&2
    exit 1
}

[[ -s "$PROCESSING_SUMMARY" ]] || {
    echo "ERROR: processing summary was not generated" >&2
    exit 1
}

python - \
    "$PROCESSING_RUNS" \
    "$PROCESSING_SUMMARY" <<'PY' \
    > "$CHIP_SAVE/processing_plan_validation.log"
import csv
import json
import sys
from collections import Counter

runs_path = sys.argv[1]
summary_path = sys.argv[2]

with open(
    runs_path,
    encoding="utf-8",
    newline="",
) as handle:
    rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

with open(
    summary_path,
    encoding="utf-8",
) as handle:
    summary = json.load(handle)

if not rows:
    raise ValueError(
        "Dynamic ChIP processing cohort is empty"
    )

accessions = [
    row["run_accession"]
    for row in rows
]

if len(accessions) != len(set(accessions)):
    raise ValueError(
        "Duplicate run accessions in processing cohort"
    )

roles = Counter(
    row["library_role"]
    for row in rows
)

if set(roles) - {"ip", "input"}:
    raise ValueError(
        f"Unexpected processing roles: {dict(roles)}"
    )

if roles["ip"] < 1:
    raise ValueError(
        "Processing cohort contains no IP runs"
    )

if roles["input"] < 1:
    raise ValueError(
        "Processing cohort contains no Input runs"
    )

if summary["unique_processing_runs"] != len(rows):
    raise ValueError(
        "Processing summary/run-table count mismatch"
    )

if summary["ready_analysis_count"] < 1:
    raise ValueError(
        "No ready ChIP analyses in processing summary"
    )

print(
    f"ready_analyses={summary['ready_analysis_count']}"
)

print(
    f"unique_processing_runs={len(rows)}"
)

print(
    f"ip_runs={roles['ip']}"
)

print(
    f"input_runs={roles['input']}"
)

print(
    "DYNAMIC_PROCESSING_PLAN_VALIDATION=PASS"
)
PY

# ------------------------------------------------------------
# Stage 2: dynamic processing cohort -> FASTQ manifest ->
#          download / FastQC / fastp / post-QC / MultiQC
# ------------------------------------------------------------

snakemake \
    results/metadata/validation_fastq_manifest.tsv \
    "${PREP_ARGS[@]}" \
    -p \
    > "$CHIP_SAVE/manifest_setup.log" \
    2>&1

python - <<'PY' \
    > "$CHIP_SAVE/manifest_validation.log"
import csv
from pathlib import Path

processing_path = Path(
    "results/metadata/chipseq/"
    "ena_snapshot_001/analysis_plan/"
    "chipseq_processing_runs.tsv"
)

manifest_path = Path(
    "results/metadata/"
    "validation_fastq_manifest.tsv"
)

with processing_path.open(
    encoding="utf-8",
    newline="",
) as handle:
    processing = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

with manifest_path.open(
    encoding="utf-8",
    newline="",
) as handle:
    manifest = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

processing_runs = {
    row["run_accession"]
    for row in processing
}

manifest_runs = {
    row["run_accession"]
    for row in manifest
}

if processing_runs != manifest_runs:
    raise ValueError(
        "FASTQ manifest run set differs from dynamic processing cohort"
    )

print(
    f"processing_runs={len(processing_runs)}"
)

print(
    f"manifest_runs={len(manifest_runs)}"
)

print(
    f"physical_fastq_files={len(manifest)}"
)

print(
    "DYNAMIC_MANIFEST_VALIDATION=PASS"
)
PY

snakemake \
    chipseq_preprocessing_all \
    "${PREP_ARGS[@]}" \
    -n \
    > "$CHIP_SAVE/preprocessing_dry_run.log" \
    2>&1

/usr/bin/time \
    -v \
    -o "$CHIP_SAVE/workflow_resources.txt" \
    snakemake \
        chipseq_preprocessing_all \
        "${PREP_ARGS[@]}" \
        -p \
        > "$CHIP_SAVE/snakemake_execution.log" \
        2>&1

python - <<'PY' \
    > "$CHIP_SAVE/qc_cohort_validation.log"
import csv
import yaml
from pathlib import Path

processing = Path(
    "results/metadata/chipseq/"
    "ena_snapshot_001/analysis_plan/"
    "chipseq_processing_runs.tsv"
)

with processing.open(
    encoding="utf-8",
    newline="",
) as handle:
    expected = {
        row["run_accession"]
        for row in csv.DictReader(
            handle,
            delimiter="\t",
        )
    }

config = yaml.safe_load(
    Path(
        "config/config.yaml"
    ).read_text()
)

for key in (
    "qc_by_fastq",
    "qc_by_job",
):
    path = Path(
        config[
            "preprocessing"
        ][key]
    )

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        observed = {
            row["run_accession"]
            for row in csv.DictReader(
                handle,
                delimiter="\t",
            )
        }

    if observed != expected:
        raise ValueError(
            f"QC run set mismatch: {path}"
        )

    print(
        f"{path}: {len(observed)} runs"
    )

print(
    "DYNAMIC_QC_COHORT_VALIDATION=PASS"
)
PY

find \
    data \
    results \
    -type f \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "$CHIP_SAVE/output.sha256"

echo "[OK] Dynamic ChIP preprocessing workflow completed."
