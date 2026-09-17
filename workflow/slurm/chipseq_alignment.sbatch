#!/bin/bash
#SBATCH --job-name=chipseq_alignment
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00

set -euo pipefail

: "${TMPDIR:?Node-local TMPDIR required}"

CHIP_SOURCE="$PWD"
CHIP_SAVE="$(dirname "$CHIP_SOURCE")/job_${SLURM_JOB_ID}"
CHIP_WORK="$TMPDIR/chipseq-alignment"
CHIP_STAGE=source_validation

mkdir -p \
    "$CHIP_SAVE" \
    "$CHIP_WORK"

finish() {
    CHIP_EXIT=$?

    trap - EXIT
    set +e

    if [[ -d "$CHIP_WORK/results/chipseq/alignment" && "$CHIP_EXIT" != 0 ]]; then

        rsync -a \
            --exclude='*.bam' \
            --exclude='*.csi' \
            --exclude='*.gz' \
            "$CHIP_WORK/results/chipseq/alignment/" \
            "$CHIP_SAVE/failure_reports/" \
            > "$CHIP_SAVE/failure_copy.log" \
            2>&1
    fi

    if [[ -d "$CHIP_WORK/.snakemake/log" ]]; then
        rsync -a \
            "$CHIP_WORK/.snakemake/log/" \
            "$CHIP_SAVE/snakemake_logs/"
    fi

    printf \
        'stage\t%s\nexit_status\t%s\n' \
        "$CHIP_STAGE" \
        "$CHIP_EXIT" \
        > "$CHIP_SAVE/job_status.tsv"

    echo \
        "Alignment exit=$CHIP_EXIT; stage=$CHIP_STAGE; logs=$CHIP_SAVE"

    exit "$CHIP_EXIT"
}

trap finish EXIT

sha256sum \
    --check ../source.sha256 \
    > "$CHIP_SAVE/source_validation.log"

rsync -a \
    "$CHIP_SOURCE/" \
    "$CHIP_WORK/"

cd "$CHIP_WORK"

module purge
module load Miniforge3/24.1.2-0

export CONDA_PKGS_DIRS="$TMPDIR/chipseq-alignment-conda-pkgs"
export XDG_CACHE_HOME="$TMPDIR/chipseq-alignment-cache"
export PYTHONDONTWRITEBYTECODE=1

mkdir -p \
    "$CONDA_PKGS_DIRS" \
    "$XDG_CACHE_HOME"

CHIP_ENV="$TMPDIR/chipseq-alignment-env"
CHIP_RULE_ENVS="$TMPDIR/chipseq-alignment-rule-envs"

CHIP_STAGE=environment_creation

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

conda list \
    --prefix "$CHIP_ENV" \
    --explicit \
    > "$CHIP_SAVE/environment_explicit.txt"

CHIP_ARGS=(
    --snakefile workflow/rules/chipseq_alignment.smk
    --cores "$SLURM_CPUS_PER_TASK"
    --use-conda
    --conda-prefix "$CHIP_RULE_ENVS"
    --rerun-incomplete
)

CHIP_STAGE=alignment_plan_validation

python - <<'PY' \
    > "$CHIP_SAVE/alignment_plan_validation.log"
import json
from collections import Counter
from pathlib import Path

path = Path(
    "config/chipseq_alignment_inputs.json"
)

obj = json.loads(
    path.read_text(
        encoding="utf-8"
    )
)

if obj.get("schema_version") != 1:
    raise ValueError(
        "Unsupported alignment-plan schema"
    )

runs = obj.get(
    "runs"
)

if not isinstance(
    runs,
    list,
) or not runs:
    raise ValueError(
        "Alignment plan contains no runs"
    )

accessions = [
    row["run_accession"]
    for row in runs
]

if len(accessions) != len(
    set(accessions)
):
    raise ValueError(
        "Duplicate run accessions in alignment plan"
    )

roles = Counter(
    row["role"]
    for row in runs
)

if set(roles) - {
    "ip",
    "input",
}:
    raise ValueError(
        f"Unexpected roles: {dict(roles)}"
    )

if roles["ip"] < 1:
    raise ValueError(
        "Alignment plan contains no IP runs"
    )

if roles["input"] < 1:
    raise ValueError(
        "Alignment plan contains no Input runs"
    )

if obj.get(
    "library_layout"
) != "SINGLE":
    raise ValueError(
        "Operational alignment supports SINGLE only"
    )

if obj.get(
    "instrument_platform"
) != "ILLUMINA":
    raise ValueError(
        "Operational alignment supports ILLUMINA only"
    )

print(
    f"alignment_runs={len(runs)}"
)

print(
    f"ip_runs={roles['ip']}"
)

print(
    f"input_runs={roles['input']}"
)

print(
    "DYNAMIC_ALIGNMENT_PLAN_VALIDATION=PASS"
)
PY

CHIP_STAGE=verified_fastq_copy

python \
    workflow/scripts/chipseq_alignment_support.py \
    stage \
    > "$CHIP_SAVE/fastq_copy_validation.log" \
    2>&1

CHIP_STAGE=dry_run

snakemake \
    chipseq_alignment_all \
    "${CHIP_ARGS[@]}" \
    -n \
    -p \
    > "$CHIP_SAVE/dry_run.log" \
    2>&1

CHIP_STAGE=alignment

/usr/bin/time \
    -v \
    -o "$CHIP_SAVE/workflow_resources.txt" \
    snakemake \
        chipseq_alignment_all \
        "${CHIP_ARGS[@]}" \
        -p \
        > "$CHIP_SAVE/snakemake_execution.log" \
        2>&1

CHIP_STAGE=environment_provenance

mkdir -p \
    "$CHIP_SAVE/rule_environments"

for CHIP_RULE_ENV in "$CHIP_RULE_ENVS"/*
do
    if [[ -d "$CHIP_RULE_ENV/conda-meta" ]]; then

        conda list \
            --prefix "$CHIP_RULE_ENV" \
            --explicit \
            > "$CHIP_SAVE/rule_environments/$(basename "$CHIP_RULE_ENV").txt"
    fi
done

CHIP_STAGE=output_copy

python \
    workflow/scripts/chipseq_alignment_support.py \
    publish \
    --save "$CHIP_SAVE" \
    > "$CHIP_SAVE/output_copy_validation.log" \
    2>&1

CHIP_STAGE=completed

echo \
    "[OK] Dynamic ChIP alignment and QC completed; outputs copied and verified."
