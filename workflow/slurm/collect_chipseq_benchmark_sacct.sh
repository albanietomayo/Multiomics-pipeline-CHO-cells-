#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 || ! "$1" =~ ^[0-9]+$ ]]; then
    echo "Usage: $0 JOB_ID OUTPUT_DIRECTORY" >&2
    exit 2
fi

CHIP_JOB_ID="$1"
CHIP_OUTPUT="$2"

if [[ -e "$CHIP_OUTPUT/slurm_accounting.tsv" || -e "$CHIP_OUTPUT/sacct.raw.psv" ]]; then
    echo "ERROR: accounting output collision in $CHIP_OUTPUT" >&2
    exit 2
fi

mkdir -p "$CHIP_OUTPUT"
sacct -j "$CHIP_JOB_ID" \
    --units=K \
    --parsable2 \
    --format=JobID,JobName,State,Elapsed,TotalCPU,AllocCPUS,MaxRSS,AveRSS,ExitCode \
    > "$CHIP_OUTPUT/sacct.raw.psv"

python workflow/scripts/chipseq_benchmark.py parse-sacct \
    --input "$CHIP_OUTPUT/sacct.raw.psv" \
    --output "$CHIP_OUTPUT/slurm_accounting.tsv"
