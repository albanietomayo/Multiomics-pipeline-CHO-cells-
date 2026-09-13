#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --job-name=rnaseq_batch
#SBATCH --output=rnaseq_batch_%j.out

set -euo pipefail

: "${BATCH_ID:?BATCH_ID must be exported by the submission command}"

if [[ ! "$BATCH_ID" =~ ^(SMALL|MEDIUM|LARGE|XLARGE|XXLARGE)_[0-9]{3}$ ]]; then
    echo "ERROR: invalid BATCH_ID: $BATCH_ID" >&2
    exit 2
fi

PROJECT_HOME="${PROJECT_HOME:-$HOME/TFM_multiomics_pipeline_git}"
PLAN_REL="benchmarks/rnaseq/2026-09-10/rnaseq_batches.tsv"
SAVE_ROOT="${RNASEQ_SAVE_ROOT:-$HOME/hpc_rnaseq_production_results}"
JOB_TAG="${SLURM_JOB_ID:-manual}"
WORK="$TMPDIR/TFM_multiomics_pipeline_${BATCH_ID}"
BATCH_SAVE="$SAVE_ROOT/batches/$BATCH_ID/job_$JOB_TAG"
RUNS_SAVE="$SAVE_ROOT/runs"
MAIN_ENV="$TMPDIR/cho-multiomics"
CONDA_ENVS="$TMPDIR/snakemake-conda"
CONDA_PKGS="$TMPDIR/conda-pkgs"
FASTQ_MANIFEST_REL="results/metadata/fastq_manifest.tsv"
VALIDATION_MANIFEST_REL="results/metadata/validation_fastq_manifest.tsv"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
MEM_MB="${MEM_MB:-40000}"
MONITOR_SECONDS="${MONITOR_SECONDS:-30}"
WARMUP_RUN="SRR9621013"
MONITOR_PID=""

mkdir -p "$WORK" "$BATCH_SAVE" "$RUNS_SAVE" "$CONDA_ENVS" "$CONDA_PKGS"

apparent_bytes() {
    local path="$1"
    if [[ -e "$path" ]]; then
        du -sb "$path" 2>/dev/null | awk '{print $1}'
    else
        printf '0\n'
    fi
}

filesystem_used_bytes() {
    df -B1 --output=used "$TMPDIR" | awk 'NR == 2 {gsub(/[[:space:]]/, "", $1); print $1+0}'
}

stop_monitor() {
    if [[ -n "${MONITOR_PID:-}" ]]; then
        kill "$MONITOR_PID" 2>/dev/null || true
        wait "$MONITOR_PID" 2>/dev/null || true
        MONITOR_PID=""
    fi
}

trap stop_monitor EXIT
trap 'stop_monitor; exit 130' INT
trap 'stop_monitor; exit 143' TERM

git -C "$PROJECT_HOME" rev-parse HEAD > "$BATCH_SAVE/pipeline_commit.txt"
git -C "$PROJECT_HOME" status --short > "$BATCH_SAVE/git_status_at_start.txt"

rsync -a \
    --exclude='.git/' \
    --exclude='.snakemake/' \
    --exclude='data/' \
    --exclude='results/' \
    --exclude='resources/' \
    "$PROJECT_HOME/" "$WORK/"

cd "$WORK"

if [[ ! -f "$PLAN_REL" ]]; then
    echo "ERROR: missing batch plan: $PLAN_REL" >&2
    exit 2
fi

mapfile -t RUNS < <(
    awk -F '\t' -v batch="$BATCH_ID" \
        'NR > 1 && $1 == batch {print $2 "\t" $4}' "$PLAN_REL" \
        | sort -n -k1,1 \
        | cut -f2
)

if (( ${#RUNS[@]} == 0 )); then
    echo "ERROR: no runs found for $BATCH_ID" >&2
    exit 2
fi

if (( $(printf '%s\n' "${RUNS[@]}" | sort -u | wc -l) != ${#RUNS[@]} )); then
    echo "ERROR: duplicate runs in $BATCH_ID" >&2
    exit 2
fi

if [[ -n "${PILOT_RUN:-}" ]]; then
    if ! printf '%s\n' "${RUNS[@]}" | grep -Fxq -- "$PILOT_RUN"; then
        echo "ERROR: PILOT_RUN $PILOT_RUN does not belong to $BATCH_ID" >&2
        exit 2
    fi
    RUNS=("$PILOT_RUN")
fi

{
    printf 'run_accession\n'
    printf '%s\n' "$WARMUP_RUN"
    printf '%s\n' "${RUNS[@]}"
} | awk 'NR == 1 || !seen[$0]++' > config/validation_runs.tsv

cp "$PLAN_REL" "$BATCH_SAVE/"
cp config/validation_runs.tsv "$BATCH_SAVE/validation_runs_for_job.tsv"

module purge
module load Miniforge3/24.1.2-0

export CONDA_PKGS_DIRS="$CONDA_PKGS"
export XDG_CACHE_HOME="$TMPDIR/.cache"
mkdir -p "$XDG_CACHE_HOME"

setup_start=$(date +%s)

conda env create \
    --prefix "$MAIN_ENV" \
    --file environment.yml \
    > "$BATCH_SAVE/main_environment_creation.log" 2>&1

export PATH="$MAIN_ENV/bin:$PATH"
SNAKEMAKE="$MAIN_ENV/bin/snakemake"
"$SNAKEMAKE" --version > "$BATCH_SAVE/snakemake_version.txt"

"$SNAKEMAKE" \
    "$FASTQ_MANIFEST_REL" \
    --cores 1 \
    --rerun-incomplete \
    -p \
    > "$BATCH_SAVE/fastq_manifest_setup.log" 2>&1

"$SNAKEMAKE" \
    "$VALIDATION_MANIFEST_REL" \
    --cores 1 \
    --rerun-incomplete \
    -p \
    > "$BATCH_SAVE/validation_manifest_setup.log" 2>&1

cp "$VALIDATION_MANIFEST_REL" "$BATCH_SAVE/validation_fastq_manifest.tsv"

WARMUP_TARGET="results/rnaseq/counts_by_run/$WARMUP_RUN/gene_counts.tsv"
"$SNAKEMAKE" \
    "$WARMUP_TARGET" \
    --use-conda \
    --conda-prefix "$CONDA_ENVS" \
    --cores "$THREADS" \
    --resources "mem_mb=$MEM_MB" \
    --rerun-incomplete \
    --notemp \
    -p \
    > "$BATCH_SAVE/common_setup_warmup.log" 2>&1

rm -rf data/raw/fastq results/qc results/preprocessing results/rnaseq
rm -rf .snakemake/log

"$SNAKEMAKE" \
    star_index \
    --use-conda \
    --conda-prefix "$CONDA_ENVS" \
    --cores "$THREADS" \
    --resources "mem_mb=$MEM_MB" \
    --rerun-incomplete \
    -p \
    > "$BATCH_SAVE/common_reference_star_setup.log" 2>&1

TARGETS=()
for run_accession in "${RUNS[@]}"; do
    TARGETS+=("results/rnaseq/counts_by_run/$run_accession/gene_counts.tsv")
done

"$SNAKEMAKE" \
    -n \
    "${TARGETS[@]}" \
    --use-conda \
    --conda-prefix "$CONDA_ENVS" \
    --cores "$THREADS" \
    --resources "mem_mb=$MEM_MB" \
    --rerun-incomplete \
    -p \
    > "$BATCH_SAVE/preflight_dry_run.txt" 2>&1

setup_end=$(date +%s)
{
    printf 'batch_id=%s\n' "$BATCH_ID"
    printf 'job_id=%s\n' "$JOB_TAG"
    printf 'host=%s\n' "$(hostname)"
    printf 'pipeline_commit=%s\n' "$(cat "$BATCH_SAVE/pipeline_commit.txt")"
    printf 'run_count=%s\n' "${#RUNS[@]}"
    printf 'pilot_run=%s\n' "${PILOT_RUN:-}"
    printf 'execution_cpus=%s\n' "$THREADS"
    printf 'requested_memory_mb=%s\n' "$MEM_MB"
    printf 'setup_elapsed_seconds=%s\n' "$((setup_end - setup_start))"
    printf 'setup_apparent_tmpdir_bytes=%s\n' "$(apparent_bytes "$TMPDIR")"
    printf 'setup_filesystem_used_bytes=%s\n' "$(filesystem_used_bytes)"
} > "$BATCH_SAVE/setup_summary.txt"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    run_accession status elapsed_seconds maximum_rss_kbytes \
    peak_apparent_tmpdir_bytes peak_filesystem_used_bytes \
    gene_counts_bytes attempt_directory \
    > "$BATCH_SAVE/run_metrics.tsv"

FAILED_RUNS=()

copy_essential_outputs() {
    local run_accession="$1"
    local destination="$2"
    local source_dir

    mkdir -p "$destination"

    source_dir="results/preprocessing/fastp/reports"
    if [[ -d "$source_dir" ]]; then
        find "$source_dir" -type f -path "*/$run_accession/*" \
            \( -name '*.json' -o -name '*.html' -o -name '*.log' \) \
            -exec cp --parents -t "$destination" {} +
    fi

    source_dir="results/rnaseq/qc/strandedness/$run_accession"
    if [[ -d "$source_dir" ]]; then
        find "$source_dir" -type f \
            -exec cp --parents -t "$destination" {} +
    fi

    source_dir="results/rnaseq/star/$run_accession"
    if [[ -d "$source_dir" ]]; then
        find "$source_dir" -type f -name 'Log.final.out' \
            -exec cp --parents -t "$destination" {} +
    fi

    source_dir="results/rnaseq/counts/$run_accession"
    if [[ -d "$source_dir" ]]; then
        find "$source_dir" -type f \
            \( -name '*.summary' -o -name '*strandedness.txt' -o -name '*.json' -o -name '*.log' \) \
            -exec cp --parents -t "$destination" {} +
    fi
}

run_is_complete() {
    local run_save="$1"
    local expected_original actual_original

    [[ -f "$run_save/complete.ok" ]] || return 1
    [[ -f "$run_save/gene_counts.tsv.gz" ]] || return 1
    [[ -f "$run_save/gene_counts.tsv.gz.sha256" ]] || return 1
    [[ -f "$run_save/gene_counts.tsv.sha256" ]] || return 1

    (cd "$run_save" && sha256sum --check --status gene_counts.tsv.gz.sha256) \
        || return 1
    gzip -t "$run_save/gene_counts.tsv.gz" || return 1

    expected_original=$(awk '{print $1}' "$run_save/gene_counts.tsv.sha256")
    actual_original=$(gzip -cd "$run_save/gene_counts.tsv.gz" | sha256sum | awk '{print $1}')
    [[ "$expected_original" == "$actual_original" ]]
}

for run_accession in "${RUNS[@]}"; do
    target="results/rnaseq/counts_by_run/$run_accession/gene_counts.tsv"
    RUN_SAVE="$RUNS_SAVE/$run_accession"
    ATTEMPT_SAVE="$RUN_SAVE/attempts/job_$JOB_TAG"
    mkdir -p "$ATTEMPT_SAVE"

    if run_is_complete "$RUN_SAVE"; then
        printf '%s\t%s\t0\t0\t0\t0\t0\t%s\n' \
            "$run_accession" SKIPPED_VALIDATED "$ATTEMPT_SAVE" \
            >> "$BATCH_SAVE/run_metrics.tsv"
        continue
    fi

    rm -rf data/raw/fastq results/qc results/preprocessing results/rnaseq
    rm -rf .snakemake/log
    mkdir -p data/raw/fastq

    scratch_log="$ATTEMPT_SAVE/scratch_usage.tsv"
    printf 'timestamp\tapparent_tmpdir_bytes\tfilesystem_used_bytes\n' > "$scratch_log"

    (
        while true; do
            printf '%s\t%s\t%s\n' \
                "$(date +%s)" \
                "$(apparent_bytes "$TMPDIR")" \
                "$(filesystem_used_bytes)" \
                >> "$scratch_log"
            sleep "$MONITOR_SECONDS"
        done
    ) &
    MONITOR_PID=$!

    start_epoch=$(date +%s)
    set +e
    /usr/bin/time -v \
        -o "$ATTEMPT_SAVE/time_verbose.txt" \
        "$SNAKEMAKE" \
        "$target" \
        --use-conda \
        --conda-prefix "$CONDA_ENVS" \
        --cores "$THREADS" \
        --resources "mem_mb=$MEM_MB" \
        --rerun-incomplete \
        --notemp \
        -p \
        > "$ATTEMPT_SAVE/snakemake_execution.log" 2>&1
    workflow_status=$?
    set -e
    end_epoch=$(date +%s)

    stop_monitor
    printf '%s\t%s\t%s\n' \
        "$(date +%s)" \
        "$(apparent_bytes "$TMPDIR")" \
        "$(filesystem_used_bytes)" \
        >> "$scratch_log"

    elapsed_seconds=$((end_epoch - start_epoch))
    max_rss_kbytes=$(awk -F ': ' \
        '/Maximum resident set size/ {gsub(/[[:space:]]/, "", $2); print $2}' \
        "$ATTEMPT_SAVE/time_verbose.txt")
    max_rss_kbytes=${max_rss_kbytes:-0}
    peak_apparent=$(awk -F '\t' 'NR > 1 && $2 > max {max=$2} END {print max+0}' "$scratch_log")
    peak_filesystem=$(awk -F '\t' 'NR > 1 && $3 > max {max=$3} END {print max+0}' "$scratch_log")
    gene_counts_bytes=$(apparent_bytes "$target")

    find data results -type f -printf '%p\t%s\n' 2>/dev/null \
        | sort > "$ATTEMPT_SAVE/output_inventory.tsv"

    {
        printf 'run_accession=%s\n' "$run_accession"
        printf 'batch_id=%s\n' "$BATCH_ID"
        printf 'job_id=%s\n' "$JOB_TAG"
        printf 'host=%s\n' "$(hostname)"
        printf 'pipeline_commit=%s\n' "$(cat "$BATCH_SAVE/pipeline_commit.txt")"
        printf 'workflow_exit_status=%s\n' "$workflow_status"
        printf 'elapsed_seconds=%s\n' "$elapsed_seconds"
        printf 'maximum_rss_kbytes=%s\n' "$max_rss_kbytes"
        printf 'peak_apparent_tmpdir_bytes=%s\n' "$peak_apparent"
        printf 'peak_filesystem_used_bytes=%s\n' "$peak_filesystem"
        printf 'gene_counts_bytes=%s\n' "$gene_counts_bytes"
        printf 'end_time=%s\n' "$(date --iso-8601=seconds)"
    } > "$ATTEMPT_SAVE/resource_summary.txt"

    if (( workflow_status != 0 )); then
        copy_essential_outputs "$run_accession" "$ATTEMPT_SAVE/structured_outputs"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$run_accession" FAILED "$elapsed_seconds" "$max_rss_kbytes" \
            "$peak_apparent" "$peak_filesystem" "$gene_counts_bytes" \
            "$ATTEMPT_SAVE" >> "$BATCH_SAVE/run_metrics.tsv"
        echo "ERROR: $run_accession failed with status $workflow_status" >&2
        FAILED_RUNS+=("$run_accession")
        rm -rf data/raw/fastq results/qc results/preprocessing results/rnaseq
        rm -rf .snakemake/log
        continue
    fi

    if [[ ! -s "$target" ]]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$run_accession" FAILED_MISSING_TARGET "$elapsed_seconds" "$max_rss_kbytes" \
            "$peak_apparent" "$peak_filesystem" "$gene_counts_bytes" \
            "$ATTEMPT_SAVE" >> "$BATCH_SAVE/run_metrics.tsv"
        echo "ERROR: missing or empty target after successful workflow: $target" >&2
        FAILED_RUNS+=("$run_accession")
        rm -rf data/raw/fastq results/qc results/preprocessing results/rnaseq
        rm -rf .snakemake/log
        continue
    fi

    PERSIST_STAGE="$RUN_SAVE/.persist_job_$JOB_TAG"
    mkdir -p "$PERSIST_STAGE/structured_outputs"
    copy_essential_outputs "$run_accession" "$PERSIST_STAGE/structured_outputs"

    original_hash=$(sha256sum "$target" | awk '{print $1}')
    gzip -c "$target" > "$PERSIST_STAGE/gene_counts.tsv.gz"
    gzip -t "$PERSIST_STAGE/gene_counts.tsv.gz"
    compressed_hash=$(sha256sum "$PERSIST_STAGE/gene_counts.tsv.gz" | awk '{print $1}')

    printf '%s  gene_counts.tsv\n' "$original_hash" \
        > "$PERSIST_STAGE/gene_counts.tsv.sha256"
    printf '%s  gene_counts.tsv.gz\n' "$compressed_hash" \
        > "$PERSIST_STAGE/gene_counts.tsv.gz.sha256"

    mkdir -p "$RUN_SAVE/structured_outputs"
    cp -a "$PERSIST_STAGE/structured_outputs/." "$RUN_SAVE/structured_outputs/"
    mv "$PERSIST_STAGE/gene_counts.tsv.gz" "$RUN_SAVE/gene_counts.tsv.gz"
    mv "$PERSIST_STAGE/gene_counts.tsv.sha256" "$RUN_SAVE/gene_counts.tsv.sha256"
    mv "$PERSIST_STAGE/gene_counts.tsv.gz.sha256" "$RUN_SAVE/gene_counts.tsv.gz.sha256"

    (cd "$RUN_SAVE" && sha256sum --check --status gene_counts.tsv.gz.sha256)
    actual_original=$(gzip -cd "$RUN_SAVE/gene_counts.tsv.gz" | sha256sum | awk '{print $1}')
    [[ "$actual_original" == "$original_hash" ]]
    printf 'completed_at=%s\njob_id=%s\nbatch_id=%s\n' \
        "$(date --iso-8601=seconds)" "$JOB_TAG" "$BATCH_ID" \
        > "$RUN_SAVE/complete.ok"
    rm -rf "$PERSIST_STAGE"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$run_accession" COMPLETED "$elapsed_seconds" "$max_rss_kbytes" \
        "$peak_apparent" "$peak_filesystem" "$gene_counts_bytes" \
        "$ATTEMPT_SAVE" >> "$BATCH_SAVE/run_metrics.tsv"

    rm -rf data/raw/fastq results/qc results/preprocessing results/rnaseq
    rm -rf .snakemake/log
done

if (( ${#FAILED_RUNS[@]} > 0 )); then
    {
        printf 'finished_at=%s\n' "$(date --iso-8601=seconds)"
        printf 'run_count=%s\n' "${#RUNS[@]}"
        printf 'failed_run_count=%s\n' "${#FAILED_RUNS[@]}"
        printf 'failed_run_accession=%s\n' "${FAILED_RUNS[@]}"
    } > "$BATCH_SAVE/finished_with_failures.txt"

    echo "Batch finished with ${#FAILED_RUNS[@]} failed run(s): ${FAILED_RUNS[*]}" >&2
    echo "Persistent batch logs: $BATCH_SAVE"
    echo "Persistent run outputs: $RUNS_SAVE"
    exit 1
fi

printf 'completed_at=%s\nrun_count=%s\n' \
    "$(date --iso-8601=seconds)" "${#RUNS[@]}" \
    > "$BATCH_SAVE/complete.ok"

echo "Batch completed successfully: $BATCH_ID"
echo "Persistent batch logs: $BATCH_SAVE"
echo "Persistent run outputs: $RUNS_SAVE"
