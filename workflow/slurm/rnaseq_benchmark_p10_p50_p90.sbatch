#!/usr/bin/env bash
#SBATCH -A C3SE2026-1-42
#SBATCH -p vera
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --job-name=rnaseq_bench_3pct
#SBATCH --output=rnaseq_benchmark_%j.out

set -euo pipefail

PROJECT_HOME="$HOME/TFM_multiomics_pipeline_git"
WORK="$TMPDIR/TFM_multiomics_pipeline_benchmark"
SAVE="$HOME/hpc_rnaseq_benchmark_results/$SLURM_JOB_ID"
MAIN_ENV="$TMPDIR/cho-multiomics"
CONDA_ENVS="$TMPDIR/snakemake-conda"
CONDA_PKGS="$TMPDIR/conda-pkgs"
REFERENCE_REL="resources/reference/CriGri-PICRH-1.0"
FASTQ_MANIFEST_REL="results/metadata/fastq_manifest.tsv"
VALIDATION_MANIFEST_REL="results/metadata/validation_fastq_manifest.tsv"
THREADS=8
MEM_MB=90000
MONITOR_SECONDS=15

mkdir -p "$WORK" "$SAVE" "$CONDA_ENVS" "$CONDA_PKGS"

bytes_of() {
    local path="$1"
    if [[ -e "$path" ]]; then
        du -sb "$path" 2>/dev/null | awk '{print $1}'
    else
        printf '0\n'
    fi
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

SETUP_START=$(date +%s)

git -C "$PROJECT_HOME" rev-parse HEAD > "$SAVE/pipeline_commit.txt"
git -C "$PROJECT_HOME" status --short > "$SAVE/git_status_at_start.txt"

rsync -a \
    --exclude='.git/' \
    --exclude='.snakemake/' \
    --exclude='data/' \
    --exclude='results/' \
    --exclude='resources/' \
    "$PROJECT_HOME/" "$WORK/"

mkdir -p "$WORK/results/metadata" "$WORK/resources/reference"

cd "$WORK"

cp config/rnaseq_hpc_benchmark_runs.tsv "$SAVE/"
cp config/validation_runs.tsv "$SAVE/"

module purge
module load Miniforge3/24.1.2-0

export CONDA_PKGS_DIRS="$CONDA_PKGS"
export XDG_CACHE_HOME="$TMPDIR/.cache"
mkdir -p "$XDG_CACHE_HOME"

conda env create \
    --prefix "$MAIN_ENV" \
    --file environment.yml \
    > "$SAVE/main_environment_creation.log" 2>&1

export PATH="$MAIN_ENV/bin:$PATH"
SNAKEMAKE="$MAIN_ENV/bin/snakemake"
"$SNAKEMAKE" --version > "$SAVE/snakemake_version.txt"

"$SNAKEMAKE" \
    "$FASTQ_MANIFEST_REL" \
    --cores 1 \
    --rerun-incomplete \
    -p \
    > "$SAVE/fastq_manifest_setup.log" 2>&1

declare -A RUN_FOR_CLASS
declare -A FASTQ_BYTES_FOR_CLASS
declare -A FASTQ_GB_FOR_CLASS

while IFS=$'\t' read -r benchmark_class run_accession study_accession omics fastq_structure processing_layout n_fastq fastq_bytes fastq_gb; do
    [[ "$benchmark_class" == "benchmark_class" ]] && continue
    case "$benchmark_class" in
        P10|P50|P90)
            RUN_FOR_CLASS["$benchmark_class"]="$run_accession"
            FASTQ_BYTES_FOR_CLASS["$benchmark_class"]="$fastq_bytes"
            FASTQ_GB_FOR_CLASS["$benchmark_class"]="$fastq_gb"
            ;;
    esac
done < config/rnaseq_hpc_benchmark_runs.tsv

CLASSES=(P10 P50 P90)
TARGETS=()

for benchmark_class in "${CLASSES[@]}"; do
    if [[ -z "${RUN_FOR_CLASS[$benchmark_class]:-}" ]]; then
        echo "ERROR: missing $benchmark_class in config/rnaseq_hpc_benchmark_runs.tsv" >&2
        exit 1
    fi
    run_accession="${RUN_FOR_CLASS[$benchmark_class]}"
    TARGETS+=("results/rnaseq/counts_by_run/$run_accession/gene_counts.tsv")
done

"$SNAKEMAKE" \
    "$VALIDATION_MANIFEST_REL" \
    --cores 1 \
    --rerun-incomplete \
    -p \
    > "$SAVE/validation_manifest_setup.log" 2>&1

WARMUP_RUN="${RUN_FOR_CLASS[P10]}"
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
    > "$SAVE/common_setup_warmup_P10.log" 2>&1

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
    > "$SAVE/common_reference_star_setup.log" 2>&1

"$SNAKEMAKE" \
    -n \
    "${TARGETS[@]}" \
    --use-conda \
    --conda-prefix "$CONDA_ENVS" \
    --cores "$THREADS" \
    --resources "mem_mb=$MEM_MB" \
    --rerun-incomplete \
    -p \
    > "$SAVE/preflight_dry_run.txt" 2>&1

SETUP_END=$(date +%s)
SETUP_SECONDS=$((SETUP_END - SETUP_START))
COMMON_BASELINE_BYTES=$(bytes_of "$TMPDIR")

{
    printf 'setup_elapsed_seconds=%s\n' "$SETUP_SECONDS"
    printf 'common_tmpdir_bytes=%s\n' "$COMMON_BASELINE_BYTES"
    printf 'reference_bytes=%s\n' "$(bytes_of "$WORK/$REFERENCE_REL")"
    printf 'star_index_bytes=%s\n' "$(bytes_of "$WORK/$REFERENCE_REL/star_index")"
    printf 'main_environment_bytes=%s\n' "$(bytes_of "$MAIN_ENV")"
    printf 'rule_environments_bytes=%s\n' "$(bytes_of "$CONDA_ENVS")"
    printf 'conda_package_cache_bytes=%s\n' "$(bytes_of "$CONDA_PKGS")"
} > "$SAVE/common_setup_summary.txt"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    benchmark_class run_accession fastq_bytes fastq_GB baseline_tmpdir_bytes \
    sampled_peak_tmpdir_bytes incremental_peak_tmpdir_bytes elapsed_seconds \
    max_rss_kbytes target_bytes raw_fastq_bytes preprocessing_bytes \
    rnaseq_bytes qc_bytes workflow_exit_status cleanup_tmpdir_bytes \
    > "$SAVE/benchmark_metrics.tsv"

for benchmark_class in "${CLASSES[@]}"; do
    run_accession="${RUN_FOR_CLASS[$benchmark_class]}"
    fastq_bytes="${FASTQ_BYTES_FOR_CLASS[$benchmark_class]}"
    fastq_gb="${FASTQ_GB_FOR_CLASS[$benchmark_class]}"
    target="results/rnaseq/counts_by_run/$run_accession/gene_counts.tsv"
    RUN_SAVE="$SAVE/${benchmark_class}_${run_accession}"

    mkdir -p "$RUN_SAVE"

    rm -rf data/raw/fastq results/qc results/preprocessing results/rnaseq
    rm -rf .snakemake/log
    mkdir -p data/raw/fastq

    baseline_bytes=$(bytes_of "$TMPDIR")
    scratch_log="$RUN_SAVE/scratch_usage.tsv"
    printf 'timestamp\ttmpdir_bytes\n' > "$scratch_log"

    (
        while true; do
            printf '%s\t%s\n' \
                "$(date +%s)" \
                "$(bytes_of "$TMPDIR")" \
                >> "$scratch_log"
            sleep "$MONITOR_SECONDS"
        done
    ) &
    MONITOR_PID=$!

    start_epoch=$(date +%s)

    set +e
    /usr/bin/time -v \
        -o "$RUN_SAVE/time_verbose.txt" \
        "$SNAKEMAKE" \
        "$target" \
        --use-conda \
        --conda-prefix "$CONDA_ENVS" \
        --cores "$THREADS" \
        --resources "mem_mb=$MEM_MB" \
        --rerun-incomplete \
        --notemp \
        -p \
        2>&1 | tee "$RUN_SAVE/snakemake_execution.log"
    workflow_status=${PIPESTATUS[0]}
    set -e

    end_epoch=$(date +%s)
    elapsed_seconds=$((end_epoch - start_epoch))

    stop_monitor
    printf '%s\t%s\n' \
        "$(date +%s)" \
        "$(bytes_of "$TMPDIR")" \
        >> "$scratch_log"

    sampled_peak_bytes=$(awk -F'\t' 'NR > 1 && $2 > max {max=$2} END {print max+0}' "$scratch_log")
    if (( sampled_peak_bytes > baseline_bytes )); then
        incremental_peak_bytes=$((sampled_peak_bytes - baseline_bytes))
    else
        incremental_peak_bytes=0
    fi

    max_rss_kbytes=$(awk -F': ' \
        '/Maximum resident set size/ {gsub(/[[:space:]]/, "", $2); print $2}' \
        "$RUN_SAVE/time_verbose.txt")
    max_rss_kbytes=${max_rss_kbytes:-0}

    target_bytes=$(bytes_of "$target")
    raw_fastq_bytes=$(bytes_of data/raw/fastq)
    preprocessing_bytes=$(bytes_of results/preprocessing)
    rnaseq_bytes=$(bytes_of results/rnaseq)
    qc_bytes=$(bytes_of results/qc)

    find data results \
        -type f \
        -printf '%p\t%s\n' \
        2>/dev/null \
        | sort > "$RUN_SAVE/output_inventory.tsv"

    mkdir -p "$RUN_SAVE/small_outputs"
    find data results \
        -type f \
        -size -20M \
        -print0 \
        2>/dev/null \
        | xargs -0 -r cp --parents -t "$RUN_SAVE/small_outputs"

    if [[ -d .snakemake/log ]]; then
        cp -a .snakemake/log "$RUN_SAVE/snakemake_logs"
    fi

    {
        printf 'benchmark_class=%s\n' "$benchmark_class"
        printf 'run_accession=%s\n' "$run_accession"
        printf 'job_id=%s\n' "$SLURM_JOB_ID"
        printf 'host=%s\n' "$(hostname)"
        printf 'pipeline_commit=%s\n' "$(cat "$SAVE/pipeline_commit.txt")"
        printf 'target=%s\n' "$target"
        printf 'execution_cpus=%s\n' "$THREADS"
        printf 'requested_memory_mb=%s\n' "$MEM_MB"
        printf 'workflow_exit_status=%s\n' "$workflow_status"
        printf 'elapsed_seconds=%s\n' "$elapsed_seconds"
        printf 'maximum_rss_kbytes=%s\n' "$max_rss_kbytes"
        printf 'baseline_tmpdir_bytes=%s\n' "$baseline_bytes"
        printf 'sampled_peak_tmpdir_bytes=%s\n' "$sampled_peak_bytes"
        printf 'incremental_peak_tmpdir_bytes=%s\n' "$incremental_peak_bytes"
        printf 'target_bytes=%s\n' "$target_bytes"
        printf 'raw_fastq_bytes=%s\n' "$raw_fastq_bytes"
        printf 'preprocessing_bytes=%s\n' "$preprocessing_bytes"
        printf 'rnaseq_bytes=%s\n' "$rnaseq_bytes"
        printf 'qc_bytes=%s\n' "$qc_bytes"
        printf 'end_time=%s\n' "$(date)"
    } > "$RUN_SAVE/resource_summary.txt"

    rm -rf data/raw/fastq results/qc results/preprocessing results/rnaseq
    rm -rf .snakemake/log
    cleanup_tmpdir_bytes=$(bytes_of "$TMPDIR")

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$benchmark_class" "$run_accession" "$fastq_bytes" "$fastq_gb" \
        "$baseline_bytes" "$sampled_peak_bytes" "$incremental_peak_bytes" \
        "$elapsed_seconds" "$max_rss_kbytes" "$target_bytes" \
        "$raw_fastq_bytes" "$preprocessing_bytes" "$rnaseq_bytes" \
        "$qc_bytes" "$workflow_status" "$cleanup_tmpdir_bytes" \
        >> "$SAVE/benchmark_metrics.tsv"

    if (( workflow_status != 0 )); then
        echo "ERROR: $benchmark_class $run_accession failed with status $workflow_status" >&2
        exit "$workflow_status"
    fi
done

echo "Benchmark completed successfully."
echo "Persistent results: $SAVE"
