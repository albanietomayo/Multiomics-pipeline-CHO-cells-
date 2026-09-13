#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --job-name=rnaseq_batch
#SBATCH --output=rnaseq_batch_%j.out

set -euo pipefail

BATCH_ID="SMALL_007"
PILOT_RUN="DRR091439"
PROJECT_HOME="/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_git"
RNASEQ_SAVE_ROOT="/cephyr/users/mayoa/Vera/hpc_rnaseq_diagnostics/DRR091439"

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


DIAG="$BATCH_SAVE/diagnostic"
mkdir -p "$DIAG"

BAM="results/rnaseq/star/DRR091439/SINGLE/Aligned.sortedByCoord.out.bam"
STRAND="results/rnaseq/qc/strandedness/DRR091439/resolved_strandedness.tsv"
GTF="resources/reference/CriGri-PICRH-1.0/annotation.gtf"

"$SNAKEMAKE" "$BAM" "$STRAND" \
    --use-conda --conda-prefix "$CONDA_ENVS" \
    --cores "$THREADS" --resources "mem_mb=$MEM_MB" \
    --rerun-incomplete --notemp -p \
    > "$DIAG/prepare_bam.log" 2>&1

# Conservar los datos antes de probar featureCounts.
cp "$BAM" "$DIAG/input.bam"
cp "$STRAND" "$DIAG/resolved_strandedness.tsv"
cp "$GTF" "$DIAG/annotation.gtf"
cp workflow/scripts/run_featurecounts.py "$DIAG/"
sha256sum "$DIAG/input.bam" "$DIAG/annotation.gtf" \
    > "$DIAG/input_checksums.txt"

mapfile -t FC_BINS < <(
    find "$CONDA_ENVS" -type f -path '*/bin/featureCounts'
)
if (( ${#FC_BINS[@]} != 1 )); then
    printf 'Expected one featureCounts binary, found %s\n' \
        "${#FC_BINS[@]}" >&2
    exit 2
fi

RNA_ENV="${FC_BINS[0]%/bin/featureCounts}"
export PATH="$RNA_ENV/bin:$PATH"
featureCounts -v > "$DIAG/featurecounts_version.txt" 2>&1
samtools --version > "$DIAG/samtools_version.txt"
conda list --prefix "$RNA_ENV" --explicit > "$DIAG/conda_explicit.txt"

# Registrar por separado la comprobación rápida y la lectura completa.
set +e
samtools quickcheck -v "$DIAG/input.bam" \
    > "$DIAG/bam_quickcheck.log" 2>&1
quick_status=$?
samtools view -c "$DIAG/input.bam" \
    > "$DIAG/bam_alignment_count.txt" 2> "$DIAG/bam_read_errors.log"
read_status=$?
set -e

printf 'quickcheck_exit=%s\nfull_read_exit=%s\n' \
    "$quick_status" "$read_status" > "$DIAG/bam_validation.txt"

if (( quick_status != 0 || read_status != 0 )); then
    echo "BAM validation failed. Diagnostic files: $DIAG"
    exit 2
fi

printf 'threads\texit_status\n' > "$DIAG/featurecounts_tests.tsv"
for nthreads in 1 4; do
    set +e
    /usr/bin/time -v -o "$DIAG/time_T${nthreads}.txt" \
        "$RNA_ENV/bin/python" "$DIAG/run_featurecounts.py" \
        "$DIAG/resolved_strandedness.tsv" DRR091439 SINGLE \
        "$DIAG/annotation.gtf" "$DIAG/input.bam" "$nthreads" \
        "$DIAG/counts_T${nthreads}.txt" \
        "$DIAG/strandedness_T${nthreads}.txt" \
        > "$DIAG/featurecounts_T${nthreads}.log" 2>&1
    test_status=$?
    set -e
    printf '%s\t%s\n' "$nthreads" "$test_status" \
        >> "$DIAG/featurecounts_tests.tsv"
done

echo "Diagnostic tests finished. Inspect exit statuses in: $DIAG"
cat "$DIAG/bam_validation.txt" "$DIAG/featurecounts_tests.tsv"
