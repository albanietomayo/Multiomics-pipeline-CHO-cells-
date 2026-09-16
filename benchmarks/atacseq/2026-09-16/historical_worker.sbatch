#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --job-name=atac_bulk_CHO_TSS

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
LAUNCH="${ATAC_LAUNCH_DIR:?ATAC_LAUNCH_DIR is required}"
JOB="${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
SCRATCH="${TMPDIR:?TMPDIR is required}"
WORK="$SCRATCH/atac_bulk_${JOB}"
SAVE="$LAUNCH/job_${JOB}"
THREADS=16
MEM_MB=90000
[[ ! -e "$WORK" && ! -e "$SAVE" ]] || { echo 'Work/save path already exists' >&2; exit 2; }
mkdir -p "$WORK" "$SAVE"

QUOTA_ROOT=/cephyr/users/mayoa
check_capacity() {
    local required="$1" stats used entries
    stats=$(cat "$QUOTA_ROOT") || return 1
    used=$(awk '$1=="rbytes:" {print $2}' <<< "$stats")
    entries=$(awk '$1=="rentries:" {print $2}' <<< "$stats")
    [[ "$used" =~ ^[0-9]+$ && "$entries" =~ ^[0-9]+$ ]] || return 1
    printf 'home_bytes=%s required_bytes=%s home_entries=%s\n' "$used" "$required" "$entries"
    (( used + required + 3*1024*1024*1024 <= 30*1024*1024*1024 && entries+500 < 60000 ))
}

finish() {
    local result=$? diagnostics_status=0
    trap - EXIT
    set +e
    local paths=() item
    for item in .snakemake/log results/metadata results/atacseq/qc results/preprocessing/fastp/reports logs; do
        [[ ! -e "$WORK/$item" ]] || paths+=("$item")
    done
    if ((${#paths[@]})); then
        tar -C "$WORK" -czf "$SAVE/diagnostics.tar.gz" -- "${paths[@]}"
        diagnostics_status=$?
    fi
    if ((result==0 && diagnostics_status!=0)); then result=$diagnostics_status; fi
    if ((result==0)) && [[ ! -f "$SAVE/all_runs_verified.ok" ]]; then result=9; fi
    printf 'exit_status\t%s\nfinished_at\t%s\n' "$result" "$(date --iso-8601=seconds)" > "$SAVE/job_status.tsv"
    if ((result==0)); then
        printf 'Both eligible single-end runs completed; BAM and TSS outputs verified.\n' > "$SAVE/complete.ok" || result=9
    fi
    printf 'ATAC job exit=%s; results=%s\n' "$result" "$SAVE"
    exit "$result"
}
trap finish EXIT
trap 'exit 143' TERM
trap 'exit 130' INT

check_capacity 0 > "$SAVE/home_capacity_at_start.txt" || { echo 'Insufficient HOME margin' >&2; exit 5; }
scontrol show job "$JOB" > "$SAVE/slurm_job_at_start.txt"
df -h "$SCRATCH" > "$SAVE/scratch_at_start.txt"
cp "$LAUNCH/worker.sbatch" "$SAVE/submitted_worker.sbatch"
(cd "$LAUNCH" && sha256sum --check snapshot.sha256) > "$SAVE/snapshot_verification.txt"
tar -xzf "$LAUNCH/snapshot.tar.gz" -C "$WORK"
cd "$WORK"

module purge
module load Miniforge3/24.1.2-0
export CONDA_PKGS_DIRS="$SCRATCH/atac-conda-pkgs"
export XDG_CACHE_HOME="$SCRATCH/atac-cache"
mkdir -p "$CONDA_PKGS_DIRS" "$XDG_CACHE_HOME"
MAIN_ENV="$SCRATCH/atac-main-env"
CONDA_ENVS="$SCRATCH/atac-rule-envs"
conda env create --prefix "$MAIN_ENV" --file environment.yml > "$SAVE/main_environment_creation.log" 2>&1
export PATH="$MAIN_ENV/bin:$PATH"
SNAKEMAKE="$MAIN_ENV/bin/snakemake"
"$SNAKEMAKE" --version > "$SAVE/snakemake_version.txt"

python workflow/scripts/selection_gate.py --repo "$WORK" assert SRR12774931 SRR12774932 --omics ATAC-seq
python workflow/scripts/build_fastq_manifest.py --samples config/samples.tsv --output results/metadata/fastq_manifest.tsv > "$SAVE/build_manifest.log" 2>&1
python workflow/scripts/build_validation_manifest.py --manifest results/metadata/fastq_manifest.tsv --validation-runs config/validation_runs.tsv --output results/metadata/validation_fastq_manifest.tsv >> "$SAVE/build_manifest.log" 2>&1
python atac_job_checks.py manifest > "$SAVE/selected_inputs.tsv"
cp config/validation_runs.tsv "$SAVE/validation_runs_for_job.tsv"

# Same files, same rules, only explicit ATAC targets. One run at a time.
for RUN in SRR12774931 SRR12774932; do
    RUN_SAVE="$SAVE/runs/$RUN"
    mkdir -p "$RUN_SAVE"
    BAM="results/atacseq/filtered_bam/$RUN/filtered.bam"
    QC="results/atacseq/qc/$RUN"
    TARGETS=("$BAM" "${BAM}.bai" "$QC/tss_profile.tsv" "$QC/tss_enrichment_summary.tsv")
    python workflow/scripts/selection_gate.py --repo "$WORK" assert "$RUN" --omics ATAC-seq
    printf 'START %s %s\n' "$RUN" "$(date --iso-8601=seconds)"
    SM_ARGS=(--snakefile Snakefile.atac_bulk --use-conda --conda-prefix "$CONDA_ENVS" --cores "$THREADS" --resources "mem_mb=$MEM_MB" --default-resources "tmpdir=$SCRATCH" --rerun-incomplete --latency-wait 60 -p)
    "$SNAKEMAKE" "${TARGETS[@]}" "${SM_ARGS[@]}" -n > "$RUN_SAVE/dry_run.log" 2>&1
    set +e
    /usr/bin/time -v -o "$RUN_SAVE/time_verbose.txt" "$SNAKEMAKE" "${TARGETS[@]}" "${SM_ARGS[@]}" > "$RUN_SAVE/snakemake_execution.log" 2>&1
    WORKFLOW_STATUS=$?
    set -e
    printf 'workflow_exit_status\t%s\n' "$WORKFLOW_STATUS" > "$RUN_SAVE/workflow_status.tsv"
    ((WORKFLOW_STATUS==0)) || exit "$WORKFLOW_STATUS"

    BOWTIE2=$(find "$CONDA_ENVS" -path '*/bin/bowtie2' -print -quit)
    [[ -n "$BOWTIE2" ]] || { echo 'Bowtie2 environment not found' >&2; exit 3; }
    SAMTOOLS="$(dirname "$BOWTIE2")/samtools"
    [[ -x "$SAMTOOLS" ]] || exit 3
    "$SAMTOOLS" --version > "$RUN_SAVE/samtools_version.txt"
    "$BOWTIE2" --version > "$RUN_SAVE/bowtie2_version.txt"
    "$SAMTOOLS" quickcheck -v "$BAM"
    "$SAMTOOLS" flagstat "$BAM" > "$RUN_SAVE/filtered.flagstat.txt"
    "$SAMTOOLS" idxstats "$BAM" > "$RUN_SAVE/filtered.idxstats.tsv"
    TOTAL=$("$SAMTOOLS" view -c "$BAM")
    VALID=$("$SAMTOOLS" view -c -q 30 -F 3844 "$BAM")
    MITO=$(awk '$1=="NC_007936.1" {n+=$3} END {print n+0}' "$RUN_SAVE/filtered.idxstats.tsv")
    printf 'metric\tvalue\ntotal_records\t%s\npassing_filter_records\t%s\nmitochondrial_records\t%s\n' "$TOTAL" "$VALID" "$MITO" > "$RUN_SAVE/bam_validation.tsv"
    ((TOTAL>0 && TOTAL==VALID && MITO==0)) || { echo 'Filtered BAM validation failed' >&2; exit 4; }
    python atac_job_checks.py tss "$RUN" > "$RUN_SAVE/tss_validation.txt"

    ITEMS=("results/atacseq/filtered_bam/$RUN" "$QC" "results/atacseq/bowtie2/$RUN/bowtie2.log" "results/preprocessing/fastp/reports/$RUN")
    BYTES=$(du -s -B1 -- "${ITEMS[@]}" | awk '{n+=$1} END {printf "%.0f\n", n}')
    check_capacity "$((BYTES + 256*1024*1024))" > "$RUN_SAVE/home_capacity_before_copy.txt" || { echo 'Not enough HOME capacity to preserve this run' >&2; exit 5; }
    mkdir -p "$RUN_SAVE/artifacts"
    for ITEM in "${ITEMS[@]}"; do rsync -aR "$ITEM" "$RUN_SAVE/artifacts/"; done
    # Hash every persisted result against its scratch source, not just the BAM.
    find "${ITEMS[@]}" -type f -print0 | sort -z | xargs -0 sha256sum > "$RUN_SAVE/SHA256SUMS.txt"
    (cd "$RUN_SAVE/artifacts" && sha256sum --check "$RUN_SAVE/SHA256SUMS.txt") > "$RUN_SAVE/persistent_checksum_verification.txt"
    printf 'BAM and TSS artifacts verified; analytical suitability still requires review.\n' > "$RUN_SAVE/complete.ok"
    cat "$QC/tss_enrichment_summary.tsv"
    du -s -B1 "$SCRATCH" > "$RUN_SAVE/scratch_usage_bytes.txt"
    printf 'DONE %s %s\n' "$RUN" "$(date --iso-8601=seconds)"

    # Delete only this job's scratch products, after verifying persistent copies.
    python atac_job_checks.py clean "$RUN" "$WORK" "$SCRATCH"
done

mkdir -p "$SAVE/reference_provenance"
find resources/reference -type f \( -name '*metadata*.tsv' -o -name SHA256SUMS.txt -o -name tss_summary.tsv -o -name tss.bed \) -print0 |
    while IFS= read -r -d '' ITEM; do rsync -aR "$ITEM" "$SAVE/reference_provenance/"; done
find resources/reference -type f ! -name .snakemake_timestamp -print0 | sort -z | xargs -0 sha256sum > "$SAVE/reference_SHA256SUMS.txt"
df -h "$SCRATCH" > "$SAVE/scratch_at_end.txt"
printf 'Both runs passed technical checks and persistence verification.\n' > "$SAVE/all_runs_verified.ok"
