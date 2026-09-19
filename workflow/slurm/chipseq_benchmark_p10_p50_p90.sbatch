#!/usr/bin/env bash
#SBATCH --job-name=chipseq_resource_benchmark
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00

# The request above is a provisional bootstrap envelope inherited from the
# validated ChIP preprocessing/alignment workers. It is not a benchmark result.

set -euo pipefail

: "${TMPDIR:?Node-local TMPDIR is required}"
: "${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
: "${CHIP_BENCHMARK_CLASS:?Select exactly one of P10, P50, P90}"
: "${CHIP_BENCHMARK_ACCESSION:?Explicit benchmark accession is required}"
: "${CHIP_SHARED_REFERENCE_ROOT:?Explicit absolute shared reference root is required}"

CHIP_SOURCE="$PWD"
CHIP_SUBMISSION_ROOT="$(dirname "$CHIP_SOURCE")"
CHIP_SAVE="$CHIP_SUBMISSION_ROOT/job_${SLURM_JOB_ID}"
CHIP_WORK="$TMPDIR/chipseq-benchmark-${CHIP_BENCHMARK_CLASS}-${CHIP_BENCHMARK_ACCESSION}"
CHIP_SUPPORT="workflow/scripts/chipseq_benchmark.py"
CHIP_MANIFEST="benchmarks/chipseq/2026-09-18/chipseq_benchmark_manifest.tsv"
CHIP_STAGE="initialization"
CHIP_MONITOR_PID=""

mkdir -p "$CHIP_SAVE" "$CHIP_WORK"
if [[ -f "$CHIP_SAVE/benchmark.complete.json" ]]; then
    python "$CHIP_SUPPORT" verify-stage --stage benchmark --marker "$CHIP_SAVE/benchmark.complete.json"
    echo "ERROR: validated completed benchmark already exists: $CHIP_SAVE" >&2
    exit 2
fi

stop_monitor() {
    if [[ -n "$CHIP_MONITOR_PID" ]]; then
        kill "$CHIP_MONITOR_PID" 2>/dev/null || true
        wait "$CHIP_MONITOR_PID" 2>/dev/null || true
        CHIP_MONITOR_PID=""
    fi
}


finish() {
    local status=$?
    trap - EXIT
    set +e
    stop_monitor
    printf 'stage\t%s\nexit_status\t%s\n' "$CHIP_STAGE" "$status" > "$CHIP_SAVE/job_status.tsv"
    exit "$status"
}

trap finish EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

sha256sum --check ../source.sha256 > "$CHIP_SAVE/source_validation.log"
python "$CHIP_SUPPORT" check-manifest > "$CHIP_SAVE/manifest_validation.log"

CHIP_EXPECTED_ACCESSION="$(awk -F '\t' -v class="$CHIP_BENCHMARK_CLASS" '
    NR > 1 && $1 == class {print $2}
' "$CHIP_MANIFEST")"

if [[ -z "$CHIP_EXPECTED_ACCESSION" || "$CHIP_EXPECTED_ACCESSION" != "$CHIP_BENCHMARK_ACCESSION" ]]; then
    echo "ERROR: benchmark class/accession mismatch" >&2
    exit 2
fi

CHIP_RUN="$CHIP_BENCHMARK_ACCESSION"
CHIP_FILENAME="$CHIP_RUN.fastq.gz"
CHIP_RAW="data/raw/$CHIP_FILENAME"
CHIP_PROCESSED="results/preprocessing/$CHIP_RUN.fastq.gz"
CHIP_FASTP_JSON="results/preprocessing/$CHIP_RUN.fastp.json"
CHIP_FASTP_HTML="results/preprocessing/$CHIP_RUN.fastp.html"
CHIP_POST_DIR="results/qc/post/$CHIP_RUN"
CHIP_ALIGN_DIR="results/chipseq/alignment/$CHIP_RUN"
CHIP_FILTER_DIR="results/chipseq/filtering/$CHIP_RUN"


rsync -a --exclude='.git/' --exclude='data/' --exclude='results/' --exclude='resources/' \
    "$CHIP_SOURCE/" "$CHIP_WORK/"
cd "$CHIP_WORK"

python "$CHIP_SUPPORT" preflight \
    --benchmark-class "$CHIP_BENCHMARK_CLASS" \
    --output-root "$CHIP_SAVE/runtime_output_preflight" \
    --scratch-root "$TMPDIR" \
    > "$CHIP_SAVE/storage_preflight.json"

module purge
module load Miniforge3/24.1.2-0
export CONDA_PKGS_DIRS="$TMPDIR/chipseq-benchmark-conda-pkgs"
export XDG_CACHE_HOME="$TMPDIR/chipseq-benchmark-cache"
export CONDA_CHANNEL_PRIORITY=strict
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$CONDA_PKGS_DIRS" "$XDG_CACHE_HOME"

CHIP_MAIN_ENV="$TMPDIR/chipseq-benchmark-main-env"
CHIP_PREP_ENV="$TMPDIR/chipseq-benchmark-preprocessing-env"
CHIP_QC_ENV="$TMPDIR/chipseq-benchmark-qc-env"
CHIP_ALIGN_ENV="$TMPDIR/chipseq-benchmark-alignment-env"
CHIP_DUP_ENV="$TMPDIR/chipseq-benchmark-duplicates-env"
CHIP_FILTER_ENV="$TMPDIR/chipseq-benchmark-filtering-env"

CHIP_STAGE="environment_setup"
CHIP_SETUP_START="$(date +%s)"
conda env create --prefix "$CHIP_MAIN_ENV" --file environment.yml > "$CHIP_SAVE/main_environment.log" 2>&1
conda env create --prefix "$CHIP_PREP_ENV" --file workflow/envs/preprocessing.yaml > "$CHIP_SAVE/preprocessing_environment.log" 2>&1
conda env create --prefix "$CHIP_QC_ENV" --file workflow/envs/qc.yaml > "$CHIP_SAVE/qc_environment.log" 2>&1
conda env create --prefix "$CHIP_ALIGN_ENV" --file workflow/envs/chipseq_alignment.yaml > "$CHIP_SAVE/alignment_environment.log" 2>&1
conda env create --prefix "$CHIP_DUP_ENV" --file workflow/envs/chipseq_duplicates.yaml > "$CHIP_SAVE/duplicates_environment.log" 2>&1
conda env create --prefix "$CHIP_FILTER_ENV" --file workflow/envs/chipseq_filtering.yaml > "$CHIP_SAVE/filtering_environment.log" 2>&1

mkdir -p "$CHIP_SAVE/software_versions"
"$CHIP_MAIN_ENV/bin/python" --version > "$CHIP_SAVE/software_versions/python.txt" 2>&1
"$CHIP_MAIN_ENV/bin/snakemake" --version > "$CHIP_SAVE/software_versions/snakemake.txt"
"$CHIP_PREP_ENV/bin/fastp" --version > "$CHIP_SAVE/software_versions/fastp.txt" 2>&1
"$CHIP_QC_ENV/bin/fastqc" --version > "$CHIP_SAVE/software_versions/fastqc.txt" 2>&1
"$CHIP_QC_ENV/bin/multiqc" --version > "$CHIP_SAVE/software_versions/multiqc.txt" 2>&1
"$CHIP_ALIGN_ENV/bin/bowtie2" --version > "$CHIP_SAVE/software_versions/bowtie2.txt" 2>&1
"$CHIP_ALIGN_ENV/bin/samtools" --version > "$CHIP_SAVE/software_versions/samtools.txt" 2>&1
"$CHIP_DUP_ENV/bin/picard" MarkDuplicates --version > "$CHIP_SAVE/software_versions/picard.txt" 2>&1
conda --version > "$CHIP_SAVE/software_versions/conda.txt" 2>&1
conda list --explicit --prefix "$CHIP_MAIN_ENV" > "$CHIP_SAVE/software_versions/main_environment_explicit.txt"
conda list --explicit --prefix "$CHIP_PREP_ENV" > "$CHIP_SAVE/software_versions/preprocessing_environment_explicit.txt"
conda list --explicit --prefix "$CHIP_QC_ENV" > "$CHIP_SAVE/software_versions/qc_environment_explicit.txt"
conda list --explicit --prefix "$CHIP_ALIGN_ENV" > "$CHIP_SAVE/software_versions/alignment_environment_explicit.txt"
conda list --explicit --prefix "$CHIP_DUP_ENV" > "$CHIP_SAVE/software_versions/duplicates_environment_explicit.txt"
conda list --explicit --prefix "$CHIP_FILTER_ENV" > "$CHIP_SAVE/software_versions/filtering_environment_explicit.txt"
"$CHIP_MAIN_ENV/bin/python" -c \
    'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["source_git_commit"])' \
    "$CHIP_SUBMISSION_ROOT/submission_plan.json" > "$CHIP_SAVE/pipeline_commit.txt"

printf '%s\n' \
    'benchmark_scope=acquisition_validation,fastp,post_fastqc,alignment,filtering' \
    'excluded=phantompeakqualtools,macs3' \
    'reference_setup=shared_prerequisite_not_attributed_to_run_stage' \
    'slurm_request=eight_cpus_32G_8h_provisional_bootstrap_not_a_result' \
    > "$CHIP_SAVE/benchmark_contract.txt"

CHIP_STAGE="shared_reference_validation"
python "$CHIP_SUPPORT" check-reference --shared-reference-root "$CHIP_SHARED_REFERENCE_ROOT" \
    --expected "$CHIP_SUBMISSION_ROOT/submission_plan.json" \
    --output "$CHIP_SAVE/reference_provenance.json"
CHIP_INDEX="$CHIP_SHARED_REFERENCE_ROOT/bowtie2_index/genome_plus_mt"
mkdir -p inputs/reference
cp -- "$CHIP_SHARED_REFERENCE_ROOT/genome_plus_mt.fa.fai" \
    inputs/reference/genome_plus_mt.fa.fai.part
cmp -- "$CHIP_SHARED_REFERENCE_ROOT/genome_plus_mt.fa.fai" \
    inputs/reference/genome_plus_mt.fa.fai.part
mv inputs/reference/genome_plus_mt.fa.fai.part inputs/reference/genome_plus_mt.fa.fai

if [[ ! -e "$CHIP_SAVE/benchmark_metrics.tsv" ]]; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    benchmark_class run_accession stage status start_utc end_utc elapsed_seconds exit_status \
    max_rss_kbytes user_cpu_seconds system_cpu_seconds baseline_scratch_bytes peak_scratch_bytes \
    incremental_peak_scratch_bytes input_bytes output_bytes \
    > "$CHIP_SAVE/benchmark_metrics.tsv"
fi

run_stage() {
    local stage="$1"
    shift
    local input_bytes="$1"
    shift
    local output_count="$1"
    shift
    local outputs=("${@:1:$output_count}")
    shift "$output_count"
    local stage_dir="$CHIP_SAVE/stages/$stage"
    local marker="$stage_dir/artifacts.json"
    mkdir -p "$stage_dir"


    CHIP_STAGE="$stage"
    local baseline start_epoch end_epoch start_utc end_utc status peak incremental output_bytes
    baseline="$(du -sb "$CHIP_WORK" | awk '{print $1}')"
    start_epoch="$(date +%s)"
    start_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'epoch_seconds\twork_bytes\n' > "$stage_dir/scratch_usage.tsv"
    (
        while true; do
            printf '%s\t%s\n' "$(date +%s)" "$(du -sb "$CHIP_WORK" | awk '{print $1}')" \
                >> "$stage_dir/scratch_usage.tsv"
            sleep 15
        done
    ) &
    CHIP_MONITOR_PID=$!

    set +e
    /usr/bin/time -v -o "$stage_dir/time_verbose.txt" "$@" > "$stage_dir/execution.log" 2>&1
    status=$?
    set -e
    stop_monitor
    printf '%s\t%s\n' "$(date +%s)" "$(du -sb "$CHIP_WORK" | awk '{print $1}')" \
        >> "$stage_dir/scratch_usage.tsv"
    end_epoch="$(date +%s)"
    end_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    peak="$(awk -F '\t' 'NR > 1 && $2 > max {max=$2} END {print max+0}' "$stage_dir/scratch_usage.tsv")"
    if (( peak > baseline )); then incremental=$((peak - baseline)); else incremental=0; fi
    output_bytes=0
    for output in "${outputs[@]}"; do
        if [[ -f "$output" ]]; then output_bytes=$((output_bytes + $(stat -c %s "$output"))); fi
    done
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$CHIP_BENCHMARK_CLASS" "$CHIP_RUN" "$stage" "$([[ "$status" == 0 ]] && echo complete || echo failed)" \
        "$start_utc" "$end_utc" "$((end_epoch - start_epoch))" "$status" \
        "$(awk -F ': ' '/Maximum resident set size/ {gsub(/[[:space:]]/, "", $2); print $2}' "$stage_dir/time_verbose.txt")" \
        "$(awk -F ': ' '/User time \(seconds\)/ {gsub(/[[:space:]]/, "", $2); print $2}' "$stage_dir/time_verbose.txt")" \
        "$(awk -F ': ' '/System time \(seconds\)/ {gsub(/[[:space:]]/, "", $2); print $2}' "$stage_dir/time_verbose.txt")" \
        "$baseline" "$peak" "$incremental" "$input_bytes" "$output_bytes" \
        >> "$CHIP_SAVE/benchmark_metrics.tsv"
    printf 'stage\t%s\nstatus\t%s\nexit_status\t%s\nstart_utc\t%s\nend_utc\t%s\n' \
        "$stage" "$([[ "$status" == 0 ]] && echo complete || echo failed)" "$status" "$start_utc" "$end_utc" \
        > "$stage_dir/stage_status.tsv"
    if (( status != 0 )); then
        return "$status"
    fi
    local record_args=()
    for output in "${outputs[@]}"; do
        record_args+=(--transient "$output")
    done
    "$CHIP_MAIN_ENV/bin/python" "$CHIP_SUPPORT" record-stage --stage "$stage" --marker "$marker" "${record_args[@]}"
}

mkdir -p data/raw results/preprocessing "$CHIP_POST_DIR" "$CHIP_ALIGN_DIR" "$CHIP_FILTER_DIR" inputs

run_stage acquisition_validation 0 1 "$CHIP_RAW" \
    "$CHIP_MAIN_ENV/bin/python" workflow/scripts/download_fastq.py \
    --manifest "$CHIP_MANIFEST" --run-accession "$CHIP_RUN" --filename "$CHIP_FILENAME" --output "$CHIP_RAW"

run_stage fastp "$(stat -c %s "$CHIP_RAW")" 3 "$CHIP_PROCESSED" "$CHIP_FASTP_JSON" "$CHIP_FASTP_HTML" \
    "$CHIP_PREP_ENV/bin/fastp" --in1 "$CHIP_RAW" --out1 "$CHIP_PROCESSED" \
    --json "$CHIP_FASTP_JSON" --html "$CHIP_FASTP_HTML" --thread 4 -Q -L -G

run_stage post_fastqc "$(stat -c %s "$CHIP_PROCESSED")" 2 \
    "$CHIP_POST_DIR/${CHIP_RUN}_fastqc.html" "$CHIP_POST_DIR/${CHIP_RUN}_fastqc.zip" \
    "$CHIP_QC_ENV/bin/fastqc" --threads 1 --outdir "$CHIP_POST_DIR" "$CHIP_PROCESSED"

"$CHIP_MAIN_ENV/bin/python" "$CHIP_SUPPORT" make-alignment-plan --run "$CHIP_RUN" \
    --processed-fastq "$CHIP_PROCESSED" --fastp-json "$CHIP_FASTP_JSON" \
    --output config/chipseq_alignment_inputs.json

CHIP_ALIGN_SCRIPT="$CHIP_SAVE/alignment_stage.sh"
sed -e "s|@RUN@|$CHIP_RUN|g" -e "s|@PROCESSED@|$CHIP_PROCESSED|g" \
    -e "s|@INDEX@|$CHIP_INDEX|g" -e "s|@OUT@|$CHIP_ALIGN_DIR|g" \
    -e "s|@ALIGN_ENV@|$CHIP_ALIGN_ENV|g" > "$CHIP_ALIGN_SCRIPT" <<'ALIGN_SCRIPT'
#!/usr/bin/env bash
set -euo pipefail
@ALIGN_ENV@/bin/bowtie2 --local --very-sensitive-local --seed 0 -p 6 -x @INDEX@ -U @PROCESSED@ \
    --rg-id @RUN@ --rg SM:@RUN@ --rg LB:@RUN@ --rg PL:ILLUMINA 2> @OUT@/bowtie2.log | \
    @ALIGN_ENV@/bin/samtools sort -@ 1 -m 768M -T @OUT@/sort_tmp -o @OUT@/raw.sorted.bam - \
    2> @OUT@/sort.log
@ALIGN_ENV@/bin/samtools quickcheck -v @OUT@/raw.sorted.bam
@ALIGN_ENV@/bin/samtools index -c @OUT@/raw.sorted.bam
@ALIGN_ENV@/bin/samtools flagstat @OUT@/raw.sorted.bam > @OUT@/flagstat.txt
@ALIGN_ENV@/bin/samtools idxstats @OUT@/raw.sorted.bam > @OUT@/idxstats.tsv
@ALIGN_ENV@/bin/samtools stats @OUT@/raw.sorted.bam > @OUT@/samtools_stats.txt
PATH=@ALIGN_ENV@/bin:$PATH @ALIGN_ENV@/bin/python workflow/scripts/chipseq_alignment_support.py qc \
    --run @RUN@ --bam @OUT@/raw.sorted.bam --output @OUT@/alignment_qc.json
ALIGN_SCRIPT
chmod 700 "$CHIP_ALIGN_SCRIPT"

run_stage alignment "$(stat -c %s "$CHIP_PROCESSED")" 3 \
    "$CHIP_ALIGN_DIR/raw.sorted.bam" "$CHIP_ALIGN_DIR/raw.sorted.bam.csi" "$CHIP_ALIGN_DIR/alignment_qc.json" \
    "$CHIP_ALIGN_SCRIPT"

"$CHIP_MAIN_ENV/bin/python" "$CHIP_SUPPORT" make-filtering-plan --run "$CHIP_RUN" \
    --bam "$CHIP_ALIGN_DIR/raw.sorted.bam" --qc "$CHIP_ALIGN_DIR/alignment_qc.json" \
    --output config/chipseq_filtering_inputs.json

CHIP_FILTER_SCRIPT="$CHIP_SAVE/filtering_stage.sh"
sed -e "s|@RUN@|$CHIP_RUN|g" -e "s|@ALIGN@|$CHIP_ALIGN_DIR|g" \
    -e "s|@OUT@|$CHIP_FILTER_DIR|g" -e "s|@DUP_ENV@|$CHIP_DUP_ENV|g" \
    -e "s|@FILTER_ENV@|$CHIP_FILTER_ENV|g" > "$CHIP_FILTER_SCRIPT" <<'FILTER_SCRIPT'
#!/usr/bin/env bash
set -euo pipefail
mkdir -p inputs/@RUN@ work/chipseq_filtering/@RUN@/picard_tmp
PATH=@FILTER_ENV@/bin:$PATH @FILTER_ENV@/bin/python workflow/scripts/chipseq_filtering_support.py stage \
    --run @RUN@ --output inputs/@RUN@/raw.sorted.bam --report @OUT@/input_validation.json
@DUP_ENV@/bin/picard -Xmx12000m -XX:ActiveProcessorCount=4 MarkDuplicates \
    INPUT=inputs/@RUN@/raw.sorted.bam OUTPUT=work/chipseq_filtering/@RUN@/dupmarked.sorted.bam \
    METRICS_FILE=@OUT@/picard_metrics.txt TMP_DIR=work/chipseq_filtering/@RUN@/picard_tmp \
    REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=false READ_NAME_REGEX=null \
    DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES CREATE_INDEX=false \
    VALIDATION_STRINGENCY=STRICT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=256 > @OUT@/picard.log 2>&1
PATH=@FILTER_ENV@/bin:$PATH @FILTER_ENV@/bin/python workflow/scripts/chipseq_filtering_support.py filter \
    --run @RUN@ --bam work/chipseq_filtering/@RUN@/dupmarked.sorted.bam \
    --metrics @OUT@/picard_metrics.txt --output @OUT@/filtered.bam \
    --report @OUT@/filtering_qc.json --flow @OUT@/filtering_flow.tsv
@FILTER_ENV@/bin/samtools index -c @OUT@/filtered.bam
PATH=@FILTER_ENV@/bin:$PATH @FILTER_ENV@/bin/python workflow/scripts/chipseq_filtering_support.py verify \
    --run @RUN@ --bam @OUT@/filtered.bam --report @OUT@/filtering_qc.json \
    --output @OUT@/filtered_validation.json
FILTER_SCRIPT
chmod 700 "$CHIP_FILTER_SCRIPT"

run_stage filtering "$(stat -c %s "$CHIP_ALIGN_DIR/raw.sorted.bam")" 5 \
    "$CHIP_FILTER_DIR/filtered.bam" "$CHIP_FILTER_DIR/filtered.bam.csi" \
    "$CHIP_FILTER_DIR/filtering_qc.json" "$CHIP_FILTER_DIR/filtered_validation.json" \
    "work/chipseq_filtering/$CHIP_RUN/dupmarked.sorted.bam" \
    "$CHIP_FILTER_SCRIPT"

CHIP_STAGE="successful_persistence"
"$CHIP_MAIN_ENV/bin/python" "$CHIP_SUPPORT" persist-files \
    --destination "$CHIP_SAVE/artifacts" \
    "filtered/filtered.bam=$CHIP_FILTER_DIR/filtered.bam" \
    "filtered/filtered.bam.csi=$CHIP_FILTER_DIR/filtered.bam.csi" \
    "qc/fastp/$CHIP_RUN.fastp.json=$CHIP_FASTP_JSON" \
    "qc/fastp/$CHIP_RUN.fastp.html=$CHIP_FASTP_HTML" \
    "qc/fastqc/${CHIP_RUN}_fastqc.html=$CHIP_POST_DIR/${CHIP_RUN}_fastqc.html" \
    "qc/fastqc/${CHIP_RUN}_fastqc.zip=$CHIP_POST_DIR/${CHIP_RUN}_fastqc.zip" \
    "qc/alignment/alignment_qc.json=$CHIP_ALIGN_DIR/alignment_qc.json" \
    "qc/alignment/flagstat.txt=$CHIP_ALIGN_DIR/flagstat.txt" \
    "qc/alignment/idxstats.tsv=$CHIP_ALIGN_DIR/idxstats.tsv" \
    "qc/alignment/samtools_stats.txt=$CHIP_ALIGN_DIR/samtools_stats.txt" \
    "qc/alignment/bowtie2.log=$CHIP_ALIGN_DIR/bowtie2.log" \
    "qc/filtering/input_validation.json=$CHIP_FILTER_DIR/input_validation.json" \
    "qc/filtering/picard_metrics.txt=$CHIP_FILTER_DIR/picard_metrics.txt" \
    "qc/filtering/picard.log=$CHIP_FILTER_DIR/picard.log" \
    "qc/filtering/filtering_qc.json=$CHIP_FILTER_DIR/filtering_qc.json" \
    "qc/filtering/filtering_flow.tsv=$CHIP_FILTER_DIR/filtering_flow.tsv" \
    "qc/filtering/filtered_validation.json=$CHIP_FILTER_DIR/filtered_validation.json" \
    "provenance/chipseq_alignment_inputs.json=config/chipseq_alignment_inputs.json" \
    "provenance/chipseq_filtering_inputs.json=config/chipseq_filtering_inputs.json" \
    > "$CHIP_SAVE/persistence.log"


mapfile -t CHIP_COMPLETION_OUTPUTS < <(
    find "$CHIP_SAVE/artifacts" -type f -print | sort
    find "$CHIP_SAVE/stages" -name artifacts.json -type f -print | sort
    printf '%s\n' \
        "$CHIP_SAVE/benchmark_metrics.tsv" \
        "$CHIP_SAVE/reference_provenance.json"
)
CHIP_COMPLETION_PENDING="$CHIP_SAVE/benchmark.complete.pending.json"
"$CHIP_MAIN_ENV/bin/python" "$CHIP_SUPPORT" complete-stage --stage benchmark \
    --marker "$CHIP_COMPLETION_PENDING" "${CHIP_COMPLETION_OUTPUTS[@]}"

CHIP_STAGE="completed"
printf 'benchmark_class\t%s\nrun_accession\t%s\nstage\tcompleted\nexit_status\t0\n' \
    "$CHIP_BENCHMARK_CLASS" "$CHIP_RUN" > "$CHIP_SAVE/job_status.tsv"
CHIP_PERSISTENT_BASE_BYTES="$(find "$CHIP_SUBMISSION_ROOT" -type f \
    ! -path "$CHIP_SAVE/persistent_footprint.tsv" -printf '%s\n' | \
    awk '{total += $1} END {print total+0}')"
CHIP_PERSISTENT_BYTES="$CHIP_PERSISTENT_BASE_BYTES"
while true; do
    CHIP_FOOTPRINT_TEXT="$(printf 'scope\tlogical_bytes\nsubmission_directory_final\t%s' \
        "$CHIP_PERSISTENT_BYTES")"
    CHIP_MEASURED_BYTES="$((CHIP_PERSISTENT_BASE_BYTES + ${#CHIP_FOOTPRINT_TEXT} + 1))"
    [[ "$CHIP_MEASURED_BYTES" == "$CHIP_PERSISTENT_BYTES" ]] && break
    CHIP_PERSISTENT_BYTES="$CHIP_MEASURED_BYTES"
done
printf '%s\n' "$CHIP_FOOTPRINT_TEXT" > "$CHIP_SAVE/persistent_footprint.tsv"
mv "$CHIP_COMPLETION_PENDING" "$CHIP_SAVE/benchmark.complete.json"
trap - EXIT
echo "[OK] ChIP-seq benchmark completed: $CHIP_BENCHMARK_CLASS $CHIP_RUN"
