#!/usr/bin/env bash
#SBATCH --job-name=chipseq_lean_production
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
set -euo pipefail

: "${TMPDIR:?Node-local TMPDIR is required}"
: "${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
: "${CHIP_ANALYSIS_ID:?One explicit analysis is required}"
: "${CHIP_CONTROL_RUN:?The authoritative control accession is required}"
: "${CHIP_OUTPUT_DIR:?An explicit persistent analysis directory is required}"
: "${CHIP_PRODUCTION_ROOT:?An explicit persistent production root is required}"
: "${CHIP_PRODUCTION_SLOT:?A reserved production slot is required}"
: "${CHIP_PRODUCTION_SLOT_TOKEN:?A production slot token is required}"
[[ "$TMPDIR" = /* && -d "$TMPDIR" && -w "$TMPDIR" ]] || { echo "ERROR: invalid TMPDIR" >&2; exit 2; }
[[ "$CHIP_OUTPUT_DIR" = /* && "$CHIP_PRODUCTION_ROOT" = /* ]] || { echo "ERROR: persistent paths must be absolute" >&2; exit 2; }
[[ "$CHIP_PRODUCTION_SLOT" == 1 || "$CHIP_PRODUCTION_SLOT" == 2 ]] || { echo "ERROR: production slot must be 1 or 2" >&2; exit 2; }

SOURCE="$PWD"
PLAN="$CHIP_OUTPUT_DIR/submission_plan.json"
ART_FINAL="$CHIP_OUTPUT_DIR/artifacts"
ART="$CHIP_OUTPUT_DIR/.artifacts.part"
LOG="$ART/processing_logs"
WORK="$TMPDIR/chipseq-production-${CHIP_ANALYSIS_ID}-${SLURM_JOB_ID}"
CONTROL_DIR="$CHIP_PRODUCTION_ROOT/shared_controls/$CHIP_CONTROL_RUN"
STAGE=initialization
START_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
MONITOR_PID=""
SLOT_ACTIVE=false
SLOT_LOCK="$CHIP_PRODUCTION_ROOT/.production.slot.${CHIP_PRODUCTION_SLOT}.lock"
exec 9>"$SLOT_LOCK"
flock -n 9 || { echo "ERROR: reserved ChIP production slot is already active" >&2; exit 2; }

stop_monitor() { if [[ -n "$MONITOR_PID" ]]; then kill "$MONITOR_PID" 2>/dev/null || true; wait "$MONITOR_PID" 2>/dev/null || true; MONITOR_PID=""; fi; }
release_slot() { if [[ "$SLOT_ACTIVE" == true ]]; then python "$SOURCE/workflow/scripts/chipseq_production.py" release-slot --production-root "$CHIP_PRODUCTION_ROOT" --slot "$CHIP_PRODUCTION_SLOT" --analysis "$CHIP_ANALYSIS_ID" --token "$CHIP_PRODUCTION_SLOT_TOKEN"; SLOT_ACTIVE=false; fi; }
finish() { status=$?; trap - EXIT; set +e; stop_monitor; printf 'stage\t%s\nexit_status\t%s\n' "$STAGE" "$status" > "$CHIP_OUTPUT_DIR/job_status.tsv"; release_slot; exit "$status"; }
trap finish EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
python workflow/scripts/chipseq_production.py validate-slot --production-root "$CHIP_PRODUCTION_ROOT" --slot "$CHIP_PRODUCTION_SLOT" --analysis "$CHIP_ANALYSIS_ID" --token "$CHIP_PRODUCTION_SLOT_TOKEN"
SLOT_ACTIVE=true
[[ ! -e "$ART" && ! -e "$ART_FINAL" ]] || { echo "ERROR: persistent artifact collision" >&2; exit 2; }
mkdir -p "$WORK" "$LOG" "$ART/software_versions"
RESERVE="$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["scratch_structural_reserve_bytes"])' "$PLAN")"
FREE="$(df -PB1 "$TMPDIR" | awk 'NR==2 {print $4}')"
(( FREE >= RESERVE )) || { echo "ERROR: TMPDIR free bytes $FREE are below required reserve $RESERVE" >&2; exit 2; }
record_disk() { printf '%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" "$(du -sb "$WORK" | awk '{print $1}')" "$(df -PB1 "$TMPDIR" | awk 'NR==2 {print $4}')" >> "$ART/disk_usage.tsv"; }
printf 'timestamp_utc\tphase\tscratch_used_bytes\tscratch_filesystem_free_bytes\n' > "$ART/disk_usage.tsv"
printf 'epoch_seconds\tscratch_used_bytes\n' > "$ART/scratch_monitor.tsv"
record_disk job_start
(while true; do printf '%s\t%s\n' "$(date +%s)" "$(du -sb "$WORK" | awk '{print $1}')" >> "$ART/scratch_monitor.tsv"; sleep 10; done) & MONITOR_PID=$!

sha256sum --check ../source.sha256 > "$LOG/source_validation.log"
rsync -a --exclude='.git/' "$SOURCE/" "$WORK/project/"
cd "$WORK/project"
module purge
module load Miniforge3/24.1.2-0
export CONDA_PKGS_DIRS="$TMPDIR/chipseq-production-conda-pkgs"
export XDG_CACHE_HOME="$TMPDIR/chipseq-production-cache"
export CONDA_CHANNEL_PRIORITY=strict PYTHONDONTWRITEBYTECODE=1 JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
export R_TEMPDIR="$TMPDIR" TMP="$TMPDIR" TEMP="$TMPDIR"
mkdir -p "$CONDA_PKGS_DIRS" "$XDG_CACHE_HOME"
ENV="$TMPDIR/chipseq-production-env"
STAGE=environment_setup
conda env create --prefix "$ENV" --file workflow/envs/chipseq_production.yaml > "$LOG/environment_creation.log" 2>&1
export PATH="$ENV/bin:$PATH"
python --version > "$ART/software_versions/python.txt" 2>&1
fastp --version > "$ART/software_versions/fastp.txt" 2>&1
fastqc --version > "$ART/software_versions/fastqc.txt" 2>&1
bowtie2 --version > "$ART/software_versions/bowtie2.txt" 2>&1
samtools --version > "$ART/software_versions/samtools.txt" 2>&1
conda list --json --prefix "$ENV" \
    | "$ENV/bin/python" workflow/scripts/chipseq_benchmark.py package-provenance --package picard \
    > "$ART/software_versions/picard.txt"
macs3 --version > "$ART/software_versions/macs3.txt" 2>&1
bedtools --version > "$ART/software_versions/bedtools.txt" 2>&1
conda list --explicit --prefix "$ENV" > "$ART/software_versions/environment_explicit.txt"
record_disk environment_ready

python workflow/scripts/chipseq_production.py fastq-manifest --plan "$PLAN" --output "$WORK/fastq_manifest.tsv"
REF="$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["reference"]["root"])' "$PLAN")"
python workflow/scripts/chipseq_benchmark.py check-reference --shared-reference-root "$REF" > "$LOG/reference_validation.json"
INDEX="$REF/bowtie2_index/genome_plus_mt"
mkdir -p inputs/reference config/runtime
cp "$REF/genome_plus_mt.fa.fai" inputs/reference/genome_plus_mt.fa.fai

process_run() {
    run="$1"; role="$2"; tag="${run}_${role}"; root="$WORK/runs/$run"
    raw="$root/raw/$run.fastq.gz"; processed="$root/preprocessed/$run.fastq.gz"
    align="$root/alignment"; filtered="$root/filtered"
    mkdir -p "$root/raw" "$root/preprocessed" "$root/fastqc" "$root/picard_tmp" "$align" "$filtered"
    STAGE="${role}_${run}_download"
    python workflow/scripts/download_fastq.py --manifest "$WORK/fastq_manifest.tsv" --run-accession "$run" --filename "$run.fastq.gz" --output "$raw" > "$LOG/${tag}_download.log" 2>&1
    gzip -t "$raw"; record_disk "${role}_${run}_after_download"
    STAGE="${role}_${run}_preprocessing"
    fastp --in1 "$raw" --out1 "$processed" --json "$root/fastp.json" --html "$root/fastp.html" --thread 4 -Q -L -G > "$LOG/${tag}_fastp.log" 2>&1
    gzip -t "$processed"
    fastqc --threads 1 --outdir "$root/fastqc" "$processed" > "$LOG/${tag}_fastqc.log" 2>&1
    python workflow/scripts/chipseq_production.py alignment-plan --plan "$PLAN" --run "$run" --role "$role" --fastq "$processed" --fastp-json "$root/fastp.json" --output "config/runtime/alignment_${tag}.json"
    rm -- "$raw"
    cp "config/runtime/alignment_${tag}.json" config/chipseq_alignment_inputs.json
    STAGE="${role}_${run}_alignment"
    bowtie2 --local --very-sensitive-local --seed 0 -p 6 -x "$INDEX" -U "$processed" --rg-id "$run" --rg "SM:$run" --rg "LB:$run" --rg PL:ILLUMINA 2> "$align/bowtie2.log" | samtools sort -@ 1 -m 768M -T "$root/sort_tmp" -o "$align/raw.sorted.bam" - 2> "$align/sort.log"
    samtools quickcheck -v "$align/raw.sorted.bam"; samtools index -c "$align/raw.sorted.bam"
    samtools flagstat "$align/raw.sorted.bam" > "$align/flagstat.txt"
    samtools idxstats "$align/raw.sorted.bam" > "$align/idxstats.tsv"
    samtools stats "$align/raw.sorted.bam" > "$align/samtools_stats.txt"
    python workflow/scripts/chipseq_alignment_support.py qc --run "$run" --bam "$align/raw.sorted.bam" --output "$align/alignment_qc.json"
    rm -- "$processed"; record_disk "${role}_${run}_after_alignment"
    python workflow/scripts/chipseq_production.py filtering-plan --plan "$PLAN" --run "$run" --role "$role" --bam "$align/raw.sorted.bam" --qc "$align/alignment_qc.json" --output "config/runtime/filtering_${tag}.json"
    cp "config/runtime/filtering_${tag}.json" config/chipseq_filtering_inputs.json
    STAGE="${role}_${run}_duplicate_marking"
    picard -Xmx12000m -XX:ActiveProcessorCount=4 MarkDuplicates INPUT="$align/raw.sorted.bam" OUTPUT="$root/dupmarked.sorted.bam" METRICS_FILE="$filtered/picard_metrics.txt" TMP_DIR="$root/picard_tmp" REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=false READ_NAME_REGEX=null DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES CREATE_INDEX=false VALIDATION_STRINGENCY=STRICT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=256 > "$filtered/picard.log" 2>&1
    samtools quickcheck -v "$root/dupmarked.sorted.bam"
    rm -- "$align/raw.sorted.bam" "$align/raw.sorted.bam.csi"; record_disk "${role}_${run}_after_duplicate_marking"
    STAGE="${role}_${run}_filtering"
    python workflow/scripts/chipseq_filtering_support.py filter --run "$run" --bam "$root/dupmarked.sorted.bam" --metrics "$filtered/picard_metrics.txt" --output "$filtered/filtered.bam" --report "$filtered/filtering_qc.json" --flow "$filtered/filtering_flow.tsv"
    samtools quickcheck -v "$filtered/filtered.bam"; samtools index -c "$filtered/filtered.bam"
    python workflow/scripts/chipseq_filtering_support.py verify --run "$run" --bam "$filtered/filtered.bam" --report "$filtered/filtering_qc.json" --output "$filtered/filtered_validation.json"
    rm -- "$root/dupmarked.sorted.bam"
    cp "$align/alignment_qc.json" "$filtered/alignment_qc.json"
    cp "$root/fastp.json" "$root/fastp.html" "$filtered/"
    cp -a "$root/fastqc" "$filtered/fastqc"
    record_disk "${role}_${run}_after_filtering"
}

STAGE=control_validation_or_creation
if [[ -e "$CONTROL_DIR" ]]; then
    mkdir -p "$CHIP_PRODUCTION_ROOT/.control_locks"
    exec 8>"$CHIP_PRODUCTION_ROOT/.control_locks/${CHIP_CONTROL_RUN}.lock"
    flock -s 8
    [[ -e "$CONTROL_DIR" ]] || { echo "ERROR: shared control disappeared while acquiring its reader lock" >&2; exit 2; }
    python workflow/scripts/chipseq_production.py load-control --plan "$PLAN" --destination "$CONTROL_DIR" > "$LOG/control_reuse_validation.json"
else
    mkdir -p "$CHIP_PRODUCTION_ROOT/.control_locks"
    exec 8>"$CHIP_PRODUCTION_ROOT/.control_locks/${CHIP_CONTROL_RUN}.lock"
    flock -x 8
    if [[ -e "$CONTROL_DIR" ]]; then
        python workflow/scripts/chipseq_production.py load-control --plan "$PLAN" --destination "$CONTROL_DIR" > "$LOG/control_reuse_validation.json"
    else
        process_run "$CHIP_CONTROL_RUN" input
        python workflow/scripts/chipseq_production.py persist-control --plan "$PLAN" --source "$WORK/runs/$CHIP_CONTROL_RUN/filtered" --destination "$CONTROL_DIR"
        rm -- "$WORK/runs/$CHIP_CONTROL_RUN/filtered/filtered.bam" "$WORK/runs/$CHIP_CONTROL_RUN/filtered/filtered.bam.csi"
    fi
    flock -s 8
fi
IP="$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["analysis"]["ip_run_accession"])' "$PLAN")"
process_run "$IP" ip
IP_DIR="$WORK/runs/$IP/filtered"
python workflow/scripts/chipseq_production.py make-runtime --plan "$PLAN" --ip-dir "$IP_DIR" --control-dir "$CONTROL_DIR" --output "$WORK/runtime.json"

PEAK_ROOT="$WORK/analysis"; PEAKS="$PEAK_ROOT/peaks"; FRAGMENT="$PEAK_ROOT/fragment"; FRIP="$PEAK_ROOT/frip"
mkdir -p "$PEAKS" "$FRAGMENT" "$FRIP"
POLICY="$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["analysis"]["fragment_size_policy"])' "$PLAN")"
FRAGMENT_ARG=()
if [[ "$POLICY" == phantompeakqualtools ]]; then
    STAGE=phantompeakqualtools
    run_spp.R -c="$IP_DIR/filtered.bam" -p=4 -savp="$FRAGMENT/cross_correlation.pdf.part" -out="$FRAGMENT/phantompeakqualtools.tsv.part" > "$FRAGMENT/phantompeakqualtools.log" 2>&1
    test -s "$FRAGMENT/cross_correlation.pdf.part"; test -s "$FRAGMENT/phantompeakqualtools.tsv.part"
    mv "$FRAGMENT/cross_correlation.pdf.part" "$FRAGMENT/cross_correlation.pdf"; mv "$FRAGMENT/phantompeakqualtools.tsv.part" "$FRAGMENT/phantompeakqualtools.tsv"
    python workflow/scripts/chipseq_peak_calling_dynamic.py parse-phantom --table "$FRAGMENT/phantompeakqualtools.tsv" --runtime "$WORK/runtime.json" --analysis "$CHIP_ANALYSIS_ID" --output "$FRAGMENT/fragment_estimate.json"
    FRAGMENT_ARG=(--fragment "$FRAGMENT/fragment_estimate.json")
elif [[ "$POLICY" != fixed ]]; then echo "ERROR: unsupported fragment policy" >&2; exit 2; fi
python workflow/scripts/chipseq_peak_calling_dynamic.py analysis-parameters --runtime "$WORK/runtime.json" --analysis "$CHIP_ANALYSIS_ID" "${FRAGMENT_ARG[@]}" --output "$PEAK_ROOT/parameters.json"
python workflow/scripts/chipseq_production.py fragment-record --plan "$PLAN" --parameters "$PEAK_ROOT/parameters.json" --output "$FRAGMENT/fragment_size_provenance.json"
STAGE=peak_calling
python workflow/scripts/chipseq_peak_calling_dynamic.py call-macs3 --runtime "$WORK/runtime.json" --parameters "$PEAK_ROOT/parameters.json" --analysis "$CHIP_ANALYSIS_ID" --output-dir "$PEAKS" --log "$PEAK_ROOT/macs3.log" --version "$PEAK_ROOT/macs3_version.txt" "${FRAGMENT_ARG[@]}"
MODE="$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["analysis"]["peak_mode"])' "$PLAN")"
if [[ "$MODE" == narrow ]]; then PRIMARY="$PEAKS/${CHIP_ANALYSIS_ID}_peaks.narrowPeak"; SECONDARY="$PEAKS/${CHIP_ANALYSIS_ID}_summits.bed"; elif [[ "$MODE" == broad ]]; then PRIMARY="$PEAKS/${CHIP_ANALYSIS_ID}_peaks.broadPeak"; SECONDARY="$PEAKS/${CHIP_ANALYSIS_ID}_peaks.gappedPeak"; else exit 2; fi
samtools view -c "$IP_DIR/filtered.bam" > "$FRIP/ip_total_reads.txt"; samtools view -c "$CONTROL_DIR/filtered.bam" > "$FRIP/input_total_reads.txt"
if [[ -s "$PRIMARY" ]]; then bedtools intersect -abam "$IP_DIR/filtered.bam" -b "$PRIMARY" -u -bed | wc -l | tr -d ' ' > "$FRIP/ip_reads_in_peaks.txt"; bedtools intersect -abam "$CONTROL_DIR/filtered.bam" -b "$PRIMARY" -u -bed | wc -l | tr -d ' ' > "$FRIP/input_reads_overlapping_ip_peaks.txt"; else printf '0\n' > "$FRIP/ip_reads_in_peaks.txt"; printf '0\n' > "$FRIP/input_reads_overlapping_ip_peaks.txt"; fi
TREAT="$PEAKS/${CHIP_ANALYSIS_ID}_treat_pileup.bdg"; LAMBDA="$PEAKS/${CHIP_ANALYSIS_ID}_control_lambda.bdg"
python workflow/scripts/chipseq_peak_calling_dynamic.py summarize-analysis --runtime "$WORK/runtime.json" --parameters "$PEAK_ROOT/parameters.json" --analysis "$CHIP_ANALYSIS_ID" --primary-peak "$PRIMARY" --secondary-peak "$SECONDARY" --xls "$PEAKS/${CHIP_ANALYSIS_ID}_peaks.xls" --treat-bdg "$TREAT" --control-bdg "$LAMBDA" --log "$PEAK_ROOT/macs3.log" --version "$PEAK_ROOT/macs3_version.txt" --ip-total "$FRIP/ip_total_reads.txt" --input-total "$FRIP/input_total_reads.txt" --ip-overlap "$FRIP/ip_reads_in_peaks.txt" --input-overlap "$FRIP/input_reads_overlapping_ip_peaks.txt" --qc "$PEAK_ROOT/peak_qc.json" --provenance "$PEAK_ROOT/provenance.json" "${FRAGMENT_ARG[@]}"
record_disk after_peak_calling

gzip -n -c "$TREAT" > "$TREAT.gz.part"; gzip -n -c "$LAMBDA" > "$LAMBDA.gz.part"; gzip -t "$TREAT.gz.part" "$LAMBDA.gz.part"
mv "$TREAT.gz.part" "$TREAT.gz"; mv "$LAMBDA.gz.part" "$LAMBDA.gz"
[[ "$(gzip -dc "$TREAT.gz" | sha256sum | awk '{print $1}')" == "$(sha256sum "$TREAT" | awk '{print $1}')" ]]
[[ "$(gzip -dc "$LAMBDA.gz" | sha256sum | awk '{print $1}')" == "$(sha256sum "$LAMBDA" | awk '{print $1}')" ]]
python workflow/scripts/chipseq_production.py compression-record --input "$TREAT" --gzip "$TREAT.gz" --input "$LAMBDA" --gzip "$LAMBDA.gz" --output "$PEAK_ROOT/signal_compression.json"
rm -- "$TREAT" "$LAMBDA"
STAGE=persistent_output
mkdir -p "$ART/ip_processing_qc" "$ART/provenance"
rsync -a "$PEAK_ROOT/" "$ART/"
rsync -a --exclude='filtered.bam' --exclude='filtered.bam.csi' "$IP_DIR/" "$ART/ip_processing_qc/"
cp "$WORK/runtime.json" "$ART/provenance/runtime_manifest.json"; cp "$PLAN" "$ART/provenance/submission_plan.json"; cp "$CHIP_OUTPUT_DIR/source.sha256" "$ART/provenance/source.sha256"
record_disk final_persistence; stop_monitor
printf 'maximum_observed_scratch_bytes\t%s\n' "$(awk -F '\t' 'NR>1 && $2>m {m=$2} END {print m+0}' "$ART/scratch_monitor.tsv")" > "$ART/maximum_scratch_usage.tsv"
mv "$ART" "$ART_FINAL"
ART="$ART_FINAL"
python workflow/scripts/chipseq_production.py finalize --plan "$PLAN" --artifacts "$ART" --output "$CHIP_OUTPUT_DIR" --start-utc "$START_UTC" --slurm-job "$SLURM_JOB_ID"
STAGE=slot_release
release_slot
STAGE=completed
printf 'stage\tcompleted\nexit_status\t0\n' > "$CHIP_OUTPUT_DIR/job_status.tsv"
trap - EXIT
echo "[OK] Lean ChIP production completed: $CHIP_ANALYSIS_ID"
