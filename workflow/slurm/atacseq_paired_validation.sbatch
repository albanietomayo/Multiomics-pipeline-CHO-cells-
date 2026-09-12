#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --job-name=atac_pe_validation
#SBATCH --output=/cephyr/users/mayoa/Vera/atacseq_paired_validation_%j.out

set -euo pipefail

RUN_ACCESSION="SRR29929613"
PROJECT_HOME="/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_atacseq_paired"
SAVE_ROOT="/cephyr/users/mayoa/Vera/atacseq_paired_validation_results"
JOB_TAG="${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
MEM_MB=90000

WORK="${TMPDIR:?TMPDIR is required}/TFM_multiomics_pipeline_atacseq_${RUN_ACCESSION}"
SAVE="${SAVE_ROOT}/${RUN_ACCESSION}/job_${JOB_TAG}"
MAIN_ENV="${TMPDIR}/cho-multiomics"
CONDA_ENVS="${TMPDIR}/snakemake-conda"
CONDA_PKGS="${TMPDIR}/conda-pkgs"
TARGET="results/atacseq/filtered_bam/${RUN_ACCESSION}/filtered.bam"
FILTERED_BAM="${WORK}/${TARGET}"

if [[ -e "$WORK" ]]; then
    printf 'ERROR: work directory already exists: %s\n' "$WORK" >&2
    exit 2
fi

mkdir -p "$WORK" "$SAVE"

git -C "$PROJECT_HOME" rev-parse HEAD > "$SAVE/pipeline_commit.txt"
git -C "$PROJECT_HOME" status --short > "$SAVE/git_status_at_start.txt"
cp "$PROJECT_HOME/config/validation_runs.tsv" "$SAVE/validation_runs_for_job.tsv"

rsync -a \
    --exclude='.git' \
    --exclude='.snakemake/' \
    --exclude='data/' \
    --exclude='results/' \
    --exclude='resources/' \
    "$PROJECT_HOME/" "$WORK/"

cd "$WORK"

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

set +e
/usr/bin/time -v \
    -o "$SAVE/time_verbose.txt" \
    "$SNAKEMAKE" "$TARGET" \
        --use-conda \
        --conda-prefix "$CONDA_ENVS" \
        --cores "$THREADS" \
        --resources "mem_mb=$MEM_MB" \
        --rerun-incomplete \
        --notemp \
        --latency-wait 60 \
        -p \
        > "$SAVE/snakemake_execution.log" 2>&1
WORKFLOW_STATUS=$?
set -e

printf 'workflow_exit_status\t%s\n' "$WORKFLOW_STATUS" \
    > "$SAVE/workflow_status.tsv"

if (( WORKFLOW_STATUS != 0 )); then
    if [[ -d "$WORK/.snakemake/log" ]]; then
        rsync -a "$WORK/.snakemake/log/" "$SAVE/snakemake_logs/"
    fi
    exit "$WORKFLOW_STATUS"
fi

SAMTOOLS=$(
    find "$CONDA_ENVS" \
        -path '*/bin/samtools' \
        -print \
        -quit
)

if [[ -z "$SAMTOOLS" || ! -x "$SAMTOOLS" ]]; then
    printf 'ERROR: samtools executable not found\n' >&2
    exit 3
fi

"$SAMTOOLS" quickcheck -v "$FILTERED_BAM"
"$SAMTOOLS" flagstat "$FILTERED_BAM" > "$SAVE/filtered.flagstat.txt"
"$SAMTOOLS" idxstats "$FILTERED_BAM" > "$SAVE/filtered.idxstats.tsv"

TOTAL=$("$SAMTOOLS" view -c "$FILTERED_BAM")
PROPER=$("$SAMTOOLS" view -c -f 3 "$FILTERED_BAM")
UNMAPPED=$("$SAMTOOLS" view -c -f 4 "$FILTERED_BAM")
MATE_UNMAPPED=$("$SAMTOOLS" view -c -f 8 "$FILTERED_BAM")
SECONDARY=$("$SAMTOOLS" view -c -f 256 "$FILTERED_BAM")
QCFAIL=$("$SAMTOOLS" view -c -f 512 "$FILTERED_BAM")
DUPLICATE=$("$SAMTOOLS" view -c -f 1024 "$FILTERED_BAM")
SUPPLEMENTARY=$("$SAMTOOLS" view -c -f 2048 "$FILTERED_BAM")
MITOCHONDRIAL=$(
    "$SAMTOOLS" view -c "$FILTERED_BAM" NC_007936.1
)
LOW_MAPQ=$(
    "$SAMTOOLS" view "$FILTERED_BAM" |
    awk '$5 < 30 {count++} END {print count+0}'
)

PAIR_DISTRIBUTION="$SAVE/pair_multiplicity.tsv"
printf 'records_per_qname\tqname_groups\n' > "$PAIR_DISTRIBUTION"

"$SAMTOOLS" sort \
    -n \
    -@ "$THREADS" \
    -m 1G \
    -T "$TMPDIR/atac_pair_validation" \
    -u \
    "$FILTERED_BAM" \
| "$SAMTOOLS" view - \
| cut -f1 \
| uniq -c \
| awk '{groups[$1]++} END {for (number in groups) print number, groups[number]}' \
| sort -n \
>> "$PAIR_DISTRIBUTION"

BAD_QNAME_GROUPS=$(
    awk 'NR > 1 && $1 != 2 {sum += $2} END {print sum+0}' \
        "$PAIR_DISTRIBUTION"
)

{
    printf 'metric\tvalue\n'
    printf 'total_records\t%s\n' "$TOTAL"
    printf 'proper_pair_records\t%s\n' "$PROPER"
    printf 'unmapped_records\t%s\n' "$UNMAPPED"
    printf 'mate_unmapped_records\t%s\n' "$MATE_UNMAPPED"
    printf 'secondary_records\t%s\n' "$SECONDARY"
    printf 'qcfail_records\t%s\n' "$QCFAIL"
    printf 'duplicate_records\t%s\n' "$DUPLICATE"
    printf 'supplementary_records\t%s\n' "$SUPPLEMENTARY"
    printf 'mitochondrial_records\t%s\n' "$MITOCHONDRIAL"
    printf 'mapq_below_30_records\t%s\n' "$LOW_MAPQ"
    printf 'non_pair_qname_groups\t%s\n' "$BAD_QNAME_GROUPS"
} > "$SAVE/filtered_validation_summary.tsv"

if (( TOTAL == 0 ||
      TOTAL != PROPER ||
      UNMAPPED != 0 ||
      MATE_UNMAPPED != 0 ||
      SECONDARY != 0 ||
      QCFAIL != 0 ||
      DUPLICATE != 0 ||
      SUPPLEMENTARY != 0 ||
      MITOCHONDRIAL != 0 ||
      LOW_MAPQ != 0 ||
      BAD_QNAME_GROUPS != 0 )); then
    printf 'ERROR: paired BAM validation failed\n' >&2
    exit 4
fi

mkdir -p "$SAVE/artifacts"

for item in \
    "resources/reference/CriGri-PICRH-1.0" \
    "results/metadata" \
    "results/preprocessing/fastp/reports/${RUN_ACCESSION}" \
    "results/atacseq/bowtie2/${RUN_ACCESSION}/bowtie2.log" \
    "results/atacseq/qc/${RUN_ACCESSION}" \
    "results/atacseq/filtered_bam/${RUN_ACCESSION}"
do
    if [[ -e "$item" ]]; then
        rsync -aR "$item" "$SAVE/artifacts/"
    fi
done

if [[ -d "$WORK/.snakemake/log" ]]; then
    rsync -a "$WORK/.snakemake/log/" "$SAVE/snakemake_logs/"
fi

sha256sum \
    "$TARGET" \
    "${TARGET}.bai" \
    > "$SAVE/filtered_sha256.txt"

du -sh "$WORK" > "$SAVE/final_work_directory_size.txt"
df -h "$TMPDIR" > "$SAVE/tmpdir_filesystem.txt"

printf 'ATAC-seq paired validation completed successfully\n'
