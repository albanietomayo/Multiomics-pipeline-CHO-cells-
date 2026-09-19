#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=atac_manifest_check
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO="/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_git"
SAVE="$HOME/atac_filter_integration_validation/manifest_job_${SLURM_JOB_ID}"
SCRATCH="${TMPDIR:?TMPDIR is required}/atac_manifest_validation_${SLURM_JOB_ID}"

TARGET="results/metadata/validation_fastq_manifest.tsv"

mkdir -p "$SAVE" "$SCRATCH"

cd "$REPO"

echo "============================================================"
echo "REPOSITORY CHECK"
echo "============================================================"

[[ "$(git branch --show-current)" == "atacseq-downstream-integration" ]] || {
    echo "ERROR: wrong branch" >&2
    exit 10
}

git status --short
git diff --check

echo
echo "============================================================"
echo "CREATE SNAKEMAKE ENVIRONMENT"
echo "============================================================"

module purge
module load Miniforge3/24.1.2-0

export CONDA_PKGS_DIRS="$SCRATCH/conda-pkgs"
export XDG_CACHE_HOME="$SCRATCH/cache"

mkdir -p "$CONDA_PKGS_DIRS" "$XDG_CACHE_HOME"

MAIN_ENV="$SCRATCH/main-env"

conda env create \
  --prefix "$MAIN_ENV" \
  --file environment.yml \
  > "$SAVE/environment_creation.log" 2>&1

SNAKEMAKE="$MAIN_ENV/bin/snakemake"

echo "Snakemake: $("$SNAKEMAKE" --version)"

echo
echo "============================================================"
echo "DRY-RUN MANIFEST TARGET"
echo "============================================================"

"$SNAKEMAKE" \
  --snakefile Snakefile \
  --cores 1 \
  --dry-run \
  --printshellcmds \
  "$TARGET" \
  > "$SAVE/manifest.dryrun.log" 2>&1

cat "$SAVE/manifest.dryrun.log"

echo
echo "============================================================"
echo "CHECK NO HEAVY ATAC PROCESSING IS SCHEDULED"
echo "============================================================"

if grep -E \
  'rule (download_validation_fastq|fastp_|atacseq_align_run|atacseq_mark_duplicates|atacseq_filter_bam)' \
  "$SAVE/manifest.dryrun.log"
then
    echo "ERROR: manifest target unexpectedly schedules heavy processing" >&2
    exit 20
fi

echo "No FASTQ download/alignment/filtering jobs detected."

echo
echo "============================================================"
echo "MATERIALIZE VALIDATION MANIFEST"
echo "============================================================"

"$SNAKEMAKE" \
  --snakefile Snakefile \
  --cores 1 \
  "$TARGET" \
  > "$SAVE/manifest.execution.log" 2>&1

test -s "$TARGET"

cp "$TARGET" "$SAVE/validation_fastq_manifest.tsv"

echo "Manifest created:"
ls -lh "$TARGET"

echo
echo "============================================================"
echo "VALIDATE TARGET RUN STRUCTURES"
echo "============================================================"

python3 - "$TARGET" "$SAVE/run_structure_validation.tsv" <<'PY'
import csv
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
output = Path(sys.argv[2])

targets = ["SRR12774931", "SRR12774932"]

with manifest.open(newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    required = {
        "run_accession",
        "fastq_structure",
        "fastq_role",
    }

    missing = required - set(reader.fieldnames or [])

    if missing:
        raise SystemExit(
            f"ERROR: required columns missing: {sorted(missing)}"
        )

    rows = list(reader)

results = []

for run in targets:
    rr = [
        row for row in rows
        if row["run_accession"] == run
    ]

    if not rr:
        raise SystemExit(
            f"ERROR: {run} not found in validation manifest"
        )

    structures = {
        row["fastq_structure"]
        for row in rr
    }

    roles = {
        row["fastq_role"]
        for row in rr
    }

    valid = (
        structures == {"SINGLE"}
        and roles == {"SINGLE"}
    )

    print(f"{run}")
    print(f"  rows: {len(rr)}")
    print(f"  fastq_structure: {sorted(structures)}")
    print(f"  fastq_role: {sorted(roles)}")
    print(f"  classification: {'SINGLE VALID' if valid else 'INVALID'}")

    if not valid:
        raise SystemExit(
            f"ERROR: {run} does not resolve unambiguously as SINGLE"
        )

    results.append(
        (
            run,
            len(rr),
            ",".join(sorted(structures)),
            ",".join(sorted(roles)),
            "PASS",
        )
    )

with output.open("w") as fh:
    fh.write(
        "run_accession\trows\tfastq_structure\tfastq_role\tstatus\n"
    )

    for row in results:
        fh.write("\t".join(map(str, row)) + "\n")

print(
    "\nBoth eligible bulk CHO ATAC runs resolve "
    "unambiguously as SINGLE."
)
PY

echo
echo "============================================================"
echo "VALIDATION TABLE"
echo "============================================================"

cat "$SAVE/run_structure_validation.tsv"

echo
echo "============================================================"
echo "CHECKSUM"
echo "============================================================"

(
    cd "$SAVE"

    sha256sum \
      validation_fastq_manifest.tsv \
      run_structure_validation.tsv \
      > SHA256SUMS.txt

    sha256sum -c SHA256SUMS.txt
)

printf '%s\n' \
  "ATAC validation manifest materialized and both target CHO runs resolve as SINGLE." \
  > "$SAVE/complete.ok"

echo
echo "============================================================"
echo "COMPLETE"
echo "============================================================"

cat "$SAVE/complete.ok"

echo
echo "RESULTS:"
echo "$SAVE"
