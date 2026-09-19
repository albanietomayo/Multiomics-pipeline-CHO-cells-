#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:30:00
#SBATCH --job-name=atac_manifest_v2
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

SOURCE="/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_git"
SAVE="$HOME/atac_filter_integration_validation/manifest_v2_job_${SLURM_JOB_ID}"

SCRATCH="${TMPDIR:?TMPDIR is required}/atac_manifest_v2_${SLURM_JOB_ID}"
WORK="$SCRATCH/repo"
MAIN_ENV="$SCRATCH/main-env"
RULE_ENVS="$SCRATCH/rule-envs"

TARGET="results/metadata/validation_fastq_manifest.tsv"

mkdir -p "$SAVE" "$WORK" "$RULE_ENVS"

echo "============================================================"
echo "SOURCE REPOSITORY CHECK"
echo "============================================================"

cd "$SOURCE"

[[ "$(git branch --show-current)" == "atacseq-downstream-integration" ]] || {
    echo "ERROR: wrong source branch" >&2
    exit 10
}

git status --short
git diff --check

git diff -- \
  config/config.yaml \
  workflow/rules/atacseq.smk \
  > "$SAVE/source_uncommitted_atac.diff"

echo
echo "============================================================"
echo "CREATE ISOLATED WORKING COPY"
echo "============================================================"

git archive HEAD | tar -x -C "$WORK"

# Overlay exactly the two current uncommitted ATAC files.
cp \
  "$SOURCE/config/config.yaml" \
  "$WORK/config/config.yaml"

cp \
  "$SOURCE/workflow/rules/atacseq.smk" \
  "$WORK/workflow/rules/atacseq.smk"

echo "Isolated repo:"
echo "$WORK"

echo
echo "============================================================"
echo "CREATE MAIN SNAKEMAKE ENVIRONMENT"
echo "============================================================"

module purge
module load Miniforge3/24.1.2-0

export CONDA_PKGS_DIRS="$SCRATCH/conda-pkgs"
export XDG_CACHE_HOME="$SCRATCH/cache"

mkdir -p "$CONDA_PKGS_DIRS" "$XDG_CACHE_HOME"

conda env create \
  --prefix "$MAIN_ENV" \
  --file "$WORK/environment.yml" \
  > "$SAVE/main_environment_creation.log" 2>&1

SNAKEMAKE="$MAIN_ENV/bin/snakemake"

echo "Snakemake: $("$SNAKEMAKE" --version)"

cd "$WORK"

echo
echo "============================================================"
echo "DRY-RUN MANIFEST"
echo "============================================================"

"$SNAKEMAKE" \
  --snakefile Snakefile \
  --cores 1 \
  --dry-run \
  --printshellcmds \
  --use-conda \
  --conda-prefix "$RULE_ENVS" \
  "$TARGET" \
  > "$SAVE/manifest.dryrun.log" 2>&1

cat "$SAVE/manifest.dryrun.log"

echo
echo "============================================================"
echo "VERIFY NO HEAVY ATAC JOBS"
echo "============================================================"

if grep -E \
  'rule (download_validation_fastq|fastp_|atacseq_align_run|atacseq_mark_duplicates|atacseq_filter_bam)' \
  "$SAVE/manifest.dryrun.log"
then
    echo "ERROR: heavy ATAC processing scheduled" >&2
    exit 20
fi

echo "Heavy ATAC processing: NONE"

echo
echo "============================================================"
echo "EXECUTE MANIFEST WITH RULE CONDA ENVS"
echo "============================================================"

"$SNAKEMAKE" \
  --snakefile Snakefile \
  --cores 1 \
  --use-conda \
  --conda-prefix "$RULE_ENVS" \
  "$TARGET" \
  > "$SAVE/manifest.execution.log" 2>&1

test -s "$TARGET"

cp \
  "$TARGET" \
  "$SAVE/validation_fastq_manifest.tsv"

echo "Manifest successfully created."

echo
echo "============================================================"
echo "VALIDATE SRR12774931 / SRR12774932"
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

    fields = set(reader.fieldnames or [])
    missing = required - fields

    if missing:
        raise SystemExit(
            f"ERROR: missing columns: {sorted(missing)}"
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
            f"ERROR: {run} absent from validation manifest"
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

    print()
    print(run)
    print(" rows:", len(rr))
    print(" fastq_structure:", sorted(structures))
    print(" fastq_role:", sorted(roles))
    print(" status:", "PASS" if valid else "FAIL")

    if not valid:
        raise SystemExit(
            f"ERROR: {run} does not resolve as SINGLE"
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

print()
print(
    "Both eligible bulk CHO ATAC runs "
    "resolve unambiguously as SINGLE."
)
PY

echo
echo "============================================================"
echo "VALIDATION TABLE"
echo "============================================================"

cat "$SAVE/run_structure_validation.tsv"

echo
echo "============================================================"
echo "CHECKSUMS"
echo "============================================================"

(
    cd "$SAVE"

    sha256sum \
      validation_fastq_manifest.tsv \
      run_structure_validation.tsv \
      source_uncommitted_atac.diff \
      > SHA256SUMS.txt

    sha256sum -c SHA256SUMS.txt
)

printf '%s\n' \
  "Isolated ATAC manifest validation completed: SRR12774931 and SRR12774932 are SINGLE." \
  > "$SAVE/complete.ok"

echo
echo "============================================================"
echo "SOURCE REPOSITORY FINAL CHECK"
echo "============================================================"

cd "$SOURCE"

git branch --show-current
git status --short
git diff --check

echo
echo "============================================================"
echo "COMPLETE"
echo "============================================================"

cat "$SAVE/complete.ok"

echo
echo "RESULTS:"
echo "$SAVE"
