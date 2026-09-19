#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=atac_manifest_v3
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

SOURCE="/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_git"
SAVE="$HOME/atac_filter_integration_validation/manifest_v3_job_${SLURM_JOB_ID}"

SCRATCH="${TMPDIR:?TMPDIR is required}/atac_manifest_v3_${SLURM_JOB_ID}"
WORK="$SCRATCH/repo"
PYENV="$SCRATCH/python-env"

mkdir -p "$SAVE" "$WORK"

echo "============================================================"
echo "SOURCE REPOSITORY CHECK"
echo "============================================================"

cd "$SOURCE"

[[ "$(git branch --show-current)" == "atacseq-downstream-integration" ]] || {
    echo "ERROR: wrong branch" >&2
    exit 10
}

STATUS="$(git status --short)"
printf '%s\n' "$STATUS"

EXPECTED=$' M config/config.yaml\n M workflow/rules/atacseq.smk'

[[ "$STATUS" == "$EXPECTED" ]] || {
    echo "ERROR: unexpected source repository state" >&2
    exit 11
}

git diff --check

echo
echo "============================================================"
echo "CREATE ISOLATED SNAPSHOT"
echo "============================================================"

git archive HEAD | tar -x -C "$WORK"

# Overlay only the two current uncommitted ATAC files.
cp \
  "$SOURCE/config/config.yaml" \
  "$WORK/config/config.yaml"

cp \
  "$SOURCE/workflow/rules/atacseq.smk" \
  "$WORK/workflow/rules/atacseq.smk"

test -s "$WORK/config/samples.tsv"

echo "Versioned samples snapshot:"
ls -lh "$WORK/config/samples.tsv"
wc -l "$WORK/config/samples.tsv"

echo
echo "============================================================"
echo "CREATE MINIMAL PYTHON ENVIRONMENT"
echo "============================================================"

module purge
module load Miniforge3/24.1.2-0

export CONDA_PKGS_DIRS="$SCRATCH/conda-pkgs"
export XDG_CACHE_HOME="$SCRATCH/cache"

mkdir -p "$CONDA_PKGS_DIRS" "$XDG_CACHE_HOME"

conda create \
  -y \
  --prefix "$PYENV" \
  -c conda-forge \
  python=3.11 \
  pandas \
  pyyaml \
  > "$SAVE/python_environment_creation.log" 2>&1

PYTHON="$PYENV/bin/python"

"$PYTHON" - <<'PY'
import sys
import pandas

print("Python:", sys.version.split()[0])
print("pandas:", pandas.__version__)
PY

cd "$WORK"

mkdir -p results/metadata

echo
echo "============================================================"
echo "BUILD FASTQ MANIFEST FROM VERSIONED SAMPLES"
echo "============================================================"

"$PYTHON" \
  workflow/scripts/build_fastq_manifest.py \
  --samples config/samples.tsv \
  --output results/metadata/fastq_manifest.tsv \
  > "$SAVE/build_fastq_manifest.log" 2>&1

test -s results/metadata/fastq_manifest.tsv

echo "fastq_manifest:"
wc -l results/metadata/fastq_manifest.tsv

echo
echo "============================================================"
echo "BUILD VALIDATION MANIFEST"
echo "============================================================"

"$PYTHON" \
  workflow/scripts/build_validation_manifest.py \
  --manifest results/metadata/fastq_manifest.tsv \
  --validation-runs config/validation_runs.tsv \
  --output results/metadata/validation_fastq_manifest.tsv \
  > "$SAVE/build_validation_manifest.log" 2>&1

test -s results/metadata/validation_fastq_manifest.tsv

cp \
  results/metadata/validation_fastq_manifest.tsv \
  "$SAVE/validation_fastq_manifest.tsv"

if [[ -f results/metadata/validation_fastq_manifest.tsv.selection.tsv ]]; then
    cp \
      results/metadata/validation_fastq_manifest.tsv.selection.tsv \
      "$SAVE/"
fi

echo
echo "============================================================"
echo "VALIDATE TARGET CHO RUNS"
echo "============================================================"

"$PYTHON" - \
  results/metadata/validation_fastq_manifest.tsv \
  "$SAVE/run_structure_validation.tsv" <<'PY'
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
    print("  rows:", len(rr))
    print("  fastq_structure:", sorted(structures))
    print("  fastq_role:", sorted(roles))
    print("  status:", "PASS" if valid else "FAIL")

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
      > SHA256SUMS.txt

    sha256sum -c SHA256SUMS.txt
)

echo
echo "============================================================"
echo "SOURCE REPOSITORY FINAL CHECK"
echo "============================================================"

cd "$SOURCE"

STATUS_AFTER="$(git status --short)"
printf '%s\n' "$STATUS_AFTER"

[[ "$STATUS_AFTER" == "$EXPECTED" ]] || {
    echo "ERROR: source repository changed during validation" >&2
    exit 30
}

printf '%s\n' \
  "ATAC manifest validation PASS: SRR12774931 and SRR12774932 resolve as SINGLE." \
  > "$SAVE/complete.ok"

echo
echo "============================================================"
echo "COMPLETE"
echo "============================================================"

cat "$SAVE/complete.ok"

echo
echo "RESULTS:"
echo "$SAVE"
