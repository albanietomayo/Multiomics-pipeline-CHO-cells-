#!/usr/bin/env bash
#SBATCH --job-name=chipseq_reference_provision
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00

# Historical evidence: bowtie2-build used approximately 6.7 GB RSS and 15
# minutes. This one-time envelope provides >2x memory and 4x build walltime.

set -euo pipefail

: "${TMPDIR:?Compute-node TMPDIR is required}"
: "${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
: "${CHIP_PROVISION_EXPECTED_COMMIT:?A reviewed source commit is required}"

CHIP_ROOT="$PWD"
CHIP_SUPPORT="workflow/scripts/chipseq_benchmark.py"
CHIP_PROVISION="workflow/scripts/chipseq_reference_provision.py"
CHIP_SOURCE="/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_chipseq/results/chipseq/alignment/slurm/submission__niyvyr2/job_10297460/outputs/reference"
CHIP_DEST="/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_benchmark/resources/reference/CriGri-PICRH-1.0/chipseq"
CHIP_PARENT="$(dirname "$CHIP_DEST")"
CHIP_STAGE="$CHIP_PARENT/.chipseq.provisioning.${SLURM_JOB_ID}.staging"
CHIP_WORK="$TMPDIR/chipseq-reference-${SLURM_JOB_ID}"
CHIP_VALIDATE="$CHIP_WORK/validation"
CHIP_INDEX="$CHIP_WORK/bowtie2_index"
CHIP_FASTA="$CHIP_WORK/genome_plus_mt.fa"
CHIP_ENV="$CHIP_WORK/environment"
CHIP_PERSISTENT_STAGE_CREATED=0

cleanup() {
    local status=$?
    trap - EXIT
    set +e
    if [[ "$status" != 0 && "$CHIP_PERSISTENT_STAGE_CREATED" == 1 && -d "$CHIP_STAGE" ]]; then
        chmod -R u+w "$CHIP_STAGE"
        rm -rf -- "$CHIP_STAGE"
    fi
    exit "$status"
}
trap cleanup EXIT

[[ "$(git branch --show-current)" == "chipseq-incremental-generalization" ]] || {
    echo "ERROR: wrong provisioning branch" >&2
    exit 2
}
[[ "$(git rev-parse HEAD)" == "$CHIP_PROVISION_EXPECTED_COMMIT" ]] || {
    echo "ERROR: provisioning commit changed after submission" >&2
    exit 2
}
[[ -z "$(git status --short)" ]] || {
    echo "ERROR: provisioning worktree is dirty" >&2
    exit 2
}

python "$CHIP_PROVISION" validate-source --source-root "$CHIP_SOURCE"
python "$CHIP_PROVISION" validate-destination --source-root "$CHIP_SOURCE" --destination "$CHIP_DEST"

mkdir -p "$CHIP_WORK" "$CHIP_INDEX" "$CHIP_VALIDATE/bowtie2_index"
gzip -dc -- "$CHIP_SOURCE/genome_plus_mt.fa.gz" > "$CHIP_FASTA"
printf '%s  %s\n' dfd445e136c4bc9c11f9616d251b86732d1aab0cbdbf9d1ababb8093f7430e94 "$CHIP_FASTA" | sha256sum --check --strict

module purge
module load Miniforge3/24.1.2-0
export CONDA_PKGS_DIRS="$CHIP_WORK/conda-pkgs"
export XDG_CACHE_HOME="$CHIP_WORK/cache"
export PYTHONDONTWRITEBYTECODE=1
conda env create --prefix "$CHIP_ENV" --file workflow/envs/chipseq_alignment.yaml
export PATH="$CHIP_ENV/bin:$PATH"
bowtie2-build --version > "$CHIP_WORK/bowtie2_build_version.txt"
grep -F -m 1 "version 2.5.5" "$CHIP_WORK/bowtie2_build_version.txt"

bowtie2-build --large-index --threads "$SLURM_CPUS_PER_TASK" "$CHIP_FASTA" "$CHIP_INDEX/genome_plus_mt"

for suffix in 1 2 3 4 rev.1 rev.2; do
    component="$CHIP_INDEX/genome_plus_mt.${suffix}.bt2l"
    [[ -s "$component" ]] || {
        echo "ERROR: missing or empty index component: $component" >&2
        exit 2
    }
    ln "$component" "$CHIP_VALIDATE/bowtie2_index/$(basename "$component")"
done
for name in genome_plus_mt.fa.gz genome_plus_mt.fa.fai reference_provenance.json; do
    ln -s "$CHIP_SOURCE/$name" "$CHIP_VALIDATE/$name"
done

python "$CHIP_SUPPORT" create-reference-inventory --shared-reference-root "$CHIP_VALIDATE"
python "$CHIP_SUPPORT" check-reference --shared-reference-root "$CHIP_VALIDATE"
python "$CHIP_PROVISION" component-report --root "$CHIP_VALIDATE" --output "$CHIP_WORK/provisioning_report.json"
cat "$CHIP_WORK/provisioning_report.json"

# No persistent payload is created until the exact generated index and quota
# calculation above pass. The hidden sibling is removed on every failed exit.
mkdir -p "$CHIP_PARENT"
[[ ! -e "$CHIP_STAGE" && ! -L "$CHIP_STAGE" ]] || {
    echo "ERROR: staging collision: $CHIP_STAGE" >&2
    exit 2
}
mkdir "$CHIP_STAGE"
CHIP_PERSISTENT_STAGE_CREATED=1
mkdir "$CHIP_STAGE/bowtie2_index"

for name in genome_plus_mt.fa.gz genome_plus_mt.fa.fai reference_provenance.json; do
    ln "$CHIP_SOURCE/$name" "$CHIP_STAGE/$name"
done
for suffix in 1 2 3 4 rev.1 rev.2; do
    source_component="$CHIP_INDEX/genome_plus_mt.${suffix}.bt2l"
    destination_component="$CHIP_STAGE/bowtie2_index/genome_plus_mt.${suffix}.bt2l"
    cp -- "$source_component" "$destination_component.part"
    cmp -- "$source_component" "$destination_component.part"
    mv -- "$destination_component.part" "$destination_component"
done
cp -- "$CHIP_VALIDATE/reference_inventory.json" "$CHIP_STAGE/reference_inventory.json.part"
cmp -- "$CHIP_VALIDATE/reference_inventory.json" "$CHIP_STAGE/reference_inventory.json.part"
mv -- "$CHIP_STAGE/reference_inventory.json.part" "$CHIP_STAGE/reference_inventory.json"

python "$CHIP_SUPPORT" check-reference --shared-reference-root "$CHIP_STAGE"
python "$CHIP_PROVISION" verify-current-reserve

# Hard-linked validated source objects become immutable through both names.
find "$CHIP_STAGE" -type f -exec chmod 0444 {} +
chmod 0555 "$CHIP_STAGE/bowtie2_index" "$CHIP_STAGE"
[[ ! -e "$CHIP_DEST" && ! -L "$CHIP_DEST" ]] || {
    echo "ERROR: destination collision immediately before publication" >&2
    exit 2
}
mv -T -- "$CHIP_STAGE" "$CHIP_DEST"
CHIP_PERSISTENT_STAGE_CREATED=0

python "$CHIP_SUPPORT" check-reference --shared-reference-root "$CHIP_DEST"
[[ ! -e "$CHIP_DEST/genome_plus_mt.fa" ]] || {
    echo "ERROR: uncompressed FASTA was published" >&2
    exit 2
}
echo "REFERENCE_PROVISIONING=COMPLETE"
echo "CANONICAL_REFERENCE_ROOT=$CHIP_DEST"
