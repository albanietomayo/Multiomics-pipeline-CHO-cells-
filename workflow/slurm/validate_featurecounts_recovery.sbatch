#!/usr/bin/env bash
#SBATCH --account=C3SE2026-1-42
#SBATCH --partition=vera
#SBATCH --job-name=validate_fc_recovery
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=/cephyr/users/mayoa/Vera/validate_fc_recovery_%j.out
set -euo pipefail

PROJECT=/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_git
SOURCE=/cephyr/users/mayoa/Vera/hpc_rnaseq_diagnostics/DRR091439/batches/SMALL_007/job_10270798/diagnostic
SAVE=/cephyr/users/mayoa/Vera/hpc_rnaseq_diagnostics/DRR091439/automation_validation_${SLURM_JOB_ID}
WORK="$TMPDIR/rnaseq_recovery_validation"
mkdir "$SAVE"
mkdir -p "$WORK/workflow/scripts" "$WORK/workflow/envs" "$WORK/input/star/DRR091439/SINGLE" "$WORK/input/strand/DRR091439"
export PROJECT SOURCE SAVE WORK

preserve() {
    result=$?
    set +e
    if [[ -d "$WORK/results" ]]; then cp -a "$WORK/results" "$SAVE/"; fi
    if [[ -d "$WORK/.snakemake/log" ]]; then cp -a "$WORK/.snakemake/log" "$SAVE/snakemake_logs"; fi
    exit "$result"
}
trap preserve EXIT

module purge
module load Miniforge3/24.1.2-0
export CONDA_PKGS_DIRS="$TMPDIR/conda-pkgs"
conda create -y --prefix "$TMPDIR/validation-main" \
    --override-channels --strict-channel-priority -c conda-forge -c bioconda \
    python=3.13 snakemake=9.25.2 "conda>=24.7.1" > "$SAVE/main_environment.log" 2>&1
export PATH="$TMPDIR/validation-main/bin:$PATH"

cp "$PROJECT/workflow/scripts/run_featurecounts.py" "$PROJECT/workflow/scripts/featurecounts_recovery.py" "$WORK/workflow/scripts/"
cp "$PROJECT/workflow/envs/rnaseq.yaml" "$PROJECT/workflow/envs/featurecounts_fallback.yaml" "$WORK/workflow/envs/"
ln -s "$SOURCE/input.bam" "$WORK/input/star/DRR091439/SINGLE/Aligned.sortedByCoord.out.bam"
ln -s "$SOURCE/resolved_strandedness.tsv" "$WORK/input/strand/DRR091439/resolved_strandedness.tsv"
ln -s "$SOURCE/annotation.gtf" "$WORK/input/annotation.gtf"
sha256sum --check "$SOURCE/input_checksums.txt" > "$SAVE/input_check.log"

python - <<'PY'
import os
from pathlib import Path
project, work, save = [Path(os.environ[x]) for x in ('PROJECT', 'WORK', 'SAVE')]
rules = (project / 'workflow/rules/rnaseq.smk').read_text()
rules = rules[rules.index('rule featurecounts_primary:'):rules.index('rule rnaseq_featurecounts_validation:')]
# Keep the rule file at its original depth so relative conda paths resolve identically.
(work / 'workflow/rules').mkdir()
(work / 'workflow/rules/recovery.smk').write_text(rules)
config = {'rnaseq': {'star_align_dir': 'input/star', 'strandedness_dir': 'input/strand', 'counts_dir': 'results/counts'}, 'reference': {'gtf': 'input/annotation.gtf'}}
(work / 'Snakefile').write_text('config = ' + repr(config) + '\ninclude: "workflow/rules/recovery.smk"\n')
(save / 'tested_rules.smk').write_text(rules)
PY

cd "$WORK"
snakemake results/counts/DRR091439/SINGLE/featurecounts.txt \
    --use-conda --conda-prefix "$TMPDIR/validation-conda" \
    --cores 4 --rerun-incomplete -p > "$SAVE/workflow.log" 2>&1

python - <<'PY'
import json
from pathlib import Path
p = Path('results/counts/DRR091439/SINGLE/featurecounts_provenance.json')
r = json.loads(p.read_text())
assert r['genes'] == 27775 and r['assigned'] == 8791549, r
assert r['fallback_used'], 'Primary succeeded in this environment: inspect provenance before rollout'
assert r['attempts'][0]['featurecounts_exit'] == -11, r
assert r['attempts'][1]['featurecounts_exit'] == 0, r
print('REAL SNAKEMAKE RECOVERY VALIDATED')
print(json.dumps(r, indent=2))
PY
printf 'Evidence: %s\n' "$SAVE"
