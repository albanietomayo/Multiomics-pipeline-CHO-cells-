#!/usr/bin/env python3
"""Freeze the current eligible ATAC workflow and optionally submit one Vera job.

Usage: python3 launch_atac_bulk.py --submit
Without --submit: prepare only. No edits are made to the source repository.
"""
import argparse
import csv
import datetime
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

sys.dont_write_bytecode = True
RUNS = ('SRR12774931','SRR12774932')
EXPECTED = {
 'workflow/rules/atacseq.smk': 'ecfa6aca8bb1a51c934b1b157a322af0ae553d02ad1eceb3adc5e4dd7f8e56d3',
 'workflow/scripts/build_validation_manifest.py': 'b6a56d1338635f08d88f4a786250bba63d744289f40a12a5520c7f6252c7a756',
 'workflow/scripts/selection_gate.py': '066fc79bfda7877f2640b8ba23fd47412f61ca5a793979b72c8a48c1c82a8fbf',
 'workflow/scripts/calculate_atac_tss_enrichment.py': 'bb21fe73add62b6ee11595e36ae71eaf2329aeb4cefdf72b05ef73c63f3783b5',
}
WORKER = '#!/usr/bin/env bash\n#SBATCH --account=C3SE2026-1-42\n#SBATCH --partition=vera\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=64\n#SBATCH --mem=96G\n#SBATCH --time=24:00:00\n#SBATCH --job-name=atac_bulk_CHO_TSS\n\nset -euo pipefail\nexport PYTHONDONTWRITEBYTECODE=1\nLAUNCH="${ATAC_LAUNCH_DIR:?ATAC_LAUNCH_DIR is required}"\nJOB="${SLURM_JOB_ID:?SLURM_JOB_ID is required}"\nSCRATCH="${TMPDIR:?TMPDIR is required}"\nWORK="$SCRATCH/atac_bulk_${JOB}"\nSAVE="$LAUNCH/job_${JOB}"\nTHREADS=16\nMEM_MB=90000\n[[ ! -e "$WORK" && ! -e "$SAVE" ]] || { echo \'Work/save path already exists\' >&2; exit 2; }\nmkdir -p "$WORK" "$SAVE"\n\nQUOTA_ROOT=/cephyr/users/mayoa\ncheck_capacity() {\n    local required="$1" stats used entries\n    stats=$(cat "$QUOTA_ROOT") || return 1\n    used=$(awk \'$1=="rbytes:" {print $2}\' <<< "$stats")\n    entries=$(awk \'$1=="rentries:" {print $2}\' <<< "$stats")\n    [[ "$used" =~ ^[0-9]+$ && "$entries" =~ ^[0-9]+$ ]] || return 1\n    printf \'home_bytes=%s required_bytes=%s home_entries=%s\\n\' "$used" "$required" "$entries"\n    (( used + required + 3*1024*1024*1024 <= 30*1024*1024*1024 && entries+500 < 60000 ))\n}\n\nfinish() {\n    local result=$? diagnostics_status=0\n    trap - EXIT\n    set +e\n    local paths=() item\n    for item in .snakemake/log results/metadata results/atacseq/qc results/preprocessing/fastp/reports logs; do\n        [[ ! -e "$WORK/$item" ]] || paths+=("$item")\n    done\n    if ((${#paths[@]})); then\n        tar -C "$WORK" -czf "$SAVE/diagnostics.tar.gz" -- "${paths[@]}"\n        diagnostics_status=$?\n    fi\n    if ((result==0 && diagnostics_status!=0)); then result=$diagnostics_status; fi\n    if ((result==0)) && [[ ! -f "$SAVE/all_runs_verified.ok" ]]; then result=9; fi\n    printf \'exit_status\\t%s\\nfinished_at\\t%s\\n\' "$result" "$(date --iso-8601=seconds)" > "$SAVE/job_status.tsv"\n    if ((result==0)); then\n        printf \'Both eligible single-end runs completed; BAM and TSS outputs verified.\\n\' > "$SAVE/complete.ok" || result=9\n    fi\n    printf \'ATAC job exit=%s; results=%s\\n\' "$result" "$SAVE"\n    exit "$result"\n}\ntrap finish EXIT\ntrap \'exit 143\' TERM\ntrap \'exit 130\' INT\n\ncheck_capacity 0 > "$SAVE/home_capacity_at_start.txt" || { echo \'Insufficient HOME margin\' >&2; exit 5; }\nscontrol show job "$JOB" > "$SAVE/slurm_job_at_start.txt"\ndf -h "$SCRATCH" > "$SAVE/scratch_at_start.txt"\ncp "$LAUNCH/worker.sbatch" "$SAVE/submitted_worker.sbatch"\n(cd "$LAUNCH" && sha256sum --check snapshot.sha256) > "$SAVE/snapshot_verification.txt"\ntar -xzf "$LAUNCH/snapshot.tar.gz" -C "$WORK"\ncd "$WORK"\n\nmodule purge\nmodule load Miniforge3/24.1.2-0\nexport CONDA_PKGS_DIRS="$SCRATCH/atac-conda-pkgs"\nexport XDG_CACHE_HOME="$SCRATCH/atac-cache"\nmkdir -p "$CONDA_PKGS_DIRS" "$XDG_CACHE_HOME"\nMAIN_ENV="$SCRATCH/atac-main-env"\nCONDA_ENVS="$SCRATCH/atac-rule-envs"\nconda env create --prefix "$MAIN_ENV" --file environment.yml > "$SAVE/main_environment_creation.log" 2>&1\nexport PATH="$MAIN_ENV/bin:$PATH"\nSNAKEMAKE="$MAIN_ENV/bin/snakemake"\n"$SNAKEMAKE" --version > "$SAVE/snakemake_version.txt"\n\npython workflow/scripts/selection_gate.py --repo "$WORK" assert SRR12774931 SRR12774932 --omics ATAC-seq\npython workflow/scripts/build_fastq_manifest.py --samples config/samples.tsv --output results/metadata/fastq_manifest.tsv > "$SAVE/build_manifest.log" 2>&1\npython workflow/scripts/build_validation_manifest.py --manifest results/metadata/fastq_manifest.tsv --validation-runs config/validation_runs.tsv --output results/metadata/validation_fastq_manifest.tsv >> "$SAVE/build_manifest.log" 2>&1\npython atac_job_checks.py manifest > "$SAVE/selected_inputs.tsv"\ncp config/validation_runs.tsv "$SAVE/validation_runs_for_job.tsv"\n\n# Same files, same rules, only explicit ATAC targets. One run at a time.\nfor RUN in SRR12774931 SRR12774932; do\n    RUN_SAVE="$SAVE/runs/$RUN"\n    mkdir -p "$RUN_SAVE"\n    BAM="results/atacseq/filtered_bam/$RUN/filtered.bam"\n    QC="results/atacseq/qc/$RUN"\n    TARGETS=("$BAM" "${BAM}.bai" "$QC/tss_profile.tsv" "$QC/tss_enrichment_summary.tsv")\n    python workflow/scripts/selection_gate.py --repo "$WORK" assert "$RUN" --omics ATAC-seq\n    printf \'START %s %s\\n\' "$RUN" "$(date --iso-8601=seconds)"\n    SM_ARGS=(--snakefile Snakefile.atac_bulk --use-conda --conda-prefix "$CONDA_ENVS" --cores "$THREADS" --resources "mem_mb=$MEM_MB" --default-resources "tmpdir=$SCRATCH" --rerun-incomplete --latency-wait 60 -p)\n    "$SNAKEMAKE" "${TARGETS[@]}" "${SM_ARGS[@]}" -n > "$RUN_SAVE/dry_run.log" 2>&1\n    set +e\n    /usr/bin/time -v -o "$RUN_SAVE/time_verbose.txt" "$SNAKEMAKE" "${TARGETS[@]}" "${SM_ARGS[@]}" > "$RUN_SAVE/snakemake_execution.log" 2>&1\n    WORKFLOW_STATUS=$?\n    set -e\n    printf \'workflow_exit_status\\t%s\\n\' "$WORKFLOW_STATUS" > "$RUN_SAVE/workflow_status.tsv"\n    ((WORKFLOW_STATUS==0)) || exit "$WORKFLOW_STATUS"\n\n    BOWTIE2=$(find "$CONDA_ENVS" -path \'*/bin/bowtie2\' -print -quit)\n    [[ -n "$BOWTIE2" ]] || { echo \'Bowtie2 environment not found\' >&2; exit 3; }\n    SAMTOOLS="$(dirname "$BOWTIE2")/samtools"\n    [[ -x "$SAMTOOLS" ]] || exit 3\n    "$SAMTOOLS" --version > "$RUN_SAVE/samtools_version.txt"\n    "$BOWTIE2" --version > "$RUN_SAVE/bowtie2_version.txt"\n    "$SAMTOOLS" quickcheck -v "$BAM"\n    "$SAMTOOLS" flagstat "$BAM" > "$RUN_SAVE/filtered.flagstat.txt"\n    "$SAMTOOLS" idxstats "$BAM" > "$RUN_SAVE/filtered.idxstats.tsv"\n    TOTAL=$("$SAMTOOLS" view -c "$BAM")\n    VALID=$("$SAMTOOLS" view -c -q 30 -F 3844 "$BAM")\n    MITO=$(awk \'$1=="NC_007936.1" {n+=$3} END {print n+0}\' "$RUN_SAVE/filtered.idxstats.tsv")\n    printf \'metric\\tvalue\\ntotal_records\\t%s\\npassing_filter_records\\t%s\\nmitochondrial_records\\t%s\\n\' "$TOTAL" "$VALID" "$MITO" > "$RUN_SAVE/bam_validation.tsv"\n    ((TOTAL>0 && TOTAL==VALID && MITO==0)) || { echo \'Filtered BAM validation failed\' >&2; exit 4; }\n    python atac_job_checks.py tss "$RUN" > "$RUN_SAVE/tss_validation.txt"\n\n    ITEMS=("results/atacseq/filtered_bam/$RUN" "$QC" "results/atacseq/bowtie2/$RUN/bowtie2.log" "results/preprocessing/fastp/reports/$RUN")\n    BYTES=$(du -s -B1 -- "${ITEMS[@]}" | awk \'{n+=$1} END {printf "%.0f\\n", n}\')\n    check_capacity "$((BYTES + 256*1024*1024))" > "$RUN_SAVE/home_capacity_before_copy.txt" || { echo \'Not enough HOME capacity to preserve this run\' >&2; exit 5; }\n    mkdir -p "$RUN_SAVE/artifacts"\n    for ITEM in "${ITEMS[@]}"; do rsync -aR "$ITEM" "$RUN_SAVE/artifacts/"; done\n    # Hash every persisted result against its scratch source, not just the BAM.\n    find "${ITEMS[@]}" -type f -print0 | sort -z | xargs -0 sha256sum > "$RUN_SAVE/SHA256SUMS.txt"\n    (cd "$RUN_SAVE/artifacts" && sha256sum --check "$RUN_SAVE/SHA256SUMS.txt") > "$RUN_SAVE/persistent_checksum_verification.txt"\n    printf \'BAM and TSS artifacts verified; analytical suitability still requires review.\\n\' > "$RUN_SAVE/complete.ok"\n    cat "$QC/tss_enrichment_summary.tsv"\n    du -s -B1 "$SCRATCH" > "$RUN_SAVE/scratch_usage_bytes.txt"\n    printf \'DONE %s %s\\n\' "$RUN" "$(date --iso-8601=seconds)"\n\n    # Delete only this job\'s scratch products, after verifying persistent copies.\n    python atac_job_checks.py clean "$RUN" "$WORK" "$SCRATCH"\ndone\n\nmkdir -p "$SAVE/reference_provenance"\nfind resources/reference -type f \\( -name \'*metadata*.tsv\' -o -name SHA256SUMS.txt -o -name tss_summary.tsv -o -name tss.bed \\) -print0 |\n    while IFS= read -r -d \'\' ITEM; do rsync -aR "$ITEM" "$SAVE/reference_provenance/"; done\nfind resources/reference -type f ! -name .snakemake_timestamp -print0 | sort -z | xargs -0 sha256sum > "$SAVE/reference_SHA256SUMS.txt"\ndf -h "$SCRATCH" > "$SAVE/scratch_at_end.txt"\nprintf \'Both runs passed technical checks and persistence verification.\\n\' > "$SAVE/all_runs_verified.ok"\n'
CHECKS = "import csv\nimport math\nimport shutil\nimport sys\nfrom pathlib import Path\n\nRUNS = {'SRR12774931', 'SRR12774932'}\n\ndef read(path):\n    with open(path, newline='') as f:\n        return list(csv.DictReader(f, delimiter='\\t'))\n\ndef manifest():\n    rows = read('results/metadata/validation_fastq_manifest.tsv')\n    if len(rows) != 2 or {r['run_accession'] for r in rows} != RUNS:\n        raise ValueError('Expected exactly one FASTQ for each of the two CHO runs')\n    print('run_accession\\tfastq_bytes\\tfastq_structure\\tfastq_role')\n    for r in rows:\n        if r['omics'] != 'ATAC-seq' or any(r[k] != 'SINGLE' for k in ('processing_layout', 'fastq_structure', 'fastq_role')):\n            raise ValueError('Unexpected modality or FASTQ structure: ' + r['run_accession'])\n        if int(r['fastq_bytes']) <= 0:\n            raise ValueError('Invalid FASTQ size')\n        print('\\t'.join(r[k] for k in ('run_accession','fastq_bytes','fastq_structure','fastq_role')))\n\ndef tss(run):\n    if run not in RUNS:\n        raise ValueError('Unexpected run')\n    qc = Path('results/atacseq/qc') / run\n    for name in ('tss_enrichment_summary.tsv','tss_profile.tsv'):\n        rows = read(qc / name)\n        if not rows:\n            raise ValueError('Empty TSS output: ' + name)\n        for row in rows:\n            for key, value in row.items():\n                if value is None or not value.strip():\n                    raise ValueError('Missing field in TSS output: ' + key)\n                try:\n                    number = float(value)\n                except ValueError:\n                    continue\n                if not math.isfinite(number):\n                    raise ValueError('Nonfinite TSS value: ' + key)\n    print('TSS files are nonempty and contain no missing/nonfinite numeric values; no biological threshold applied.')\n\ndef clean(run, work, scratch):\n    root = Path(work).resolve()\n    tmp = Path(scratch).resolve()\n    if run not in RUNS or root.parent != tmp or Path.cwd().resolve() != root or not root.name.startswith('atac_bulk_'):\n        raise ValueError('Scratch cleanup scope mismatch')\n    for relative in ('data/raw/fastq','results/preprocessing/fastp','results/atacseq/bowtie2','results/atacseq/filtered_bam'):\n        path = root / relative / run\n        if path.is_symlink() or root not in path.resolve().parents:\n            raise ValueError('Unsafe scratch path')\n        if path.exists():\n            shutil.rmtree(path)\n\nif __name__ == '__main__':\n    if sys.argv[1] == 'manifest': manifest()\n    elif sys.argv[1] == 'tss': tss(sys.argv[2])\n    elif sys.argv[1] == 'clean': clean(*sys.argv[2:])\n    else: raise ValueError('Unknown operation')\n"
SNAKEFILE = '''configfile: "config/config.yaml"
include: "workflow/rules/reference.smk"
include: "workflow/rules/common.smk"
include: "workflow/rules/atacseq.smk"
'''

def capture(args, cwd=None):
    return subprocess.check_output(args, cwd=cwd, universal_newlines=True)

def validate_source(repo):
    for name, expected in EXPECTED.items():
        if hashlib.sha256((repo/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed source has changed: ' + name)
    # Stdlib-only selection gate; never install a new login-node environment.
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    subprocess.check_call([sys.executable, str(repo/'workflow/scripts/selection_gate.py'),
        '--repo', str(repo), 'assert', *RUNS, '--omics', 'ATAC-seq'], env=env)
    with (repo/'config/samples.tsv').open(newline='') as f:
        rows = [r for r in csv.DictReader(f, delimiter='\t') if r['run_accession'] in RUNS]
    if len(rows)!=2 or {r['run_accession'] for r in rows} != set(RUNS):
        raise ValueError('Missing or duplicated requested run')
    for row in rows:
        if row['library_layout']!='SINGLE' or row['omics']!='ATAC-seq' or row['instrument_platform']!='ILLUMINA':
            raise ValueError('Unexpected library metadata: ' + row['run_accession'])
        if ';' in row['fastq_ftp'] or not row['fastq_ftp'].endswith('/'+row['run_accession']+'.fastq.gz'):
            raise ValueError('Unexpected FASTQ structure')
        if int(row['fastq_bytes'])<=0 or len(row['fastq_md5'])!=32:
            raise ValueError('Missing FASTQ validation data')
        print('%s: eligible SINGLE, %.3f GiB compressed' % (row['run_accession'],int(row['fastq_bytes'])/1024**3))
    return rows

def prepare(repo, output_root):
    rows = validate_source(repo)
    payload = {}
    for folder in ('config','workflow'):
        for p in sorted((repo/folder).rglob('*')):
            rel = p.relative_to(repo)
            if any(part.startswith('.') or part=='__pycache__' for part in rel.parts) or p.suffix=='.pyc':
                continue
            if p.is_symlink():
                raise ValueError('Source symlink needs review: '+str(rel))
            if p.is_file(): payload[rel.as_posix()] = p.read_bytes()
    payload['environment.yml'] = (repo/'environment.yml').read_bytes()
    payload['Snakefile.original'] = (repo/'Snakefile').read_bytes()
    payload['config/validation_runs.original.tsv'] = payload['config/validation_runs.tsv']
    output = io.StringIO()
    writer = csv.writer(output, delimiter='\t', lineterminator='\n')
    writer.writerow(['study_accession','run_accession','omics','library_layout','processing_layout','fastq_structure'])
    for r in sorted(rows, key=lambda x:x['run_accession']):
        writer.writerow([r['study_accession'],r['run_accession'],'ATAC-seq','SINGLE','SINGLE','SINGLE'])
    payload['config/validation_runs.tsv'] = output.getvalue().encode()
    payload['Snakefile.atac_bulk'] = SNAKEFILE.encode()
    payload['atac_job_checks.py'] = CHECKS.encode()
    source_hashes = {n:hashlib.sha256(b).hexdigest() for n,b in sorted(payload.items())}
    provenance = {'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'source_repo':str(repo), 'source_commit':capture(['git','rev-parse','HEAD'], repo).strip(),
        'source_status':capture(['git','status','--short'],repo),
        'run_accessions':list(RUNS), 'snapshot_file_hashes':source_hashes,
        'scope':'Single-end ATAC through TSS. Original metadata acquisition and RNA rules not included in job entry point.'}
    output_root.mkdir(parents=True, exist_ok=True)
    prefix = 'launch_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'
    launch = Path(tempfile.mkdtemp(prefix=prefix, dir=str(output_root)))
    with tarfile.open(str(launch/'snapshot.tar.gz'), 'w:gz') as archive:
        for name, data in sorted(payload.items()):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mode = 0o644
            archive.addfile(info, io.BytesIO(data))
    checksum = hashlib.sha256((launch/'snapshot.tar.gz').read_bytes()).hexdigest()
    (launch/'snapshot.sha256').write_text(checksum+'  snapshot.tar.gz\n')
    (launch/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (launch/'worker.sbatch').write_text(WORKER)
    subprocess.check_call(['bash','-n',str(launch/'worker.sbatch')])
    print('PREPARED: '+str(launch))
    return launch

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path('/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_git'))
    parser.add_argument('--output-root',type=Path,default=Path('/cephyr/users/mayoa/Vera/atacseq_bulk_validation_results'))
    parser.add_argument('--submit',action='store_true')
    args = parser.parse_args()
    repo=args.repo.resolve()
    root=args.output_root.resolve()
    if args.submit:
        active=capture(['squeue','-h','-u',os.environ['USER'],'-o','%j']).splitlines()
        if any(name.strip()=='atac_bulk_CHO_TSS' for name in active):
            raise ValueError('ATAC bulk job already active; no duplicate submitted')
        if list(root.glob('launch_*/job_*/complete.ok')):
            raise ValueError('A completed ATAC bulk job already exists; review its results before repeating')
    launch = prepare(repo,root)
    if args.submit:
        cmd=['sbatch','--parsable','--output='+str(launch/'slurm_%j.out'),
             '--export=ALL,ATAC_LAUNCH_DIR='+str(launch), str(launch/'worker.sbatch')]
        job=capture(cmd).strip()
        (launch/'submitted_job.txt').write_text(job+'\n')
        print('SUBMITTED_JOB: '+job)
        print('LOG: '+str(launch/('slurm_'+job.split(';')[0]+'.out')))
    else:
        print('PREPARE ONLY: no job submitted. Use --submit to prepare a fresh snapshot and submit.')

if __name__=='__main__':
    main()
