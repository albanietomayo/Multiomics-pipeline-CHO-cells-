#!/usr/bin/env python3
"""Freeze the current eligible ATAC workflow and optionally submit one Vera job.

Usage: python3 workflow/scripts/launch_atac_bulk.py --submit
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

def capture(args, cwd=None):
    return subprocess.check_output(args, cwd=cwd, universal_newlines=True)

def validate_source(repo):
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
    payload['Snakefile'] = (repo/'Snakefile').read_bytes()
    payload['Snakefile.atac'] = (repo/'Snakefile.atac').read_bytes()
    with (repo/'config/atacseq_validation_runs.tsv').open(newline='') as handle:
        configured = list(csv.DictReader(handle, delimiter='\t'))
    if len(configured) != len(RUNS) or {r['run_accession'] for r in configured} != set(RUNS):
        raise ValueError('This Vera launcher supports only the documented two-run CHO cohort')
    for r in configured:
        if r['omics'] != 'ATAC-seq' or r['fastq_structure'] != 'SINGLE':
            raise ValueError('Unexpected cohort structure')
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
    (launch/'worker.sbatch').write_bytes((repo/'workflow/slurm/atacseq_bulk_validation.sbatch').read_bytes())
    subprocess.check_call(['bash','-n',str(launch/'worker.sbatch')])
    print('PREPARED: '+str(launch))
    return launch

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[2])
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
