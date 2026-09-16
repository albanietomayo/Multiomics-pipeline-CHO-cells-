"""Shared offline eligibility gate for plans, workers, validation and matrices.

Historical rows supply explanations for exclusions missing from samples.tsv.
Only runs present in the current samples.tsv can be eligible. Decisions are
recomputed from versioned rules; a cached old 'target' is never sufficient.
"""
import argparse
import csv
import functools
import gzip
import hashlib
import json
from pathlib import Path
from selection_policy import audit, read_tsv, load_decisions

ROOT = Path(__file__).resolve().parents[2]
FILES = ('config/selection_audit_input_20260914.tsv', 'config/samples.tsv',
         'config/classification_rules.tsv', 'config/curated_studies.tsv',
         'config/selection_decisions.tsv', 'workflow/scripts/selection_policy.py',
         'workflow/scripts/selection_gate.py')
SUPPORTED = {'ILLUMINA', 'DNBSEQ'}

def fingerprint(repo=ROOT):
    repo=Path(repo)
    return {f:hashlib.sha256((repo/f).read_bytes()).hexdigest() for f in FILES}

@functools.lru_cache(maxsize=4)
def _load(repo, signature):
    repo=Path(repo)
    old=read_tsv(repo/FILES[0]);current=read_tsv(repo/FILES[1])
    for rows in (old,current):
        if len(rows)!=len({r['run_accession'] for r in rows}):
            raise ValueError('Duplicate run in selection source')
    merged={r['run_accession']:r for r in old}
    merged.update({r['run_accession']:r for r in current})
    rows=audit(list(merged.values()),read_tsv(repo/FILES[2]),read_tsv(repo/FILES[3]),load_decisions(repo/FILES[4]))
    present={r['run_accession'] for r in current}
    for r in rows:r['in_current_catalogue']=r['run_accession'] in present
    return {r['run_accession']:r for r in rows}

def catalogue(repo=ROOT):
    repo=Path(repo).resolve()
    sig=tuple(((repo/f).stat().st_mtime_ns,(repo/f).stat().st_size) for f in FILES)
    return _load(str(repo),sig)

def decision(run, repo=ROOT, omics=None):
    r=catalogue(repo).get(run)
    if r is None:return False,'unknown_run'
    if r['selection_status']!='retained_by_rules':
        return False,r['selection_status']+':'+r['selection_reason']
    if not r['in_current_catalogue']:return False,'not_in_current_catalogue'
    if omics and r['omics']!=omics:return False,'different_omics'
    if r['omics']=='RNA-seq' and r.get('instrument_platform','').upper() not in SUPPORTED:
        return False,'unsupported_platform'
    return True,'eligible'

def require(run, repo=ROOT, omics=None):
    ok,reason=decision(run,repo,omics)
    if not ok:raise ValueError(f'Run {run} is not eligible: {reason}')

def write_rows(path, rows, columns):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def counts_complete(run_dir):
    """Same integrity contract as the production worker, including gzip payload."""
    d=Path(run_dir)
    try:
        if not all((d/n).is_file() for n in ('complete.ok','gene_counts.tsv.gz','gene_counts.tsv.gz.sha256','gene_counts.tsv.sha256')):return False
        expected=(d/'gene_counts.tsv.gz.sha256').read_text().split()[0]
        if hashlib.sha256((d/'gene_counts.tsv.gz').read_bytes()).hexdigest()!=expected:return False
        h=hashlib.sha256()
        with gzip.open(d/'gene_counts.tsv.gz','rb') as f:
            for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
        return h.hexdigest()==(d/'gene_counts.tsv.sha256').read_text().split()[0]
    except (OSError,ValueError,IndexError,EOFError):return False

def select_batch(plan, batch, repo, runs_root, pilot=None):
    rows=[r for r in read_tsv(plan) if r['batch_id']==batch]
    if not rows:raise ValueError('Batch not found: '+batch)
    if len(rows)!=len({r['run_accession'] for r in rows}):raise ValueError('Duplicate runs in batch')
    if pilot:
        if pilot not in {r['run_accession'] for r in rows}:raise ValueError('Pilot is not in this batch')
        rows=[r for r in rows if r['run_accession']==pilot]
    decisions=[];selected=[]
    for r in sorted(rows,key=lambda r:int(r['batch_order'])):
        acc=r['run_accession'];ok,why=decision(acc,repo,'RNA-seq')
        state='BLOCKED_SELECTION'
        if ok:
            state='SKIPPED_VALIDATED' if counts_complete(Path(runs_root)/acc) else 'TO_PROCESS'
        decisions.append(dict(run_accession=acc,status=state,reason=why))
        if state=='TO_PROCESS':selected.append(dict(run_accession=acc))
    return decisions,selected

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,default=ROOT)
    sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('assert');a.add_argument('runs',nargs='+');a.add_argument('--omics')
    b=sub.add_parser('batch');b.add_argument('--plan',required=True,type=Path);b.add_argument('--batch',required=True)
    b.add_argument('--runs-root',required=True,type=Path);b.add_argument('--outdir',required=True,type=Path);b.add_argument('--pilot')
    a=p.parse_args()
    if a.command=='assert':
        for acc in a.runs:require(acc,a.repo,a.omics)
    else:
        if a.outdir.exists():raise ValueError('Use a new selection output directory')
        ds,rs=select_batch(a.plan,a.batch,a.repo,a.runs_root,a.pilot)
        write_rows(a.outdir/'decisions.tsv',ds,['run_accession','status','reason'])
        write_rows(a.outdir/'runs.tsv',rs,['run_accession'])
        (a.outdir/'provenance.json').write_text(json.dumps({'policy_inputs':fingerprint(a.repo),
            'plan_sha256':hashlib.sha256(a.plan.read_bytes()).hexdigest(),'batch_id':a.batch},indent=2)+'\n')
        print(f'{len(rs)} runs need processing; {len(ds)-len(rs)} blocked or already validated')

if __name__=='__main__':main()
