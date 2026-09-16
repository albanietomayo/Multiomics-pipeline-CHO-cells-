"""Read-only reconciliation of eligibility and saved results. Does not delete or submit jobs."""
import argparse
import json
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from selection_gate import ROOT, catalogue, decision, counts_complete, fingerprint, write_rows

RUN = re.compile(r'^[SED]RR[0-9]+$')
LARGE = ('.bam','.sam','.cram','.fastq','.fastq.gz','.fq','.fq.gz')

def inventory(repo, runs_root, artifact_roots, outdir):
    repo=Path(repo);runs_root=Path(runs_root);outdir=Path(outdir)
    if not runs_root.is_dir():raise ValueError('Missing production runs directory: '+str(runs_root))
    roots=[runs_root]+[Path(p) for p in artifact_roots]
    if any(not p.is_dir() for p in roots):raise ValueError('Every artifact root must exist')
    outdir.mkdir(parents=True,exist_ok=False)
    cat=catalogue(repo);states=[];ready=[];files=[];visited=set()
    # Include unknown result directories so an old output cannot disappear from reconciliation.
    accessions=set(cat)|{p.name for p in runs_root.iterdir() if p.is_dir() and RUN.fullmatch(p.name)}
    for acc in sorted(accessions):
        row=cat.get(acc,{});ok,why=decision(acc,repo,'RNA-seq')
        saved=counts_complete(runs_root/acc)
        status='BLOCKED_SELECTION' if not ok else ('VERIFIED_COUNTS' if saved else 'MISSING_OR_UNVERIFIED')
        states.append(dict(run_accession=acc,omics=row.get('omics',''),selection_status=row.get('selection_status','unknown'),status=status,reason=why,verified_counts=saved))
        if ok and saved:ready.append(str((runs_root/acc/'gene_counts.tsv.gz').resolve()))
    for root in roots:
        for parent,dirs,names in os.walk(root,followlinks=False):
            dirs[:]=[d for d in dirs if not (Path(parent)/d).is_symlink()]
            for name in names:
                p=Path(parent)/name
                if p.is_symlink():continue
                resolved=str(p.resolve())
                if resolved in visited:continue
                visited.add(resolved)
                ids={part for part in p.parts if RUN.fullmatch(part)}
                if len(ids)!=1:continue
                acc=next(iter(ids));row=cat.get(acc,{})
                if row.get('selection_status')!='excluded':continue
                stat=p.stat()
                # Candidates require a separate usage review on Vera; this is not a deletion list.
                candidate=name.endswith(LARGE) and stat.st_size>=1024*1024
                files.append(dict(run_accession=acc,path=resolved,bytes=stat.st_size,mtime_ns=stat.st_mtime_ns,device=stat.st_dev,inode=stat.st_ino,candidate_for_usage_review=candidate,reason=row.get('selection_reason','')))
    write_rows(outdir/'run_status.tsv',states,['run_accession','omics','selection_status','status','reason','verified_counts'])
    write_rows(outdir/'excluded_files.tsv',files,['run_accession','path','bytes','mtime_ns','device','inode','candidate_for_usage_review','reason'])
    (outdir/'eligible_counts.list').write_text(''.join(p+'\n' for p in ready))
    write_rows(outdir/'missing_or_unverified.tsv',[r for r in states if r['status']=='MISSING_OR_UNVERIFIED'],list(states[0]))
    scheduler={'available':False}
    if shutil.which('squeue'):
        result=subprocess.run(['squeue','-u',os.environ.get('USER',''),'--noheader','-o','%i|%j|%T|%M|%R'],capture_output=True,text=True,timeout=30)
        scheduler={'available':True,'returncode':result.returncode}
        (outdir/'squeue.txt').write_text(result.stdout+result.stderr)
    if not scheduler['available']:(outdir/'squeue.txt').write_text('squeue unavailable\n')
    summary={'policy_inputs':fingerprint(repo),'runs_root':str(runs_root.resolve()),'artifact_roots':[str(p.resolve()) for p in roots],'statuses':dict(Counter(r['status'] for r in states)),'candidate_bytes':sum(r['bytes'] for r in files if r['candidate_for_usage_review']),'scheduler':scheduler,'limitations':'Missing does not mean ready to submit: reconcile active jobs first. File candidates have not been proven unused. Archives are not unpacked. No files deleted.'}
    (outdir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary['statuses']))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,default=ROOT)
    p.add_argument('--runs-root',required=True,type=Path);p.add_argument('--artifact-root',action='append',default=[],type=Path);p.add_argument('--outdir',required=True,type=Path)
    p.add_argument('--refresh',action='store_true',help='Replace generated reports in outdir; never modify result inputs')
    a=p.parse_args()
    if not a.refresh:
        inventory(a.repo,a.runs_root,a.artifact_root,a.outdir)
    else:
        a.outdir.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='selection_inventory_') as tmp:
            reports=Path(tmp)/'reports'
            inventory(a.repo,a.runs_root,a.artifact_root,reports)
            for src in reports.iterdir():
                # Copy to destination filesystem before atomic replacement.
                temp=a.outdir/(src.name+'.new')
                shutil.copyfile(src,temp);os.replace(temp,a.outdir/src.name)
if __name__=='__main__':main()
