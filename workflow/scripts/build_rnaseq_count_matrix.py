"""Build a matrix from currently eligible RNA runs; record every input decision."""
import argparse
import csv
import gzip
import hashlib
import json
import os
import tempfile
from pathlib import Path
from selection_gate import ROOT, decision, fingerprint, write_rows

ANNOTATION_COLUMNS = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']

def open_counts(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)

def read_run_counts(path, header_only=False):
    with open_counts(path) as f:
        reader=csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames or reader.fieldnames[:6]!=ANNOTATION_COLUMNS or len(reader.fieldnames)!=7:
            raise ValueError(f'Expected six annotation columns and one run column: {path}')
        acc=reader.fieldnames[6];rows={}
        if header_only:return acc,rows
        for row in reader:
            gene=row['Geneid'];count=int(row[acc])
            if not gene or gene in rows or count<0:raise ValueError(f'Invalid or duplicate gene/count: {path}: {gene}')
            rows[gene]={'annotation':[row[c] for c in ANNOTATION_COLUMNS],'count':count}
        if not rows:raise ValueError(f'Empty counts: {path}')
    return acc,rows

def build(output, inputs, repo=ROOT):
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    report=[];data=[];seen=set();hashes={}
    for path in inputs:
        path=Path(path);acc,_=read_run_counts(path,True)
        ok,why=decision(acc,repo,'RNA-seq')
        report.append(dict(path=str(path.resolve()),run_accession=acc,status='INCLUDED' if ok else 'BLOCKED',reason=why))
        if not ok:continue
        if acc in seen:raise ValueError('Duplicate run input: '+acc)
        seen.add(acc);data.append(read_run_counts(path))
        hashes[str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
    write_rows(str(output)+'.selection.tsv',report,['path','run_accession','status','reason'])
    if not data:raise ValueError('No eligible RNA-seq count inputs; matrix not written')
    first=data[0][1]
    for acc,rows in data[1:]:
        if rows.keys()!=first.keys():raise ValueError('Gene set differs: '+acc)
        for gene in first:
            if rows[gene]['annotation']!=first[gene]['annotation']:raise ValueError('Annotation differs: '+acc+': '+gene)
    # Never replace a valid old matrix with a partially written new matrix.
    fd,name=tempfile.mkstemp(dir=output.parent,prefix=output.name+'.')
    try:
        with os.fdopen(fd,'w') as f:
            w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['Geneid']+[a for a,_ in data])
            for gene in first:w.writerow([gene]+[rows[gene]['count'] for _,rows in data])
        os.replace(name,output)
    finally:
        if os.path.exists(name):os.unlink(name)
    Path(str(output)+'.provenance.json').write_text(json.dumps({'policy_inputs':fingerprint(repo),'count_inputs':hashes,'matrix_sha256':hashlib.sha256(output.read_bytes()).hexdigest()},indent=2)+'\n')
    return len(data)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);p.add_argument('counts',nargs='*',type=Path)
    p.add_argument('--repo',type=Path,default=ROOT);p.add_argument('--counts-list',type=Path)
    a=p.parse_args();inputs=a.counts
    if a.counts_list:inputs += [Path(s) for s in a.counts_list.read_text().splitlines() if s.strip()]
    if not inputs:p.error('Provide count paths or --counts-list')
    print('Matrix columns:',build(a.output,inputs,a.repo))

if __name__=='__main__':main()
