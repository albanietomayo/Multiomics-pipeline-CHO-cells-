#!/usr/bin/env python3
"""Standalone Snakemake DAG and synthetic filter checks; needs Snakemake/pysam/pandas.
Samtools CLI calls in fixtures use pysam samtools dispatch, not a mock filter.
No real FASTQ, alignment, or biological acceptance is tested.
"""
import os,subprocess,shutil,tempfile,csv
from pathlib import Path
repo=Path(__file__).resolve().parents[2]
work=Path(tempfile.mkdtemp(prefix='atac-integrated-dag-'))/'repo';shutil.copytree(repo,work,ignore=shutil.ignore_patterns('__pycache__','.snakemake'))
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
logs=work.parent/'logs';logs.mkdir(exist_ok=True)
def call(args,name):
 p=subprocess.run(args,cwd=work,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 (logs/(name+'.log')).write_text(p.stdout)
 if p.returncode: print(p.stdout[-5000:]);raise RuntimeError(name)
 return p.stdout
call(['python3','workflow/scripts/build_fastq_manifest.py','--samples','config/samples.tsv','--output','results/metadata/fastq_manifest.tsv'],'manifest')
call(['python3','workflow/scripts/build_validation_manifest.py','--manifest','results/metadata/fastq_manifest.tsv','--validation-runs','config/atacseq_validation_runs.tsv','--output','results/metadata/atacseq_validation_fastq_manifest.tsv'],'atac_manifest')
args=['python3','-m','snakemake','-s','Snakefile.atac','--cores','4','--default-resources','tmpdir=/tmp','-n','-p']
s=call(args+['atacseq_all'],'fresh_dag')
assert 'atacseq_validate_outputs' in s and 'atacseq_filter_bam' in s and 'SRR12774932' in s
assert 'rule metadata:' not in s and 'rnaseq' not in s.lower() and 'SRR12774934' not in s
# Existing output audit must schedule only the two validations and aggregation.
existing=work/'fixture_existing'
for run in ['SRR12774931','SRR12774932']:
 for suffix in [f'filtered_bam/{run}/filtered.bam',f'filtered_bam/{run}/filtered.bam.bai',f'qc/{run}/tss_profile.tsv',f'qc/{run}/tss_enrichment_summary.tsv']:
  p=existing/'runs'/run/'artifacts/results/atacseq'/suffix;p.parent.mkdir(parents=True,exist_ok=True);p.touch()
s=call(args+['atacseq_existing_all','--config',f'atacseq_existing_job={existing}'],'existing_dag')
assert 'rule atacseq_validate_existing:' in s and 'rule atacseq_align_run:' not in s and 'rule atacseq_tss_enrichment:' not in s
# Root Snakefile uses the same rules when the ATAC configuration is supplied.
s=call(['python3','-m','snakemake','-s','Snakefile','atacseq_all','--configfile','config/atacseq_bulk.yaml','--cores','4','-n'],'main_entry_dag')
assert 'atacseq_validate_outputs' in s
# Test actual SE and PE filtering shells against synthetic SAM/BAM with pysam's samtools.
binpath=work/'testbin';binpath.mkdir()
wrapper=binpath/'samtools';wrapper.write_text('''#!/usr/bin/env python3
import sys,pysam
try:
 binary=sys.argv[1]=='view' and any(x in sys.argv[2:] for x in ['-b','-u'])
 r=getattr(pysam.samtools,sys.argv[1])(*sys.argv[2:],catch_stdout=not binary)
 if isinstance(r,bytes):sys.stdout.buffer.write(r)
 elif isinstance(r,str):sys.stdout.write(r)
except Exception as e:
 print(e,file=sys.stderr);sys.exit(1)
''');wrapper.chmod(0o755)
env['PATH']=str(binpath)+os.pathsep+env['PATH']
import sys
import pysam
manifest=work/'results/metadata/atacseq_validation_fastq_manifest.tsv'
with manifest.open() as f: reader=csv.DictReader(f,delimiter='\t');columns=reader.fieldnames;rows=list(reader)
run='SRR12774931'
def createbam(layout):
 p=work/f'results/atacseq/bowtie2/{run}/dupmarked.sorted.bam';p.parent.mkdir(parents=True,exist_ok=True)
 h={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'nuclear','LN':10000},{'SN':'NC_007936.1','LN':10000}]}
 if layout=='SINGLE': records=[('good',0,100,0,40,0),('dup',1024,200,0,40,0),('low',0,300,0,10,0),('mito',0,100,0,40,1)]
 else:records=[('good',99,100,200,40,0),('good',147,200,100,40,0),('orphan',99,300,400,40,0),('orphan',147,400,300,10,0)]
 with pysam.AlignmentFile(str(p),'wb',header=h) as f:
  for name,flag,pos,mate,mq,chrom in records:
   a=pysam.AlignedSegment();a.query_name=name;a.query_sequence='A'*30;a.flag=flag;a.reference_id=chrom;a.reference_start=pos;a.mapping_quality=mq;a.cigarstring='30M'
   if flag&1:a.next_reference_id=chrom;a.next_reference_start=mate
   f.write(a)
 pysam.index('-o',str(p.with_suffix('.bai')),str(p))
 return p
bamtarget=f'results/atacseq/filtered_bam/{run}/filtered.bam'
for layout,expected in [('SINGLE',1),('PAIRED',2)]:
 if layout=='PAIRED':
  row=next(r for r in rows if r['run_accession']==run);new=[]
  for role in ['R1','R2']:
   q=dict(row,fastq_role=role,fastq_structure='PAIRED',processing_layout='PAIRED',filename=run+'_'+role+'.fastq.gz');new.append(q)
  with manifest.open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(new+[r for r in rows if r['run_accession']!=run])
 createbam(layout)
 call(['python3','-m','snakemake','-s','Snakefile.atac',bamtarget,'--allowed-rules','atacseq_filter_bam','--cores','2','--default-resources','tmpdir=/tmp','--forcerun','atacseq_filter_bam','-p'],'filter_'+layout)
 with pysam.AlignmentFile(str(work/bamtarget),'rb') as b:
  r=list(b.fetch(until_eof=True));assert len(r)==expected and all(x.query_name=='good' for x in r)
print('PASS: fresh DAG, main entry DAG, existing-only DAG, real synthetic SE/PE filter execution')
print(work)
