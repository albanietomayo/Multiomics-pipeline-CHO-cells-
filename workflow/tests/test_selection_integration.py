"""Offline behavioral checks for production, stale plans and matrix selection."""
import csv
import gzip
import hashlib
import json
import subprocess
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'workflow/scripts'))
import selection_gate as gate
from build_rnaseq_count_matrix import build
from inventory_selection_results import inventory

GOOD='SRR9621013'
BAD='SRR12774934'

def counts(path,acc,annotation='chr1',value=3):
    content=f'Geneid\tChr\tStart\tEnd\tStrand\tLength\t{acc}\ng1\t{annotation}\t1\t10\t+\t10\t{value}\n'.encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_bytes(gzip.compress(content) if path.suffix=='.gz' else content)
    return content

def saved(root,acc):
    d=root/acc;p=d/'gene_counts.tsv.gz';content=counts(p,acc)
    (d/'complete.ok').write_text('OK\n')
    (d/'gene_counts.tsv.sha256').write_text(hashlib.sha256(content).hexdigest()+'  gene_counts.tsv\n')
    (d/'gene_counts.tsv.gz.sha256').write_text(hashlib.sha256(p.read_bytes()).hexdigest()+'  gene_counts.tsv.gz\n')
    return d

class Integration(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def test_gate_blocks_tissues_multiome_and_protocol_pending(self):
        self.assertTrue(gate.decision(GOOD)[0])
        for a in [BAD,'SRR29929613']+[a for a,r in gate.catalogue().items() if r['selection_status']=='review_required']:
            self.assertFalse(gate.decision(a)[0],a)
    def test_unknown_and_wrong_modality(self):
        self.assertFalse(gate.decision('SRR0')[0]);self.assertFalse(gate.decision('SRR12774931',omics='RNA-seq')[0])
    def test_all_retained_rna_technical_support(self):
        rows=gate.catalogue();self.assertEqual(sum(gate.decision(a,omics='RNA-seq')[0] for a in rows),1753)
    def test_checksums_require_compressed_and_payload(self):
        d=saved(self.p,GOOD);self.assertTrue(gate.counts_complete(d))
        (d/'gene_counts.tsv.sha256').write_text('0'*64);self.assertFalse(gate.counts_complete(d))
    def test_old_batch_filters_and_reuses_verified_counts(self):
        pending=next(a for a,r in gate.catalogue().items() if r['selection_status']=='review_required')
        plan=self.p/'plan.tsv';accs=[BAD,pending,GOOD,'SRR10493862']
        gate.write_rows(plan,[dict(batch_id='MEDIUM_023',batch_order=i,run_accession=a) for i,a in enumerate(accs,1)],['batch_id','batch_order','run_accession'])
        saved(self.p/'runs',GOOD)
        ds,rs=gate.select_batch(plan,'MEDIUM_023',ROOT,self.p/'runs')
        self.assertEqual([r['status'] for r in ds],['BLOCKED_SELECTION','BLOCKED_SELECTION','SKIPPED_VALIDATED','TO_PROCESS'])
        self.assertEqual(rs,[{'run_accession':'SRR10493862'}])
        with self.assertRaises(ValueError):gate.select_batch(plan,'MEDIUM_023',ROOT,self.p/'runs','SRR0')
    def test_matrix_omits_old_excluded_and_pending_counts(self):
        pending=next(a for a,r in gate.catalogue().items() if r['selection_status']=='review_required')
        inputs=[]
        for a in [BAD,pending,GOOD]:
            p=self.p/(a+'.tsv.gz');counts(p,a);inputs.append(p)
        out=self.p/'matrix.tsv';self.assertEqual(build(out,inputs),1)
        self.assertEqual(out.read_text().splitlines()[0],'Geneid\t'+GOOD)
        self.assertIn('BLOCKED',Path(str(out)+'.selection.tsv').read_text())
    def test_matrix_duplicate_and_annotation_failures_preserve_old_output(self):
        a=self.p/'a.tsv';b=self.p/'b.tsv';counts(a,GOOD);counts(b,'SRR10493862','chr2');out=self.p/'matrix.tsv';out.write_text('old')
        with self.assertRaises(ValueError):build(out,[a,a])
        with self.assertRaises(ValueError):build(out,[a,b])
        self.assertEqual(out.read_text(),'old')
    def test_empty_eligible_set_does_not_produce_matrix(self):
        p=self.p/'bad.tsv';counts(p,BAD)
        with self.assertRaises(ValueError):build(self.p/'matrix.tsv',[p])
        self.assertFalse((self.p/'matrix.tsv').exists())
    def test_inventory_only_candidates_and_no_mutation(self):
        runs=self.p/'runs';saved(runs,GOOD);saved(runs,BAD)
        bam=runs/BAD/'raw.bam';bam.write_bytes(b'a'*(1024*1024))
        before=hashlib.sha256(bam.read_bytes()).hexdigest();out=self.p/'report'
        inventory(ROOT,runs,[],out)
        self.assertEqual(hashlib.sha256(bam.read_bytes()).hexdigest(),before)
        self.assertEqual((out/'eligible_counts.list').read_text().strip(),str((runs/GOOD/'gene_counts.tsv.gz').resolve()))
        self.assertIn(str(bam),(out/'excluded_files.tsv').read_text())
        rows=gate.read_tsv(out/'missing_or_unverified.tsv');self.assertFalse(any(r['selection_status']=='review_required' for r in rows))
    def test_rule_selectors_reject_stale_manifest(self):
        import pandas as pd
        manifest=pd.DataFrame([{'run_accession':a,'omics':o} for a,o in [(GOOD,'RNA-seq'),('SRR29929613','ATAC-seq'),('SRR12774931','ATAC-seq'),(BAD,'ATAC-seq')]])
        for name,expected in [('rnaseq',[GOOD]),('atacseq',['SRR12774931'])]:
            source=(ROOT/f'workflow/rules/{name}.smk').read_text();start=source.index(f'def load_{name}_validation_manifest');end=source.index('\ndef ',start+4)
            ns={'_load_validation_manifest_for_preprocessing':lambda:manifest,'selection_decision':gate.decision}
            exec(source[start:end],ns)
            self.assertEqual(ns[f'load_{name}_validation_manifest']()['run_accession'].tolist(),expected)
    def test_planner_writes_empty_batches_for_blocked_input(self):
        sizes=self.p/'sizes.tsv';metrics=self.p/'metrics.tsv'
        gate.write_rows(sizes,[dict(run_accession=BAD,fastq_bytes=1000)],['run_accession','fastq_bytes'])
        gate.write_rows(metrics,[dict(fastq_bytes=i*1000000000,elapsed_seconds=i*60,incremental_peak_tmpdir_bytes=i*1000000000,baseline_tmpdir_bytes=1000000000) for i in [1,2,3]],['fastq_bytes','elapsed_seconds','incremental_peak_tmpdir_bytes','baseline_tmpdir_bytes'])
        out=self.p/'plan'
        subprocess.run([sys.executable,str(ROOT/'workflow/scripts/plan_rnaseq_batches.py'),'--benchmark-metrics',str(metrics),'--run-sizes',str(sizes),'--outdir',str(out)],cwd=ROOT,check=True,capture_output=True)
        self.assertEqual(gate.read_tsv(out/'rnaseq_batches.tsv'),[])
        self.assertEqual(gate.read_tsv(out/'rnaseq_execution_plan.tsv')[0]['scheduling_status'],'BLOCKED_SELECTION')
    def test_validation_cli_filters_liver_before_missing_check(self):
        manifest=self.p/'manifest.tsv';runs=self.p/'runs.tsv';out=self.p/'selected.tsv'
        gate.write_rows(manifest,[dict(run_accession=GOOD,omics='RNA-seq',file_index='1',fastq_bytes='12',fastq_role='SINGLE')],['run_accession','omics','file_index','fastq_bytes','fastq_role'])
        gate.write_rows(runs,[{'run_accession':a} for a in [GOOD,BAD]],['run_accession'])
        subprocess.run([sys.executable,str(ROOT/'workflow/scripts/build_validation_manifest.py'),'--manifest',str(manifest),'--validation-runs',str(runs),'--output',str(out)],check=True,capture_output=True)
        self.assertEqual([r['run_accession'] for r in gate.read_tsv(out)],[GOOD])
        self.assertIn(BAD,Path(str(out)+'.selection.tsv').read_text())
    def test_historical_snapshot_cannot_reenable_omitted_run(self):
        repo=self.p/'repo'
        for rel in gate.FILES:
            target=repo/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,target)
        rows=gate.read_tsv(repo/'config/samples.tsv')
        gate.write_rows(repo/'config/samples.tsv',[r for r in rows if r['run_accession']!=GOOD],list(rows[0]))
        self.assertEqual(gate.decision(GOOD,repo),(False,'not_in_current_catalogue'))
        self.assertTrue(gate.decision(BAD,repo)[1].startswith('excluded:'))
    def test_inventory_refresh_updates_reports_and_preserves_other_files(self):
        runs=self.p/'runs';saved(runs,GOOD);out=self.p/'report';out.mkdir();note=out/'notes.txt';note.write_text('keep')
        for _ in range(2):
            subprocess.run([sys.executable,str(ROOT/'workflow/scripts/inventory_selection_results.py'),'--runs-root',str(runs),'--outdir',str(out),'--refresh'],check=True,capture_output=True)
        self.assertEqual(note.read_text(),'keep')
        self.assertTrue(gate.counts_complete(runs/GOOD))
    def test_shell_syntax(self):
        subprocess.run(['bash','-n',str(ROOT/'workflow/slurm/rnaseq_batch_worker.sbatch')],check=True)

if __name__=='__main__':unittest.main()
