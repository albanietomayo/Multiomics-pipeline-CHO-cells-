"""Small synthetic records test contracts; not biological validation of CHO PE."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'workflow/scripts'))
from validate_atac_outputs import check_bam, check_tss, atomic_json
from selection_gate import require
import pysam


class ATACContracts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.d = Path(self.tmp.name)

    def tearDown(self): self.tmp.cleanup()

    def bam(self, records):
        path = self.d / 'test.bam'
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                  'SQ': [{'SN': 'nuclear', 'LN': 10000}, {'SN': 'NC_007936.1', 'LN': 10000}]}
        with pysam.AlignmentFile(path, 'wb', header=header) as f:
            for name, flag, pos, mate, mapq, chrom in records:
                r = pysam.AlignedSegment(); r.query_name = name; r.query_sequence = 'A' * 30
                r.flag = flag; r.reference_id = chrom; r.reference_start = pos
                r.mapping_quality = mapq; r.cigarstring = '30M'
                if flag & 1: r.next_reference_id = chrom; r.next_reference_start = mate
                f.write(r)
        pysam.index(str(path)); return path

    def check(self, rows, layout):
        b = self.bam(rows)
        return check_bam(b, str(b)+'.bai', layout, 30, 'NC_007936.1', self.d)

    def test_single(self):
        self.assertEqual(self.check([('s',0,100,0,30,0)], 'SINGLE'), 1)

    def test_duplicate_lowmapq_and_mito_rejected(self):
        for flag, mapq, chrom in [(1024,40,0), (0,29,0), (0,40,1)]:
            with self.subTest(flag=flag, mapq=mapq, chrom=chrom):
                with self.assertRaises(ValueError): self.check([('s',flag,100,0,mapq,chrom)], 'SINGLE')

    def test_complete_pair(self):
        self.assertEqual(self.check([('p',99,100,200,40,0),('p',147,200,100,40,0)], 'PAIRED'), 2)

    def test_two_orphans_with_even_count(self):
        with self.assertRaises(ValueError):
            self.check([('a',99,100,200,40,0),('b',147,200,100,40,0)], 'PAIRED')

    def test_empty_and_singleton(self):
        for rows in [[], [('p',99,100,200,40,0)]]:
            with self.assertRaises(ValueError): self.check(rows,'PAIRED')

    def test_inconsistent_mates(self):
        with self.assertRaises(ValueError):
            self.check([('p',99,100,999,40,0),('p',147,200,100,40,0)],'PAIRED')

    def test_wrong_layout(self):
        with self.assertRaises(ValueError):
            self.check([('p',99,100,200,40,0),('p',147,200,100,40,0)],'SINGLE')

    def test_selection(self):
        for run in ['SRR12774931','SRR12774932']: require(run, ROOT, 'ATAC-seq')
        for run in ['SRR12774934','SRR29929613','UNKNOWN']:
            with self.assertRaises(ValueError): require(run, ROOT, 'ATAC-seq')

    def test_real_tss_and_corruptions(self):
        base = ROOT / 'benchmarks/atacseq/2026-09-16/job_10297369'
        for run, count in [('SRR12774931',7068304),('SRR12774932',5858817)]:
            q = base/run
            check_tss(q/'tss_profile.tsv', q/'tss_enrichment_summary.tsv', count,2000,10,100)
            with self.assertRaises(ValueError):
                check_tss(q/'tss_profile.tsv', q/'tss_enrichment_summary.tsv', count+1,2000,10,100)
            rows = (q/'tss_profile.tsv').read_text().splitlines()
            fields = rows[1].split('\t'); fields[-1] = 'nan'; rows[1] = '\t'.join(fields)
            bad = self.d/'bad.tsv'; bad.write_text('\n'.join(rows)+'\n')
            with self.assertRaises(ValueError): check_tss(bad,q/'tss_enrichment_summary.tsv',count,2000,10,100)
            fields[-1] = '999'; rows[1] = '\t'.join(fields); bad.write_text('\n'.join(rows)+'\n')
            with self.assertRaises(ValueError): check_tss(bad,q/'tss_enrichment_summary.tsv',count,2000,10,100)

    def test_invalid_write_preserves_previous_report(self):
        p=self.d/'report.json';p.write_text('previous\n')
        with self.assertRaises(ValueError): atomic_json(p,{'score':float('nan')})
        self.assertEqual(p.read_text(),'previous\n')


if __name__ == '__main__': unittest.main()
