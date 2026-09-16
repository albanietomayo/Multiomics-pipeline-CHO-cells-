"""Offline regression cases from the 2026-09-14 selection audit.

Run from the repository root:
python3 -m unittest discover -s workflow/tests -p test_cho_selection_policy.py -v
"""
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'workflow/scripts'))
from selection_policy import audit, classify, load_decisions, read_tsv

RULES = read_tsv(ROOT/'config/classification_rules.tsv')
CURATED = read_tsv(ROOT/'config/curated_studies.tsv')
DECISIONS = load_decisions(ROOT/'config/selection_decisions.tsv')
ROWS = json.loads((Path(__file__).parent/'fixtures/cho_selection_cases.json').read_text())
BY = {r['run_accession']:r for r in ROWS}
RESULT = {r['run_accession']:r for r in audit(ROWS,RULES,CURATED,DECISIONS)}

class SelectionRegression(unittest.TestCase):
    def test_non_cho_and_bulk_control(self):
        for acc in ('SRR12774934','DRR706774','SRR950107'):
            self.assertEqual(RESULT[acc]['selection_status'],'excluded')
        self.assertEqual(RESULT['SRR12774931']['selection_status'],'retained_by_rules')

    def test_cell_cloning_is_not_single_cell_sequencing(self):
        for acc in ('DRR126388','SRR29929613'):
            self.assertEqual(RESULT[acc]['protocol_status'],'single_cell')
        self.assertEqual(RESULT['ERR4184070']['protocol_status'],'bulk_confirmed')

    def test_specialized_rna_keeps_cho_identity(self):
        for acc in ('SRR3284463','SRR070303','SRR32329206'):
            self.assertEqual(RESULT[acc]['target_status'],'target')
            self.assertEqual(RESULT[acc]['selection_status'],'excluded')
        self.assertEqual(RESULT['SRR3284457']['selection_status'],'retained_by_rules')

    def test_ip_remains_pending(self):
        self.assertEqual(RESULT['SRR36663586']['target_status'],'target')
        self.assertEqual(RESULT['SRR36663586']['selection_status'],'review_required')

    def test_changed_identifiers_fail(self):
        row=dict(BY['SRR070303'],study_accession='DIFFERENT_PROJECT')
        with self.assertRaises(ValueError):
            classify(row,RULES,CURATED,DECISIONS,('RNA-seq','ATAC-seq','ChIP-seq'))

if __name__ == '__main__':
    unittest.main()
