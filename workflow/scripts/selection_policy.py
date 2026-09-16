#!/usr/bin/env python3
"""Conservative, offline selection policy. No network or file mutations on import.

Study titles are review context, never sufficient CHO evidence. Existing curated
and structured lineage rules precede generic text matches. Conflicting decisions
at that level are held for review (no last-match-wins). Protocol is independent.
'not_flagged' and 'retained_by_rules' do NOT mean protocol/biology validated.
"""
import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

POLICY_VERSION = '2026-09-14-proposal3'
SC = re.compile(r'single[\s_-]*cell|multiome|single[\s_-]*nucle(?:us|i)', re.I)
RNA_SPECIAL = re.compile(
    r'ribo[\s_-]*seq|ribosome profil|small[\s_-]*RNA|miRNA[\s_-]*seq|'
    r'meRIP|immunoprecip|target[\s_-]*seq|amplicon|(?:^|[\s_-])(?:IP|input)(?:$|[\s_-])|'
    r'(?:PRO|GRO|RIP)[\s_-]*seq', re.I)
ASSAY_FIELDS = ('sample_title', 'experiment_title', 'library_name', 'library_selection')
PROTOCOL_FIELDS = ('library_source', 'study_title', 'sample_title',
                   'experiment_title', 'library_name', 'cell_type')
ADDED = ('previous_target_status', 'previous_target_evidence', 'previous_omics',
         'target_status', 'target_evidence', 'identity_rule_ids', 'identity_reason',
         'assay_flag_fields', 'curation_evidence',
         'study_context_rule_ids', 'protocol_status', 'protocol_flag_fields',
         'curation_references', 'selection_status', 'selection_reason', 'selection_policy')


def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def match(row, rule):
    field = rule['field']
    if field not in row:
        raise ValueError(f'Missing classification field: {field}')
    value = str(row.get(field) or '').strip()
    pattern = rule['pattern'].strip()
    if not pattern:
        raise ValueError('Empty classification pattern')
    kind = rule['match_type'].strip().lower()
    if kind == 'regex':
        return bool(re.search(pattern, value, re.I))
    if kind == 'exact':
        return value.casefold() == pattern.casefold()
    if kind == 'prefix':
        return value.casefold().startswith(pattern.casefold())
    if kind == 'contains':
        return pattern.casefold() in value.casefold()
    raise ValueError(f'Unsupported match type: {kind}')


def outcome(hits):
    values = {h['classification'] for h in hits}
    if not values.issubset({'target', 'non_target', 'unresolved'}):
        raise ValueError(f'Unknown identity classifications: {values}')
    return next(iter(values)) if len(values) == 1 else 'unresolved'


def classify(row, rules, curated, decisions, target_omics):
    out = dict(row)
    for field in ('target_status', 'target_evidence', 'omics'):
        out['previous_' + field] = row.get('previous_' + field, row.get(field, ''))
    context, text_hits, specific = [], [], []
    for r in rules:
        if not match(row, r):
            continue
        if r['field'] == 'study_title':
            context.append(r)
        elif r['stage'] == 'explicit_evidence':
            text_hits.append(r)
        elif r['stage'] in {'cell_line', 'cell_type', 'external'}:
            specific.append(r)
        else:
            raise ValueError(f"Unknown classification stage: {r['stage']}")
    curation_hits = [r for r in curated if r['study_accession'] == row['study_accession']
                     and match(row, r)]
    decisive = specific + curation_hits
    if decisive:
        identity = outcome(decisive)
        reason = 'curated_or_lineage_rules' if identity != 'unresolved' else 'conflicting_or_unresolved_specific_rules'
    elif text_hits:
        decisive = text_hits
        identity = outcome(decisive)
        reason = 'sample_level_text' if identity != 'unresolved' else 'conflicting_sample_text'
    else:
        identity = 'unresolved'
        reason = 'study_title_only' if context else 'no_sample_identity_evidence'
    refs = [r.get('reference', '') for r in decisive if r.get('reference')]
    ids = [r.get('rule_id', r.get('curation_id', '')) for r in decisive]
    flags = [f for f in PROTOCOL_FIELDS if SC.search(str(row.get(f) or ''))]
    source_flag = bool(SC.search(str(row.get('library_source') or '')))
    protocol = 'single_cell' if source_flag else 'review_required' if flags else 'not_flagged'
    assay_flags = [f for f in ASSAY_FIELDS if RNA_SPECIAL.search(str(row.get(f) or ''))] if str(row.get('library_strategy') or '').lower() == 'rna-seq' else []
    if assay_flags and protocol != 'single_cell':
        protocol = 'review_required'
    decision = decisions.get(row['run_accession'])
    if decision:
        if decision['study_accession'] != row['study_accession'] or decision['sample_accession'] != row['sample_accession']:
            raise ValueError('Curated run identity changed: ' + row['run_accession'])
        if decision['identity_status']:
            identity = decision['identity_status']
            if identity not in {'target', 'non_target', 'unresolved'}:
                raise ValueError('Invalid curated identity')
            reason = 'reviewed_run_evidence'
        if decision['protocol_status']:
            protocol = decision['protocol_status']
            if protocol not in {'bulk_confirmed', 'single_cell', 'other_assay', 'review_required'}:
                raise ValueError('Invalid curated protocol')
        ids.append(decision['decision_id'])
        refs.append(decision['reference'])
    # Keep assay identity separate from exclusion due to library architecture.
    omics = {'rna-seq': 'RNA-seq', 'atac-seq': 'ATAC-seq', 'chip-seq': 'ChIP-seq'}.get(
        str(row.get('library_strategy') or '').strip().lower(), 'unclassified')
    if identity == 'non_target' or protocol in {'single_cell', 'other_assay'} or omics not in target_omics:
        selection = 'excluded'
    elif identity == 'target' and protocol in {'not_flagged', 'bulk_confirmed'}:
        selection = 'retained_by_rules'
    else:
        selection = 'review_required'
    selection_reason = ';'.join(x for x in (
        'non_CHO' if identity == 'non_target' else '',
        'identity_unresolved' if identity == 'unresolved' else '',
        'single_cell_protocol' if protocol == 'single_cell' else '',
        'specialized_RNA_assay' if protocol == 'other_assay' else '',
        'protocol_review' if protocol == 'review_required' else '',
        'unsupported_omics' if omics not in target_omics else '',
    ) if x) or 'sample_evidence_without_unresolved_flags'
    out.update(target_status=identity, target_evidence=bool(text_hits), omics=omics,
               identity_rule_ids=';'.join(ids), identity_reason=reason,
               study_context_rule_ids=';'.join(r['rule_id'] for r in context),
               protocol_status=protocol, protocol_flag_fields=';'.join(flags),
               assay_flag_fields=';'.join(assay_flags),
               curation_evidence=decision.get('evidence', '') if decision else '',
               curation_references=';'.join(dict.fromkeys(refs)),
               selection_status=selection, selection_reason=selection_reason,
               selection_policy=POLICY_VERSION)
    return out


def load_decisions(path):
    result = {}
    for r in read_tsv(path):
        acc = r['run_accession']
        if acc in result or not r['reference'] or not r['evidence']:
            raise ValueError('Duplicate or unsupported curation: ' + acc)
        result[acc] = r
    return result


def audit(rows, rules, curated, decisions, target_omics=('RNA-seq','ATAC-seq','ChIP-seq')):
    ids = [r['run_accession'] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate run_accession in input')
    return [classify(r, rules, curated, decisions, target_omics) for r in rows]


def reclassify_dataframe(frame, rules, curated, decisions_path, target_omics):
    import pandas as pd
    records = audit(frame.fillna('').to_dict('records'),
                    rules.fillna('').to_dict('records'), curated.fillna('').to_dict('records'),
                    load_decisions(decisions_path), target_omics)
    return pd.DataFrame(records, index=frame.index,
                        columns=list(dict.fromkeys(list(frame.columns) + list(ADDED))))


def write_tsv(path, rows, fields):
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='raise')
        w.writeheader(); w.writerows(rows)


def save_workflow_audit(frame, outdir, source_paths, target_omics):
    """Persist decisions and source fingerprints whenever metadata is rebuilt.

    Scope is the configured omics, before eligibility exclusion. Upstream FASTQ
    filtering is already represented by the fingerprinted input snapshot.
    """
    records = frame.fillna('').to_dict('records')
    records = [r for r in records if r['omics'] in target_omics]
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for status in ('excluded', 'review_required'):
        write_tsv(outdir/(status+'.tsv'),
                  [r for r in records if r['selection_status']==status], list(frame.columns))
    summary = {
        'policy': POLICY_VERSION,
        'scope': 'Configured target omics after FASTQ availability filtering, before eligibility exclusion',
        'target_omics': list(target_omics),
        'input_rows': len(records),
        'statuses': dict(Counter(r['selection_status'] for r in records)),
        'by_omics': {o: dict(Counter(r['selection_status'] for r in records if r['omics']==o))
                     for o in target_omics},
        'sources': {k: {'path': str(p), 'sha256': hashlib.sha256(Path(p).read_bytes()).hexdigest()}
                    for k,p in source_paths.items()},
    }
    (outdir/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2)+'\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--rules', type=Path, required=True)
    ap.add_argument('--curated', type=Path, required=True)
    ap.add_argument('--decisions', type=Path, required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    args = ap.parse_args()
    # Fresh output directory avoids overwriting any production or previous audit.
    if args.outdir.exists():
        ap.error('Output directory already exists; select a new audit directory')
    rows = read_tsv(args.input)
    result = audit(rows, read_tsv(args.rules), read_tsv(args.curated), load_decisions(args.decisions))
    fields = list(dict.fromkeys(list(rows[0]) + list(ADDED))) if rows else list(ADDED)
    args.outdir.mkdir(parents=True)
    write_tsv(args.outdir/'selection_audit.tsv', result, fields)
    for status in ('retained_by_rules','excluded','review_required'):
        write_tsv(args.outdir/(status+'.tsv'), [r for r in result if r['selection_status']==status], fields)
    summary = {'policy': POLICY_VERSION, 'input_rows': len(rows),
               'scope': 'Only the supplied input; previously excluded runs absent from it are not audited.',
               'statuses': dict(Counter(r['selection_status'] for r in result)),
               'by_omics': {o:dict(Counter(r['selection_status'] for r in result if r['omics']==o))
                            for o in sorted({r['omics'] for r in result})},
               'sha256_inputs': {k:hashlib.sha256(getattr(args,k).read_bytes()).hexdigest()
                                 for k in ('input','rules','curated','decisions')}}
    (args.outdir/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__ == '__main__':
    main()
