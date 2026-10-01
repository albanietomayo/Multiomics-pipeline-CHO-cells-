#!/usr/bin/env python3
"""Validate a curated ChIP eligibility decision; never infer new approvals."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

SNAPSHOT = Path('snapshots/chipseq/control_validation_001')
PROVENANCE_FILES = (
    str(SNAPSHOT / 'chipseq_conditions.tsv'),
    str(SNAPSHOT / 'control_candidate_summary.json'),
    'config/samples.tsv',
    str(SNAPSHOT / 'request_plan_provenance.json'),
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(path):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    return indexed(rows, str(path))


def indexed(rows, label):
    result = {}
    for row in rows:
        run = row['run_accession']
        if not run or run in result:
            raise ValueError(f'Missing or duplicate run in {label}: {run}')
        result[run] = row
    return result


def validate(args):
    review = json.loads(args.decisions.read_text(encoding='utf-8'))
    if review['schema_version'] != 1 or review['scope'] != 'experimental_protocol_eligibility_only':
        raise ValueError('Unsupported eligibility review schema or scope')
    provenance = review['input_sha256']
    if set(provenance) != set(PROVENANCE_FILES):
        raise ValueError('Unexpected provenance file set')
    for name in PROVENANCE_FILES:
        if sha(name) != provenance[name]:
            raise ValueError(f'Review evidence changed; a new review is required: {name}')
    if sha(args.samples) != provenance['config/samples.tsv']:
        raise ValueError('Active catalog differs from reviewed catalog')
    frozen_conditions = str(SNAPSHOT / 'chipseq_conditions.tsv')
    if sha(args.conditions) != provenance[frozen_conditions]:
        raise ValueError('Active conditions differ from reviewed conditions')

    samples = table(args.samples)
    cohort = {run: row for run, row in samples.items() if row['omics'] == 'ChIP-seq'}
    decisions = indexed(review['runs'], 'eligibility review')
    conditions = table(args.conditions)
    if not cohort or set(decisions) != set(cohort) or set(conditions) != set(cohort):
        raise ValueError('Catalog, conditions and review must cover exactly the same ChIP runs')
    states = {'retained_by_rules', 'excluded', 'review_required'}
    for run, decision in decisions.items():
        if decision['experimental_eligibility'] not in states:
            raise ValueError(f'Unknown eligibility state: {run}')
        if decision['study_accession'] != cohort[run]['study_accession']:
            raise ValueError(f'Review/catalog study mismatch: {run}')
        for key in ('study_accession', 'library_role', 'declared_target'):
            if decision[key] != conditions[run][key]:
                raise ValueError(f'Review/conditions mismatch in {key}: {run}')
        if not decision.get('reason') or not decision.get('evidence_scope'):
            raise ValueError(f'Missing decision evidence: {run}')
        if not review['sources'].get(decision['study_accession'], {}).get('url'):
            raise ValueError(f'Missing study source: {run}')

    selected = []
    if args.pilot:
        pilot = json.loads(args.pilot.read_text(encoding='utf-8'))
        if pilot['schema_version'] != 1:
            raise ValueError('Unsupported pilot schema')
        if pilot['provenance_sha256']['samples'] != sha(args.samples):
            raise ValueError('Pilot catalog provenance mismatch')
        expected = {pilot['ip_run_accession']: 'ip', pilot['input_run_accession']: 'input'}
        declared = indexed(pilot['runs'], 'pilot')
        if len(expected) != 2 or {run: row['role'] for run, row in declared.items()} != expected:
            raise ValueError('Pilot must specify one IP and one distinct input')
        for run, role in expected.items():
            decision = decisions.get(run)
            if decision is None or decision['experimental_eligibility'] != 'retained_by_rules':
                raise ValueError(f'Pilot run is not approved for this protocol: {run}')
            if decision['library_role'] != role or decision['study_accession'] != pilot['study_accession']:
                raise ValueError(f'Pilot role or study mismatch: {run}')
            if role == 'ip' and decision['declared_target'] != pilot['target']:
                raise ValueError(f'Pilot histone target mismatch: {run}')
            selected.append(run)

    report = {
        'schema_version': 1,
        'scope': 'protocol_eligibility_only',
        'reviewed_runs': len(decisions),
        'eligibility_counts': dict(sorted(Counter(
            row['experimental_eligibility'] for row in decisions.values()).items())),
        'selected_runs': sorted(selected),
        'control_pairing': 'not_approved_by_this_check',
        'technical_qc': 'not_assessed_by_this_check',
        'decisions_sha256': sha(args.decisions),
        'samples_sha256': sha(args.samples),
        'conditions_sha256': sha(args.conditions),
        'validator_sha256': sha(__file__),
    }
    if args.pilot:
        report['pilot_sha256'] = sha(args.pilot)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--decisions', type=Path, required=True)
    parser.add_argument('--samples', type=Path, required=True)
    parser.add_argument('--conditions', type=Path, required=True)
    parser.add_argument('--pilot', type=Path)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    report = validate(args)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.report.with_name(args.report.name + '.part')
        temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        temporary.replace(args.report)
    print(f"[OK] Eligibility provenance and {report['reviewed_runs']} reviewed runs verified.")
    if args.pilot:
        print('[OK] Pilot protocol eligibility verified: ' + ', '.join(report['selected_runs']))


if __name__ == '__main__':
    main()
