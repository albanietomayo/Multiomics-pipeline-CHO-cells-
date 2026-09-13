"""Two Snakemake stages: primary attempt, then validated finalization/fallback.

Only a directly reported SIGSEGV from featureCounts permits fallback.
The original wrapper remains responsible for strandedness and paired flags.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def validate(counts, summary):
    with counts.open() as handle:
        reader = csv.reader((line for line in handle if not line.startswith('#')), delimiter='\t')
        header = next(reader)
        if len(header) != 7 or header[:6] != ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']:
            raise ValueError('Unexpected count header')
        seen = set()
        total = 0
        for row in reader:
            if len(row) != 7 or row[0] in seen:
                raise ValueError('Malformed row or duplicate gene')
            seen.add(row[0])
            value = int(row[6])
            if value < 0:
                raise ValueError('Negative count')
            total += value
    if not seen:
        raise ValueError('Empty gene table')
    with summary.open() as handle:
        rows = list(csv.reader(handle, delimiter='\t'))
    assigned = [int(row[1]) for row in rows[1:] if row[0] == 'Assigned']
    if assigned != [total]:
        raise ValueError('Counts do not sum to Assigned')
    return {'genes': len(seen), 'assigned': total}


def execute(folder, wrapper, classification, run, unit, gtf, bam, threads, expected_version):
    folder.mkdir(parents=True, exist_ok=True)
    # Remove only outputs of this attempt, so stale success cannot be reused.
    for name in ['counts.txt', 'counts.txt.summary', 'strandedness.txt', 'signal.json', 'status.json']:
        (folder / name).unlink(missing_ok=True)
    version_result = subprocess.run(['featureCounts', '-v'], text=True, capture_output=True, check=True)
    version = (version_result.stdout + version_result.stderr).strip()
    if expected_version not in version.split():
        # Subread normally reports 'featureCounts v2.1.1'.
        if ('v' + expected_version) not in version.split():
            raise RuntimeError(f'Unexpected featureCounts version: {version}')
    env = os.environ.copy()
    env['FEATURECOUNTS_SIGNAL_REPORT'] = str((folder / 'signal.json').resolve())
    command = [sys.executable, wrapper, classification, run, unit, gtf, bam,
               str(threads), str(folder / 'counts.txt'), str(folder / 'strandedness.txt')]
    with (folder / 'execution.log').open('w') as log:
        result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
    signal_report = folder / 'signal.json'
    code = json.loads(signal_report.read_text())['returncode'] if signal_report.exists() else result.returncode
    record = {'wrapper_exit': result.returncode, 'featurecounts_exit': code,
              'version': version, 'threads': int(threads)}
    if result.returncode == 0:
        record.update(validate(folder / 'counts.txt', folder / 'counts.txt.summary'))
        if not (folder / 'strandedness.txt').is_file():
            raise ValueError('Missing strandedness record')
        record['status'] = 'success'
    elif signal_report.exists() and code == -11:
        record['status'] = 'sigsegv'
    else:
        print((folder / 'execution.log').read_text(), file=sys.stderr)
        raise RuntimeError(f'Non-recoverable failure {record}; see {folder / "execution.log"}')
    (folder / 'status.json').write_text(json.dumps(record, indent=2) + '\n')
    return record


def main():
    mode, primary, wrapper, classification, run, unit, gtf, bam, threads, counts, strandedness, provenance = sys.argv[1:]
    primary = Path(primary)
    if mode == 'primary':
        execute(primary, wrapper, classification, run, unit, gtf, bam, threads, '2.1.1')
        return
    if mode != 'finalize':
        raise ValueError(mode)
    initial = json.loads((primary / 'status.json').read_text())
    chosen = primary
    attempts = [initial]
    if initial['status'] == 'sigsegv':
        chosen = Path(counts).parent / '.fallback'
        # One thread: the configuration validated on DRR091439.
        fallback = execute(chosen, wrapper, classification, run, unit, gtf, bam, 1, '2.0.6')
        attempts.append(fallback)
        if fallback['status'] != 'success':
            raise RuntimeError('Fallback also failed; inspect .fallback/execution.log')
    elif initial['status'] != 'success':
        raise ValueError('Invalid primary status')
    metrics = validate(chosen / 'counts.txt', chosen / 'counts.txt.summary')
    Path(counts).parent.mkdir(parents=True, exist_ok=True)
    for src, dst in [('counts.txt', counts), ('counts.txt.summary', counts + '.summary'), ('strandedness.txt', strandedness)]:
        shutil.copy2(chosen / src, dst)
    record = {'run_accession': run, 'alignment_unit': unit,
              'fallback_used': len(attempts) == 2, 'attempts': attempts,
              'selected_version': attempts[-1]['version'], **metrics,
              'counts_sha256': hashlib.sha256(Path(counts).read_bytes()).hexdigest()}
    Path(provenance).write_text(json.dumps(record, indent=2) + '\n')


if __name__ == '__main__':
    main()
