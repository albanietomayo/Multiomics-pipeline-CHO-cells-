#!/usr/bin/env python3
"""Check filtered ATAC BAM and consistency of TSS outputs; no quality cutoff.

Checks do not recalculate a TSS profile. Provenance and input hashes identify
which outputs were checked. A successful report is not biological acceptance.
"""
import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import tempfile


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def table(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def check_tss(profile, summary, records, window, bin_size, background):
    rows = table(summary)
    if not rows or any(set(r) != {'metric', 'value'} for r in rows):
        raise ValueError('Invalid summary columns')
    if len({r['metric'] for r in rows}) != len(rows):
        raise ValueError('Duplicate summary metric')
    values = {r['metric']: float(r['value']) for r in rows}
    if any(not math.isfinite(v) or v < 0 for v in values.values()):
        raise ValueError('Nonfinite or negative TSS summary')
    integer_fields = ('total_tss', 'usable_tss', 'tss_skipped_boundary',
                      'tss_skipped_missing_reference', 'reads_seen', 'reads_used',
                      'insertion_events_in_tss_windows', 'window_bp', 'bin_size_bp',
                      'background_bp_per_side')
    if any(values[k] != int(values[k]) for k in integer_fields):
        raise ValueError('Noninteger TSS count')
    if values['reads_seen'] != records or not 0 <= values['reads_used'] <= records:
        raise ValueError('TSS read counts disagree with BAM')
    if values['usable_tss'] <= 0 or (values['usable_tss'] + values['tss_skipped_boundary']
                                  + values['tss_skipped_missing_reference'] != values['total_tss']):
        raise ValueError('Inconsistent TSS reference counts')
    if (values['window_bp'], values['bin_size_bp'], values['background_bp_per_side']) != (window, bin_size, background):
        raise ValueError('TSS parameters differ from requested workflow configuration')
    if min(window, bin_size, background) <= 0 or 2 * window % bin_size or background % bin_size or background > window:
        raise ValueError('Invalid profile geometry')
    bins = table(profile)
    if len(bins) != 2 * window // bin_size:
        raise ValueError('Unexpected TSS bin count')
    counts, signal = [], []
    for i, row in enumerate(bins):
        start, end, count = (int(row[k]) for k in ('relative_start', 'relative_end', 'insertion_count'))
        v = float(row['normalized_signal'])
        if start != -window + i * bin_size or end != start + bin_size or count < 0 or not math.isfinite(v) or v < 0:
            raise ValueError('Invalid TSS bin')
        counts.append(count); signal.append(v)
    n = background // bin_size
    mean = sum(counts[:n] + counts[-n:]) / (2 * n)
    if mean <= 0 or not math.isclose(mean, values['background_mean'], abs_tol=1e-6):
        raise ValueError('TSS background mismatch')
    if sum(counts) != values['insertion_events_in_tss_windows']:
        raise ValueError('TSS insertion total mismatch')
    if any(not math.isclose(v, count / mean, abs_tol=1.1e-6, rel_tol=1e-8) for count, v in zip(counts, signal)):
        raise ValueError('TSS normalization mismatch')
    if not math.isclose(max(signal), values['tss_enrichment_score'], abs_tol=1.1e-6, rel_tol=1e-8):
        raise ValueError('TSS maximum differs from score')
    return values


def check_bam(path, index, layout, min_mapq, mitochondrial, tmpdir):
    import pysam
    pysam.quickcheck(str(path))
    total = 0
    with pysam.AlignmentFile(str(path), 'rb', index_filename=str(index), require_index=True) as bam:
        bam.check_index()
        if bam.header.to_dict().get('HD', {}).get('SO') != 'coordinate':
            raise ValueError('BAM is not declared coordinate sorted')
        previous = (-1, -1)
        for read in bam.fetch(until_eof=True):
            total += 1
            excluded = 3852 if layout == 'PAIRED' else 3844
            if read.flag & excluded or read.mapping_quality < min_mapq or read.reference_name == mitochondrial:
                raise ValueError('BAM contains an alignment that fails filtering')
            if layout == 'SINGLE' and read.is_paired:
                raise ValueError('Paired record in SINGLE library')
            if layout == 'PAIRED' and read.flag & 3 != 3:
                raise ValueError('Record is not a proper pair')
            position = (read.reference_id, read.reference_start)
            if position < previous:
                raise ValueError('BAM coordinates are not sorted')
            previous = position
        if sum(s.mapped for s in bam.get_index_statistics()) != total:
            raise ValueError('BAM index record count mismatch')
    if total == 0:
        raise ValueError('Empty filtered BAM')
    if layout == 'PAIRED':
        if total % 2:
            raise ValueError('Odd paired record count')
        with tempfile.TemporaryDirectory(prefix='atac_pair_check_', dir=tmpdir) as folder:
            named = str(Path(folder) / 'names.bam')
            pysam.sort('-n', '-m', '256M', '-T', str(Path(folder) / 'sort'), '-o', named, str(path))
            with pysam.AlignmentFile(named, 'rb') as bam:
                for name, group in itertools.groupby(bam.fetch(until_eof=True), lambda r: r.query_name):
                    pair = list(group)
                    if len(pair) != 2 or sum(r.is_read1 for r in pair) != 1 or sum(r.is_read2 for r in pair) != 1:
                        raise ValueError('Incomplete or malformed read pair: ' + str(name))
                    a, b = pair
                    if (a.next_reference_id, a.next_reference_start) != (b.reference_id, b.reference_start) or (b.next_reference_id, b.next_reference_start) != (a.reference_id, a.reference_start):
                        raise ValueError('Inconsistent mate coordinates')
    return total


def atomic_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=path.parent, prefix='.' + path.name + '.')
    try:
        with os.fdopen(fd, 'w') as handle:
            json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def main():
    from selection_gate import require, fingerprint
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('run', 'bam', 'bai', 'profile', 'summary', 'output'):
        p.add_argument('--' + key, required=True)
    p.add_argument('--layout', choices=['SINGLE', 'PAIRED'], required=True)
    p.add_argument('--min-mapq', type=int, default=30)
    p.add_argument('--mitochondrial', default='NC_007936.1')
    p.add_argument('--window', type=int, default=2000)
    p.add_argument('--bin-size', type=int, default=10)
    p.add_argument('--background', type=int, default=100)
    p.add_argument('--tmpdir', default=None)
    a = p.parse_args()
    require(a.run, omics='ATAC-seq')
    if a.min_mapq < 0: raise ValueError('Negative MAPQ')
    count = check_bam(a.bam, a.bai, a.layout, a.min_mapq, a.mitochondrial, a.tmpdir)
    metrics = check_tss(a.profile, a.summary, count, a.window, a.bin_size, a.background)
    inputs = {k: {'path': getattr(a, k), 'sha256': sha256(getattr(a, k))}
              for k in ('bam', 'bai', 'profile', 'summary')}
    atomic_json(a.output, {'run_accession': a.run, 'layout': a.layout,
        'technical_checks': 'passed', 'biological_acceptance': 'not_assessed',
        'total_records': count, 'min_mapq': a.min_mapq, 'mitochondrial_records': 0,
        'tss_metrics': metrics, 'inputs': inputs, 'selection_inputs': fingerprint()})


if __name__ == '__main__': main()
