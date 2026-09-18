#!/usr/bin/env python3
"""Layout-aware ATAC BAM filtering with rigorous complete-pair validation."""
import argparse
import itertools
import json
import tempfile
from pathlib import Path
from file_hash import sha256_file

def rejected(read, min_mapq, exclude_flags, mitochondrial):
    return (read.flag & exclude_flags or read.mapping_quality < min_mapq or
            read.reference_name == mitochondrial)


def validate_pair(pair):
    if len(pair) != 2:
        return False
    if len({r.query_name for r in pair}) != 1:
        return False
    if sum(r.is_read1 for r in pair) != 1 or sum(r.is_read2 for r in pair) != 1:
        return False
    a, b = pair
    return ((a.next_reference_id, a.next_reference_start) == (b.reference_id, b.reference_start)
            and (b.next_reference_id, b.next_reference_start) == (a.reference_id, a.reference_start)
            and a.reference_id == b.reference_id)


def nuclear_references(references, mitochondrial):
    result = [name for name in references if name not in {mitochondrial, "*"}]
    if not result:
        raise ValueError("No valid nuclear references exist")
    return result


def filter_bam(source, output, layout, min_mapq, single_exclude, paired_require,
               paired_exclude, mitochondrial, tmpdir, metrics_path=None):
    import pysam
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(source, "rb") as bam:
        nuclear_references(bam.references, mitochondrial)

    with pysam.AlignmentFile(source, "rb") as bam:
        input_records = sum(1 for _ in bam.fetch(until_eof=True))
    kept = 0
    with tempfile.TemporaryDirectory(prefix="atac_filter_", dir=tmpdir) as folder:
        folder = Path(folder)
        unsorted = folder / "kept.bam"
        if layout == "SINGLE":
            with pysam.AlignmentFile(source, "rb") as bam, pysam.AlignmentFile(unsorted, "wb", template=bam) as out:
                for read in bam.fetch(until_eof=True):
                    if not rejected(read, min_mapq, single_exclude, mitochondrial):
                        if read.is_paired:
                            raise ValueError("Paired record encountered in SINGLE library")
                        out.write(read); kept += 1
        elif layout == "PAIRED":
            named = folder / "queryname.bam"
            pysam.sort("-n", "-T", str(folder / "namesort"), "-o", str(named), str(source))
            with pysam.AlignmentFile(named, "rb") as bam, pysam.AlignmentFile(unsorted, "wb", template=bam) as out:
                candidates = (r for r in bam.fetch(until_eof=True)
                              if not rejected(r, min_mapq, paired_exclude, mitochondrial)
                              and r.flag & paired_require == paired_require)
                for _, records in itertools.groupby(candidates, lambda r: r.query_name):
                    pair = list(records)
                    if validate_pair(pair):
                        for read in pair: out.write(read)
                        kept += 2
        else:
            raise ValueError("Unsupported layout: " + layout)
        if kept == 0:
            raise ValueError("ATAC filtering produced no valid nuclear records")
        pysam.sort("-T", str(folder / "coordsort"), "-o", str(output), str(unsorted))
    pysam.index(str(output))
    pysam.quickcheck(str(output))
    if metrics_path:
        metrics={"layout":layout,"input_records":input_records,"output_records":kept,
          "min_mapq":min_mapq,"single_exclude_flags":single_exclude,
          "paired_require_flags":paired_require,"paired_exclude_flags":paired_exclude,
          "mitochondrial_accession":mitochondrial,
          "output_bam_sha256":sha256_file(output)}
        path=Path(metrics_path);path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(metrics,indent=2,sort_keys=True)+"\n")
    return kept


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True); p.add_argument("--output", required=True)
    p.add_argument("--layout", choices=["SINGLE", "PAIRED"], required=True)
    p.add_argument("--min-mapq", type=int, default=30)
    p.add_argument("--single-exclude", type=int, default=3844)
    p.add_argument("--paired-require", type=int, default=3)
    p.add_argument("--paired-exclude", type=int, default=3852)
    p.add_argument("--mitochondrial", default="NC_007936.1")
    p.add_argument("--tmpdir", required=True)
    p.add_argument("--metrics", required=True)
    a = p.parse_args()
    filter_bam(a.input, a.output, a.layout, a.min_mapq, a.single_exclude,
               a.paired_require, a.paired_exclude, a.mitochondrial, a.tmpdir, a.metrics)


if __name__ == "__main__": main()
