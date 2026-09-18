#!/usr/bin/env python3
"""Create canonical 1-bp insertions and pilot-compatible shifted read intervals."""
import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path

from file_hash import sha256_file

def insertion_position(start, end, reverse, forward_shift=4, reverse_shift=-5):
    return end + reverse_shift if reverse else start + forward_shift


def external_bed_sort(source, destination, tmpdir):
    """Use GNU sort's bounded-memory/external merge implementation."""
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    command = [
        "sort", "--stable", "--temporary-directory", str(tmpdir),
        "-k1,1", "-k2,2n", "-k3,3n", "-k4,4", "-k6,6", "-k5,5n",
        str(source),
    ]
    with open(destination, "w", encoding="utf-8") as output:
        subprocess.run(command, stdout=output, check=True, env=environment)


def transform(bam_path, insertions_path, tagalign_path, metrics_path,
              forward_shift=4, reverse_shift=-5, tmpdir=None):
    import pysam
    metrics = {"schema_version": 1, "representation": {
        "insertions": "canonical_1bp_tn5_events",
        "tagalign": "read_intervals_with_only_the_5prime_end_tn5_shifted",
    }, "tn5_forward_shift": forward_shift, "tn5_reverse_shift": reverse_shift,
        "input_records": 0, "insertion_events": 0, "shifted_intervals": 0,
        "dropped_invalid_coordinates": 0}
    Path(insertions_path).parent.mkdir(parents=True, exist_ok=True)
    Path(tagalign_path).parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="atac_tn5_", dir=tmpdir) as folder:
        folder = Path(folder)
        raw_insertions = folder / "insertions.unsorted.bed"
        raw_tagalign = folder / "tagalign.unsorted.bed"
        with pysam.AlignmentFile(bam_path, "rb") as bam, \
                raw_insertions.open("w", encoding="utf-8") as ins, \
                raw_tagalign.open("w", encoding="utf-8") as tags:
            lengths = dict(zip(bam.references, bam.lengths))
            if not lengths: raise ValueError("BAM has no reference sequences")
            for read in bam.fetch(until_eof=True):
                metrics["input_records"] += 1
                if read.is_unmapped or read.reference_name not in lengths or read.reference_end is None:
                    raise ValueError("Filtered BAM contains an unmapped or unknown-reference record")
                chrom = read.reference_name; length = lengths[chrom]
                pos = insertion_position(read.reference_start, read.reference_end, read.is_reverse,
                                         forward_shift, reverse_shift)
                start = read.reference_start if read.is_reverse else pos
                end = pos if read.is_reverse else read.reference_end
                if pos < 0 or pos >= length or start < 0 or end > length or start >= end:
                    metrics["dropped_invalid_coordinates"] += 1
                    continue
                strand = "-" if read.is_reverse else "+"
                ins.write(f"{chrom}\t{pos}\t{pos + 1}\t{read.query_name}\t{read.mapping_quality}\t{strand}\n")
                tags.write(f"{chrom}\t{start}\t{end}\t{read.query_name}\t{read.mapping_quality}\t{strand}\n")
                metrics["insertion_events"] += 1; metrics["shifted_intervals"] += 1
        external_bed_sort(raw_insertions, insertions_path, folder)
        external_bed_sort(raw_tagalign, tagalign_path, folder)
    if metrics["input_records"] == 0 or metrics["shifted_intervals"] == 0:
        raise ValueError("Tn5 transformation produced no intervals")
    for key, path in (("input_bam", bam_path), ("insertions", insertions_path), ("tagalign", tagalign_path)):
        metrics[key] = {"path": str(path), "sha256": sha256_file(path)}
    Path(metrics_path).parent.mkdir(parents=True, exist_ok=True)
    Path(metrics_path).write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    return metrics


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam",required=True);p.add_argument("--insertions",required=True)
    p.add_argument("--tagalign",required=True);p.add_argument("--metrics",required=True)
    p.add_argument("--forward-shift",type=int,default=4);p.add_argument("--reverse-shift",type=int,default=-5)
    p.add_argument("--tmpdir")
    a=p.parse_args();transform(a.bam,a.insertions,a.tagalign,a.metrics,a.forward_shift,a.reverse_shift,a.tmpdir)
if __name__=="__main__":main()
