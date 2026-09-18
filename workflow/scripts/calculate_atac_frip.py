#!/usr/bin/env python3
"""Deterministically count shifted read intervals overlapping at least one peak."""
import argparse, json, platform
from collections import defaultdict
from pathlib import Path
from file_hash import sha256_file


def load_intervals(path, min_columns=3):
    data=defaultdict(list)
    with open(path,encoding="utf-8") as handle:
        for number,line in enumerate(handle,1):
            if not line.strip() or line.startswith("#"):continue
            fields=line.rstrip("\n").split("\t")
            if len(fields)<min_columns:raise ValueError(f"Malformed interval line {number}")
            try:start,end=int(fields[1]),int(fields[2])
            except ValueError as exc:raise ValueError(f"Noninteger coordinates at line {number}") from exc
            if start<0 or start>=end:raise ValueError(f"Invalid coordinates at line {number}")
            data[fields[0]].append((start,end))
    for chrom in data:data[chrom].sort()
    return data


def count_frip(tags,peaks):
    peak_data=load_intervals(peaks);total=overlap=0
    indexes=defaultdict(int);previous_key=None
    with open(tags,encoding="utf-8") as handle:
        for number,line in enumerate(handle,1):
            if not line.strip() or line.startswith("#"):continue
            fields=line.rstrip("\n").split("\t")
            if len(fields)<3:raise ValueError(f"Malformed interval line {number}")
            try:start,end=int(fields[1]),int(fields[2])
            except ValueError as exc:raise ValueError(f"Noninteger coordinates at line {number}") from exc
            chrom=fields[0]
            if start<0 or start>=end:raise ValueError(f"Invalid coordinates at line {number}")
            key=(chrom,start,end)
            if previous_key is not None and key<previous_key:raise ValueError("TagAlign is not coordinate sorted")
            previous_key=key;total+=1
            chrom_peaks=peak_data.get(chrom,[]);j=indexes[chrom]
            while j<len(chrom_peaks) and chrom_peaks[j][1]<=start:j+=1
            indexes[chrom]=j
            if j<len(chrom_peaks) and chrom_peaks[j][0]<end:overlap+=1
    if total==0:raise ValueError("FRiP denominator is zero")
    return total,overlap,overlap/total


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--tagalign",required=True)
    p.add_argument("--peaks",required=True);p.add_argument("--output",required=True);p.add_argument("--run",required=True)
    a=p.parse_args();total,overlap,value=count_frip(a.tagalign,a.peaks)
    result={"run_accession":a.run,"usable_records_denominator":total,
      "records_overlapping_at_least_one_peak":overlap,"frip":value,
      "representation":"tn5_shifted_read_intervals_pilot_semantics",
      "software":{"name":"pipeline-python","python":platform.python_version()},
      "inputs":{k:{"path":str(v),"sha256":sha256_file(v)}
                for k,v in (("tagalign",a.tagalign),("peak_set",a.peaks))}}
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()
