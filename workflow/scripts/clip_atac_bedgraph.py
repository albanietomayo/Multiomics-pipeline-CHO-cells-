#!/usr/bin/env python3
"""Clip MACS3 bedGraph signal to reference bounds and emit auditable metrics."""
import argparse, json
from pathlib import Path


def read_sizes(path):
    sizes={}
    with open(path,encoding="utf-8") as handle:
        for n,line in enumerate(handle,1):
            fields=line.rstrip("\n").split("\t")
            if len(fields)<2:raise ValueError(f"Malformed chromosome sizes line {n}")
            if fields[0] in sizes:raise ValueError("Duplicate chromosome: "+fields[0])
            sizes[fields[0]]=int(fields[1])
    if not sizes or any(v<=0 for v in sizes.values()):raise ValueError("Invalid chromosome sizes")
    return sizes


def clip(source,sizes_path,output,metrics_path):
    sizes=read_sizes(sizes_path);metrics={"input_intervals":0,"kept_intervals":0,
      "clipped_intervals":0,"dropped_intervals":0,"unknown_chromosome_intervals":0,"clipped_bp_total":0}
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    with open(source,encoding="utf-8") as src,open(output,"w",encoding="utf-8") as dst:
        previous={}
        for n,line in enumerate(src,1):
            if not line.strip() or line.startswith(("track","browser","#")):continue
            f=line.rstrip("\n").split("\t");metrics["input_intervals"]+=1
            if len(f)<4:raise ValueError(f"Malformed bedGraph line {n}")
            chrom=f[0]
            try:start,end=int(f[1]),int(f[2]);float(f[3])
            except ValueError as exc:raise ValueError(f"Invalid bedGraph line {n}") from exc
            if chrom not in sizes:
                metrics["unknown_chromosome_intervals"]+=1;continue
            clipped_start=max(0,start);clipped_end=min(end,sizes[chrom])
            removed=abs(clipped_start-start)+abs(end-clipped_end)
            if removed:metrics["clipped_intervals"]+=1;metrics["clipped_bp_total"]+=removed
            if clipped_start>=clipped_end:
                metrics["dropped_intervals"]+=1;continue
            if clipped_start<previous.get(chrom,0):raise ValueError("bedGraph is not sorted/nonoverlapping")
            previous[chrom]=clipped_end
            dst.write(f"{chrom}\t{clipped_start}\t{clipped_end}\t{f[3]}\n");metrics["kept_intervals"]+=1
    if metrics["unknown_chromosome_intervals"]:raise ValueError("Unknown chromosome intervals encountered")
    if not metrics["kept_intervals"]:raise ValueError("No bedGraph signal remains after clipping")
    Path(metrics_path).parent.mkdir(parents=True,exist_ok=True)
    Path(metrics_path).write_text(json.dumps(metrics,indent=2,sort_keys=True)+"\n")
    return metrics


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--input",required=True)
    p.add_argument("--chrom-sizes",required=True);p.add_argument("--output",required=True);p.add_argument("--metrics",required=True)
    a=p.parse_args();clip(a.input,a.chrom_sizes,a.output,a.metrics)
if __name__=="__main__":main()
