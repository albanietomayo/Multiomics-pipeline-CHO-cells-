#!/usr/bin/env python3
"""Preserve deterministic IDR peak sets at both requested thresholds."""
import argparse,json
from pathlib import Path
from file_hash import sha256_file


def filter_idr(source,out10,out05,metrics):
    n10=n05=total=0
    Path(out10).parent.mkdir(parents=True,exist_ok=True)
    with open(source,encoding="utf-8") as src,open(out10,"w",encoding="utf-8") as a,open(out05,"w",encoding="utf-8") as b:
        for n,line in enumerate(src,1):
            if not line.strip() or line.startswith("#"):continue
            f=line.rstrip("\n").split("\t")
            if len(f)<12:raise ValueError(f"IDR output line {n} has fewer than 12 columns")
            try:score=float(f[11])
            except ValueError as exc:raise ValueError(f"Invalid IDR score line {n}") from exc
            total+=1
            # IDR's column 12 is -log10(global IDR).
            if score>=1.0:a.write(line);n10+=1
            if score>=1.301029995664:b.write(line);n05+=1
    if total==0:raise ValueError("Empty IDR output")
    Path(metrics).write_text(json.dumps({"all_ranked_peaks":total,"idr_le_0.10":n10,
      "idr_le_0.05":n05,"thresholds_preserved":[0.10,0.05],"software":{"idr":"2.0.4.2","seed":0},
      "sha256":{"idr_ranked_input":sha256_file(source),"idr_le_0.10":sha256_file(out10),"idr_le_0.05":sha256_file(out05)}},indent=2)+"\n")


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--input",required=True)
    p.add_argument("--idr-010",required=True);p.add_argument("--idr-005",required=True);p.add_argument("--metrics",required=True)
    a=p.parse_args();filter_idr(a.input,a.idr_010,a.idr_005,a.metrics)
if __name__=="__main__":main()
