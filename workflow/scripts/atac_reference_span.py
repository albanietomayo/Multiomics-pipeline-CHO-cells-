#!/usr/bin/env python3
"""Create chromosome sizes and nuclear-reference-span approximation from FAI."""
import argparse, json
from pathlib import Path


def calculate(fai,mitochondrial):
    rows=[]
    with open(fai,encoding="utf-8") as handle:
        for n,line in enumerate(handle,1):
            f=line.rstrip("\n").split("\t")
            if len(f)<2:raise ValueError(f"Malformed FAI line {n}")
            length=int(f[1])
            if f[0] not in {mitochondrial,"*"} and length>0:rows.append((f[0],length))
    if not rows:raise ValueError("No valid nuclear references exist")
    return rows,sum(x[1] for x in rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--fai",required=True)
    p.add_argument("--mitochondrial",required=True);p.add_argument("--sizes",required=True);p.add_argument("--metadata",required=True)
    a=p.parse_args();rows,span=calculate(a.fai,a.mitochondrial)
    Path(a.sizes).parent.mkdir(parents=True,exist_ok=True)
    Path(a.sizes).write_text("".join(f"{c}\t{n}\n" for c,n in rows))
    Path(a.metadata).write_text(json.dumps({"nuclear_reference_span":span,
      "definition":"sum of nuclear reference lengths; approximation, not mappability-derived effective genome size",
      "mitochondrial_accession_excluded":a.mitochondrial,"nuclear_contigs":len(rows)},indent=2)+"\n")
if __name__=="__main__":main()
