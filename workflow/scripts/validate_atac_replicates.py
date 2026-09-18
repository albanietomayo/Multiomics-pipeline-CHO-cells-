#!/usr/bin/env python3
"""Validate an explicit ATAC replicate group and record its eligibility."""
import argparse,json
from pathlib import Path
from atac_manifest import load_groups
from file_hash import sha256_file


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--groups",required=True)
    p.add_argument("--manifest",required=True);p.add_argument("--group",required=True);p.add_argument("--output",required=True)
    a=p.parse_args();groups=load_groups(a.groups,a.manifest)
    if a.group not in groups:raise ValueError("Unknown replicate group: "+a.group)
    value={"replicate_group":a.group,"eligible":True,"members":groups[a.group],
      "replicate_type":"unspecified_comparable_replicates",
      "group_table_sha256":sha256_file(a.groups),
      "production_manifest_sha256":sha256_file(a.manifest)}
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(value,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()
