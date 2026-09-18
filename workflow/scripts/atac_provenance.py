#!/usr/bin/env python3
"""Collect portable hashes and configured ATAC semantics into one JSON record."""
import argparse,json
from pathlib import Path
from file_hash import sha256_file


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--run",required=True)
    p.add_argument("--manifest",required=True);p.add_argument("--reference-metadata",required=True)
    p.add_argument("--input",action="append",default=[]);p.add_argument("--output",required=True)
    p.add_argument("--parameters-json",required=True);a=p.parse_args()
    from atac_manifest import eligible_run
    rows=eligible_run(a.manifest,a.run)
    value={"schema_version":1,"run_accession":a.run,"selection_decision":"eligible",
      "source_fastqs":[{"filename":r["filename"],"md5":r["fastq_md5"],"bytes":r["fastq_bytes"]} for r in rows],
      "parameters":json.loads(a.parameters_json),
      "inputs":{Path(x).name:{"path":str(Path(x)),"sha256":sha256_file(x)} for x in a.input},
      "reference_metadata":{"path":str(Path(a.reference_metadata)),"sha256":sha256_file(a.reference_metadata)}}
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(value,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()
