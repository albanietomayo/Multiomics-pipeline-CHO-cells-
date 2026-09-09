#!/usr/bin/env python3
"""Save an audited complete-label subset without inferring missing cell lines."""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import sys
import pandas as pd
from collect_sample_counts import DEFAULT_OUTPUT, missing_mask

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    raw = pd.read_csv(args.metadata, sep="\t", dtype=str, keep_default_na=False)
    missing = missing_mask(raw[["sample_accession", "cell_line", "omics"]])
    excluded = missing.any(axis=1)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw.loc[~excluded].to_csv(args.output_dir / "metadata_complete.tsv", sep="\t", index=False)
    review = raw.loc[excluded].copy()
    review["exclusion_reason"] = missing.loc[excluded].apply(lambda row: "Missing: " + ", ".join(row.index[row]), axis=1)
    review.to_csv(args.output_dir / "excluded_metadata.tsv", sep="\t", index=False)
    summary = {"source": str(args.metadata.resolve()), "source_sha256": hashlib.sha256(args.metadata.read_bytes()).hexdigest(),
               "command": shlex.join([sys.executable, *sys.argv]), "source_rows": len(raw),
               "included_rows": int((~excluded).sum()), "excluded_rows": int(excluded.sum()),
               "missing_by_field": missing.sum().to_dict(), "excluded_unique_sample_ids": raw.loc[excluded, "sample_accession"].nunique(),
               "excluded_assay_rows": raw.loc[excluded, "omics"].value_counts().to_dict(),
               "policy": "Complete required fields only. No cell-line inference or label merging. Plot limited to annotated samples."}
    (args.output_dir / "metadata_coverage.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
