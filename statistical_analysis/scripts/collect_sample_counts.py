#!/usr/bin/env python3
"""Count biological samples using authoritative CSV/TSV metadata."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import shlex
import sys
import pandas as pd

DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "results"

def versions():
    result = {"python": sys.version.split()[0]}
    for name in ("pandas", "plotly", "kaleido"):
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = "not installed"
    return result

def normalize(label):
    key = re.sub(r"[^a-z0-9]", "", label.lower())
    if key in {"rna", "rnaseq", "scrna", "scrnaseq"}:
        return ["RNA"], False
    if key in {"atac", "atacseq", "scatac", "scatacseq"}:
        return ["ATAC"], False
    if key in {"multiome", "rnaatac", "atacrna", "rnaseqatacseq", "atacseqrnaseq", "scrnascatac", "scatacscrna"}:
        return ["RNA", "ATAC"], False
    return [label], True

def missing_mask(frame):
    return frame.apply(lambda s: s.fillna("").astype(str).str.strip().str.lower().isin({"", "na", "n/a", "nan", "null", "none"}))

def collect(frame):
    frame = frame.copy()
    for col in frame:
        frame[col] = frame[col].fillna("").astype(str).str.strip()
    missing = missing_mask(frame)
    if missing.any().any():
        raise ValueError(f"Missing required values: {missing.sum().to_dict()}. Correct metadata or supply an explicitly audited complete subset.")
    if frame.empty:
        raise ValueError("Metadata has no sample rows.")
    conflicts = frame.groupby("sample_id").cell_line.nunique()
    if (conflicts > 1).any():
        raise ValueError("Conflicting cell lines for sample IDs: " + ", ".join(conflicts[conflicts > 1].index[:20]))
    unknown, expanded = set(), []
    for row in frame.itertuples(index=False):
        assays, is_unknown = normalize(row.sequencing_type)
        if is_unknown:
            unknown.add(row.sequencing_type)
        expanded.extend((row.sample_id, row.cell_line, assay) for assay in assays)
    inventory = pd.DataFrame(expanded, columns=frame.columns).drop_duplicates().sort_values(list(frame.columns)).reset_index(drop=True)
    grid = pd.MultiIndex.from_product([sorted(inventory.cell_line.unique()), sorted(inventory.sequencing_type.unique())], names=["cell_line", "sequencing_type"])
    counts = inventory.groupby(["cell_line", "sequencing_type"]).sample_id.nunique().reindex(grid, fill_value=0).rename("count").reset_index()
    return inventory, counts, len(expanded) - len(inventory), sorted(unknown)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--sample-id-col", default="sample_id")
    parser.add_argument("--cell-line-col", default="cell_line")
    parser.add_argument("--sequencing-type-col", default="sequencing_type")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    try:
        raw = pd.read_csv(args.metadata, sep="\t" if args.metadata.suffix.lower() in {".tsv", ".tab"} else ",", dtype=str, keep_default_na=False)
        selected = [args.sample_id_col, args.cell_line_col, args.sequencing_type_col]
        if len(set(selected)) != 3:
            raise ValueError("Supply three different column names.")
        absent = set(selected) - set(raw.columns)
        if absent:
            raise ValueError(f"Required columns absent: {sorted(absent)}. Available: {list(raw.columns)}")
        frame = raw[selected].copy()
        frame.columns = ["sample_id", "cell_line", "sequencing_type"]
        inventory, counts, removed, unknown = collect(frame)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        inventory.to_csv(args.output_dir / "sample_assay_inventory.csv", index=False)
        counts.to_csv(args.output_dir / "counts_by_cell_line_and_assay.csv", index=False)
        counts.pivot(index="cell_line", columns="sequencing_type", values="count").to_csv(args.output_dir / "counts_pivot.csv")
        summary = {"metadata": str(args.metadata.resolve()), "metadata_sha256": hashlib.sha256(args.metadata.read_bytes()).hexdigest(),
                   "command": shlex.join([sys.executable, *sys.argv]), "versions": versions(), "input_rows": len(raw),
                   "missing_required_values": 0, "duplicate_sample_assay_rows_removed": removed, "unknown_assay_labels_preserved": unknown,
                   "unique_biological_samples": inventory.sample_id.nunique(), "sample_assay_assignments": len(inventory),
                   "cell_lines": inventory.cell_line.nunique(), "assay_totals": counts.groupby("sequencing_type")["count"].sum().to_dict(),
                   "note": "Metadata availability; downloaded/processed availability not assessed. Multiome counts in both assays. Cell-line labels preserved after trimming."}
        (args.output_dir / "validation_summary.txt").write_text(json.dumps(summary, indent=2) + "\n")
        if unknown:
            print(f"Warning: preserving unrecognized assay labels: {unknown}", file=sys.stderr)
        print(json.dumps(summary, indent=2))
    except (ValueError, OSError, pd.errors.ParserError) as exc:
        parser.exit(1, f"Error: {exc}\n")

if __name__ == "__main__":
    main()
