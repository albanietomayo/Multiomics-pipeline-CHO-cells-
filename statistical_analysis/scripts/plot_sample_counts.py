#!/usr/bin/env python3
"""Export one grouped Plotly bar chart to standalone HTML and PDF."""
import argparse
import json
from pathlib import Path
import shlex
import sys
import pandas as pd
import plotly.express as px
from collect_sample_counts import DEFAULT_OUTPUT, versions

def make_figure(counts, subtitle=""):
    lines = sorted(counts.cell_line.unique())
    assays = sorted(counts.sequencing_type.unique())
    fig = px.bar(counts, x="cell_line", y="count", color="sequencing_type", barmode="group", text="count",
                 category_orders={"cell_line": lines, "sequencing_type": assays},
                 labels={"cell_line": "Cell line (metadata label)", "count": "Biological samples", "sequencing_type": "Sequencing type"},
                 title="Sequencing samples by cell line and assay" + (f"<br><sup>{subtitle}</sup>" if subtitle else ""),
                 color_discrete_sequence=px.colors.qualitative.Safe)
    fig.update_traces(textposition="outside", cliponaxis=False, hovertemplate="Cell line: %{x}<br>Biological samples: %{y}<extra>%{fullData.name}</extra>")
    fig.update_layout(template="plotly_white", width=max(1000, 80 * len(lines)), height=850,
                      margin=dict(l=90, r=50, t=130, b=310), legend=dict(orientation="h", y=1.12, x=0), font=dict(size=14), bargap=0.22)
    fig.update_xaxes(tickangle=-55, automargin=True)
    fig.update_yaxes(rangemode="tozero", tickformat="d", range=[0, max(1, counts["count"].max()) * 1.18])
    return fig

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, default=DEFAULT_OUTPUT / "counts_by_cell_line_and_assay.csv")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--subtitle", default="Distinct biological samples; Multiome contributes to both RNA and ATAC")
    parser.add_argument("--windows-chrome", type=Path, help="Existing Windows chrome.exe exposed through /mnt/c (WSL PDF alternative)")
    args = parser.parse_args()
    try:
        counts = pd.read_csv(args.counts, dtype={"cell_line": str, "sequencing_type": str}, keep_default_na=False)
        if not {"cell_line", "sequencing_type", "count"}.issubset(counts.columns) or counts.empty:
            raise ValueError("Counts require nonempty cell_line, sequencing_type, and count columns.")
        counts["count"] = pd.to_numeric(counts["count"], errors="raise")
        if counts["count"].isna().any() or ((counts["count"] < 0) | (counts["count"] % 1 != 0)).any():
            raise ValueError("Counts must be finite nonnegative integers.")
        if counts.duplicated(["cell_line", "sequencing_type"]).any():
            raise ValueError("Duplicate cell-line/assay combinations.")
        if (counts[["cell_line", "sequencing_type"]].apply(lambda s: s.str.strip().eq("")).any().any()):
            raise ValueError("Cell-line and assay labels cannot be blank.")
        fig = make_figure(counts, args.subtitle)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        stem = args.output_dir / "sample_counts_by_cell_line_and_assay"
        fig.write_html(str(stem) + ".html", include_plotlyjs=True, full_html=True)
        if args.windows_chrome:
            from windows_chrome_export import export_pdf
            backend = export_pdf(fig, str(stem) + ".pdf", args.windows_chrome)
        else:
            fig.write_image(str(stem) + ".pdf")
            backend = {"backend": "Kaleido"}
        (args.output_dir / "plot_provenance.json").write_text(json.dumps({"command": shlex.join([sys.executable, *sys.argv]), "versions": versions(), **backend}, indent=2) + "\n")
        print(f"Created {stem}.html and {stem}.pdf")
    except Exception as exc:
        parser.exit(1, f"Plot/export failed: {exc}\nPDF requires Chrome/Chromium for Kaleido, or --windows-chrome for existing Windows Chrome from WSL.\n")

if __name__ == "__main__":
    main()
