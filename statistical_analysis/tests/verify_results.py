"""Cross-check emitted tables, plotted data and offline/export artifacts."""
import json
from pathlib import Path
import sys
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from plot_sample_counts import make_figure

def verify(directory):
    inventory = pd.read_csv(directory / "sample_assay_inventory.csv", dtype=str, keep_default_na=False)
    counts = pd.read_csv(directory / "counts_by_cell_line_and_assay.csv", dtype={"cell_line": str}, keep_default_na=False)
    pivot = pd.read_csv(directory / "counts_pivot.csv", index_col=0, keep_default_na=False)
    expected = inventory.groupby(["cell_line", "sequencing_type"]).sample_id.nunique()
    for row in counts.itertuples():
        assert row.count == expected.get((row.cell_line, row.sequencing_type), 0)
        assert row.count == pivot.loc[row.cell_line, row.sequencing_type]
    assert counts["count"].sum() == len(inventory)
    fig = make_figure(counts)
    plotted = {(x, t.name): int(y) for t in fig.data for x, y in zip(t.x, t.y)}
    assert plotted == {(r.cell_line, r.sequencing_type): r.count for r in counts.itertuples()}
    html = (directory / "sample_counts_by_cell_line_and_assay.html").read_text()
    assert "Plotly.newPlot" in html and "plotly.js" in html
    assert '<script src="https://cdn.plot.ly' not in html
    pdf = directory / "sample_counts_by_cell_line_and_assay.pdf"
    assert pdf.read_bytes().startswith(b"%PDF-") and pdf.stat().st_size > 1000
    result = {"status": "passed", "unique_samples": inventory.sample_id.nunique(),
              "sample_assay_assignments": len(inventory), "grid_rows": len(counts),
              "checks": ["counts match distinct inventory", "pivot matches tidy table", "sum equals inventory",
                         "plot traces match every count", "HTML includes Plotly inline", "PDF signature and size"]}
    (directory / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    verify(Path(sys.argv[1]))
