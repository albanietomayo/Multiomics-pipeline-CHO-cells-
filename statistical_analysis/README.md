# Sequencing sample statistics

Count distinct biological sample IDs per cell line and sequencing type using metadata. The self-contained Plotly HTML and vector PDF show the same grouped bar chart.

## Current results

The repository run used `config/samples.tsv`: 1,947 run records. Of these, 1,132 have complete required fields and represent 915 distinct biological samples across 49 cell-line labels (910 RNA, 5 ATAC). Deduplication removed 217 repeated sample-assay assignments.

**Coverage limitation:** 815 run records have missing cell lines: 686 RNA-seq, 120 ChIP-seq, and 9 ATAC-seq. These records cover 606 distinct sample IDs and are saved in `results/excluded_metadata.tsv`. All ChIP records are excluded, so no ChIP category can be assigned to a cell line. Counts describe catalogued metadata availability; local FASTQ download or processing completion is not assessed.

Cell-line spelling and case variants remain distinct, including clone names. No labels are inferred from study titles or merged into parent lines.

## Environment and PDF export

All analysis scripts run in `cho-multiomics`. Installed versions for this run: Python 3.14.7, pandas 3.0.5, Plotly 6.9.0, Kaleido 1.3.0. Installation also updated OpenSSL from 3.6.3 to 3.6.4.

```bash
conda install -n cho-multiomics -c conda-forge pandas plotly python-kaleido
```

By default PDF export uses Kaleido with Linux Chrome/Chromium. If unavailable, `conda run -n cho-multiomics plotly_get_chrome -y` installs it.

This machine already has **Windows Chrome**. The initial results use `--windows-chrome` to launch it from WSL and print the identical Plotly figure with a page sized to the chart. This avoids installing Linux Chrome. A separate temporary browser profile and print HTML are cleaned up after export; the personal browser profile is not used. The Windows path must be supplied as a WSL `/mnt/c/...` path. This alternative currently expects Chrome under the user's `AppData/Local/Google/Chrome/Application` installation. Export backend and browser path are recorded in `plot_provenance.json`.

## Reproduce this repository analysis

Run from the repository root in WSL:

```bash
export CONDA_EXE=/home/alba/miniconda3/bin/conda
"$CONDA_EXE" run -n cho-multiomics python statistical_analysis/scripts/prepare_repository_metadata.py \
  --metadata config/samples.tsv
bash statistical_analysis/scripts/run_analysis.sh \
  --metadata statistical_analysis/results/metadata_complete.tsv \
  --sample-id-col sample_accession --sequencing-type-col omics \
  --subtitle 'Annotated samples only; see metadata_coverage.json for exclusions' \
  --windows-chrome /mnt/c/Users/Usuario/AppData/Local/Google/Chrome/Application/chrome.exe
```

Preparation writes a complete-label subset, the excluded source rows with reasons, and a coverage report including the source SHA256 and command. It does not modify source metadata. The collector itself remains strict: missing required values are errors.

## Generic metadata interface

```bash
bash statistical_analysis/scripts/run_analysis.sh --metadata path/to/metadata.tsv \
  --sample-id-col sample_id --cell-line-col cell_line \
  --sequencing-type-col sequencing_type --output-dir statistical_analysis/results
```

CSV and TSV are supported. The three column names shown are defaults; `--metadata` is required. Output defaults to `statistical_analysis/results` relative to the scripts. The wrapper runs both Python scripts with `conda run -n cho-multiomics`; set `CONDA_EXE` if conda is absent from PATH. Pass `--windows-chrome` as above on this machine. `--subtitle` changes the plot annotation.

To regenerate only plots:

```bash
conda run -n cho-multiomics python statistical_analysis/scripts/plot_sample_counts.py \
  --counts statistical_analysis/results/counts_by_cell_line_and_assay.csv \
  --windows-chrome /mnt/c/Users/Usuario/AppData/Local/Google/Chrome/Application/chrome.exe
```

## Counting rules

- Trim required fields; reject empty/NA/N/A/null/none/nan markers and conflicting cell lines for one sample ID.
- Normalize RNA, RNA-seq, scRNA/scRNA-seq to RNA; ATAC equivalents to ATAC, case-insensitively.
- Expand Multiome and RNA+ATAC/ATAC+RNA labels into one assignment per assay.
- Preserve other assay labels, such as ChIP-seq, with a warning.
- Deduplicate sample ID × cell line × normalized assay. Read pairs, lanes, and repeated run records sharing a sample ID count once.
- Include zero-valued combinations of every observed cell line and observed assay in the accepted metadata.

Totals across assays count sample-assay assignments and can exceed unique biological samples. Archive sample IDs are the available proxy for biological sample identity; separate accessions are not inferred to be the same replicate.

## Outputs

- `sample_assay_inventory.csv`: deduplicated sample-level audit table.
- `counts_by_cell_line_and_assay.csv`: tidy counts; `counts_pivot.csv`: wide table.
- `validation_summary.txt`: input SHA256, command, dependency versions, duplicate removals, unknown assays, and totals.
- `sample_counts_by_cell_line_and_assay.html`: offline interactive chart, with hover and legend controls.
- `sample_counts_by_cell_line_and_assay.pdf`: vector figure; intentionally wide to accommodate 49 labels.
- `plot_provenance.json`: export command, package versions, browser/backend.
- `metadata_complete.tsv`, `excluded_metadata.tsv`, `metadata_coverage.json`: preparation audit.
- `verification.json`: real-data table/figure consistency checks.

Use zoom or horizontal scrolling for the wide chart. Export failure returns nonzero and preserves tables/HTML. Use a fresh output directory for independent analyses; successful runs replace same-named generated files.

## Tests

```bash
conda run -n cho-multiomics python -m unittest discover -s statistical_analysis/tests -v
```

Tests cover technical duplicates, Multiome expansion, zero combinations, unknown assays, missing fields, conflicting assignments, empty input, aliases, and exact equality of plotted trace values and count tables. `tests/fixture.csv` supports a complete wrapper smoke test; its outputs are in `results/test_fixture/`.
