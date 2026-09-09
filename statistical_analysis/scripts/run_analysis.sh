#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/../results"
SUBTITLE="Distinct biological samples; Multiome contributes to both RNA and ATAC"
ARGS=()
PLOT_ARGS=()
while (($#)); do
    case "$1" in
        --output-dir) OUTPUT_DIR="${2:?Missing output directory}"; shift 2 ;;
        --subtitle) SUBTITLE="${2:?Missing subtitle}"; shift 2 ;;
        --windows-chrome) PLOT_ARGS+=("--windows-chrome" "${2:?Missing Chrome path}"); shift 2 ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
CONDA_BIN="${CONDA_EXE:-conda}"
"$CONDA_BIN" run --no-capture-output -n cho-multiomics python "$SCRIPT_DIR/collect_sample_counts.py" "${ARGS[@]}" --output-dir "$OUTPUT_DIR"
"$CONDA_BIN" run --no-capture-output -n cho-multiomics python "$SCRIPT_DIR/plot_sample_counts.py" --counts "$OUTPUT_DIR/counts_by_cell_line_and_assay.csv" --output-dir "$OUTPUT_DIR" --subtitle "$SUBTITLE" "${PLOT_ARGS[@]}"
