# Downstream multi-omic exploration v1

This directory contains the final descriptive downstream analysis used to close the TFM results section.

## Input

`benchmarks/multiomics/2026-09-23/feature_views_v1/benchmark_37_multiomic_core_feature_view_v1.tsv`

The frozen Core view contains 37 loci x 69 columns: `canonical_locus_id`, `benchmark_role`, and 67 model-facing predictors.

## Analysis 1: within-study sensitivity (Dhiman et al.)

The four negative loci all originate from the Dhiman study. To reduce the possibility that the global positive/negative contrast is driven only by study provenance, a descriptive sensitivity analysis was restricted to the 13 Dhiman loci:

- 6 positive
- 4 negative
- 3 support_only

The positive-vs-negative comparison is descriptive only. No inferential tests, supervised feature selection, classifier, AUC, or predictive performance estimate is produced.

Output: `dhiman_within_study_summary_v1.tsv`.

## Analysis 2: Figure 5 multi-omic heatmap

Figure 5 uses 24 biologically interpretable variables selected a priori from the Core view:

- 4 RNA-seq features
- 8 ATAC-seq features
- 12 ChIP-seq features: two per histone mark for H3K27ac, H3K4me3, H3K36me3, H3K4me1, H3K27me3, and H3K9me3

Each feature is standardized independently as a z-score for visualization only. Loci are ordered by unsupervised hierarchical clustering using Ward linkage and Euclidean distance. `benchmark_role` is not used in feature selection, scaling, distance calculation, or clustering; it is shown only as an annotation (`P`, `N`, `S`).

Outputs:

- `Figure_5_heatmap_multiomic_37_loci.png`
- `Figure_5_heatmap_multiomic_37_loci.pdf`
- `Figure_5_heatmap_matrix_estandarizada.tsv`
- `SHA256SUMS.txt`

## Reproducibility

Run from repository root:

```bash
conda env create -f workflow/envs/multiomics_downstream.yaml
conda run -n multiomics_downstream python workflow/scripts/summarize_multiomic_downstream.py \
  --core benchmarks/multiomics/2026-09-23/feature_views_v1/benchmark_37_multiomic_core_feature_view_v1.tsv \
  --outdir benchmarks/multiomics/2026-09-24/downstream_exploration_v1
```

This downstream analysis does not modify the frozen RNA, ATAC, ChIP, master, feature-manifest, Core, or Extended products.
