# ATAC-seq × biological benchmark integration freeze

Freeze date: 2026-09-21

## Reference

- Assembly: GCF_003668045.3
- Benchmark: 31 Gold loci
  - 27 positive
  - 4 negative

## ATAC-seq dataset

- Study: PRJNA667472
- Context: CHO-K1 bulk ATAC-seq
- Runs:
  - SRR12774931
  - SRR12774932

## Reproducible accessibility

Primary definition:

- IDR <= 0.10
- 44,216 raw IDR records
- 39,332 exact unique genomic intervals

Sensitivity definition:

- IDR <= 0.05
- 36,800 raw IDR records
- 33,226 exact unique genomic intervals

After exact-coordinate deduplication, no remaining ATAC
intervals overlapped one another within either IDR set.

## Locus representation

Each Gold locus retains its original experimental interval and
also receives a canonical anchor defined as the integer centre of
the 1-based closed benchmark interval.

Primary fixed anchor windows:

- +/- 1 kb
- +/- 10 kb
- +/- 50 kb

This avoids using heterogeneous source-locus interval length as the
main representation for class comparisons.

## Continuous ATAC signal

MACS3 treatment pileup SPMR BigWig tracks were summarized
independently for both biological libraries.

Consensus signal is the arithmetic mean of replicate window means.

Pearson correlation between replicate window means across the
31 Gold loci:

- +/- 1 kb: 0.997614
- +/- 10 kb: 0.994910
- +/- 50 kb: 0.995451

These values describe concordance specifically across the benchmark
locus windows and are not genome-wide correlations.

## Exploratory biological pattern

Global median consensus SPMR:

- +/- 1 kb:
  - positive: 0.0689047
  - negative: 0.0758625

- +/- 10 kb:
  - positive: 0.0606521
  - negative: 0.0743960

- +/- 50 kb:
  - positive: 0.0820933
  - negative: 0.0617772

Dhiman-only sensitivity analysis:

- +/- 1 kb:
  - positive: 0.0698046
  - negative: 0.0758625

- +/- 10 kb:
  - positive: 0.0898570
  - negative: 0.0743960

- +/- 50 kb:
  - positive: 0.1030404
  - negative: 0.0617772

The observed pattern is scale-dependent and exploratory. It must not
be interpreted as a predictive or causal relationship because the
Gold set contains only four negative loci and the ATAC-seq layer
comes from a CHO-K1 context not experimentally matched to every
benchmark locus.

## Primary output

`benchmark_atac_features_primary.tsv`

Contains one row per Gold locus and the predefined primary ATAC
features used for later multi-omic integration.

## Sensitivity output

`benchmark_atac_features_idr05_sensitivity.tsv`

Contains the stricter IDR <= 0.05 accessibility representation.

## Reproducibility

Input and output checksums and provenance metadata are stored in this
directory and in the corresponding runtime results directory.

The Python analysis completed successfully.

The Snakemake wrapper rule exists in the repository but its dry-run
remains to be validated when a persistent Snakemake environment is
available on VERA.
