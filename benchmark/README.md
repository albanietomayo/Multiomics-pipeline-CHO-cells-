# CHO genomic integration benchmark

Independent benchmark for future predictive genome-engineering models in CHO cells.

## Canonical reference

- Assembly: CriGri-PICRH-1.0
- RefSeq: GCF_003668045.3
- GenBank: GCA_003668045.2

All harmonized benchmark coordinates must ultimately refer to this assembly.

## Coordinate conventions

- Source coordinates are preserved exactly as reported and accompanied by their source coordinate convention.
- Curated target coordinates use 1-based closed coordinates.
- BED outputs use 0-based half-open coordinates.

## Evidence policy

The gold benchmark contains experimentally observed integration outcomes.

Prediction-only regions derived from RNA-seq, ChIP-seq, Hi-C or other omics analyses are retained only as supporting evidence and are not promoted to gold positive or negative labels.

An untested genomic region is never considered a negative benchmark locus.

## Data model

`benchmark_observations.tsv`
contains one row per published experimental observation.

`benchmark_loci.tsv`
is a derived table with one row per harmonized genomic locus.

The original published coordinates are never overwritten.
