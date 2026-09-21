# RNA-seq protocol curation freeze

Date: 2026-09-21

## Purpose

Protocol-level re-audit performed before construction of the
primary CHO bulk RNA-seq expression atlas.

Public-repository records labelled RNA-Seq may correspond to
specialized transcriptomic assays that are not quantitatively
comparable with conventional whole-transcriptome bulk RNA-seq.

## Newly curated exclusions

Thirteen previously eligible runs were assigned:

protocol_status = other_assay

which produces:

selection_reason = specialized_RNA_assay

Breakdown:

- PRJNA575354
  - 11 runs
  - metabolic-labeling / RNA-stability experimental design

- PRJNA378939
  - 2 runs
  - RNA-PET

## Eligibility effect

Before re-curation:

- selection_gate eligible RNA runs: 1753

After re-curation:

- selection_gate eligible RNA runs: 1740

Exact change:

- 13 eligible -> excluded
- 0 blocked/excluded -> eligible
- no unintended eligibility changes

## Current catalogue

RNA records: 1849

Selection-policy level:

- retained_by_rules: 1748
- excluded: 93
- review_required: 8

Among the 1748 retained_by_rules records, 8 use sequencing
platforms unsupported by the production RNA pipeline.

Therefore:

1748 retained_by_rules
- 8 unsupported platform
= 1740 selection_gate-eligible RNA runs

## Gate exclusions

- non_CHO: 23
- single_cell_protocol: 35
- specialized_RNA_assay: 35
- protocol_review: 8
- unsupported_platform: 8

## Effect on expression QC

Before protocol re-curation, the extreme low-count tail contained
specialized assays.

After removal of the 13 newly curated specialized assays in the
available production snapshot:

- minimum assigned gene counts: 210955
- minimum non-zero gene fraction: 0.332457
- retained runs below 100000 assigned counts: 0
- retained runs below 0.30 non-zero gene fraction: 0

No arbitrary quantitative run-level expression threshold was
introduced.

## Experimental unit

Individual sequencing runs are not treated as independent biological
expression units when multiple runs belong to the same ENA experiment.

Planned RNA integration hierarchy:

run-level raw counts
    -> sum technical runs within experiment_accession
    -> experiment-level expression QC
    -> normalization
    -> robust within-study summary
    -> balanced cross-study CHO transcriptomic consensus
    -> projection onto benchmark loci
