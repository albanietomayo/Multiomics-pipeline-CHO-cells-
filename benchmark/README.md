# Canonical CHO integration-site benchmark inputs

`curated/` contains byte-for-byte copies of the external curated benchmark at
`/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_benchmark/benchmark/curated/`.
`benchmark_loci.tsv` and `.bed` define all 37 loci: 27 positive, four negative,
and six `support_only`. `benchmark_observations.tsv` retains the observation-level
curation and evidence behind the canonical rows. The `benchmark_gold_binary`
files are the historical 31-locus positive/negative subset used by the legacy
ATAC integration rule; they are not the canonical 37-locus benchmark.

`derived/benchmark_37_nearest_gene_tss_audit.tsv` is a byte-for-byte copy of
the 37-row locus semantic-class and nearest-gene TSS audit used by the frozen
RNA and ATAC locus-context producers. Its source was
`/cephyr/users/mayoa/Vera/hpc_rnaseq_production_results/2026-09-10/reference_gtf_recovery_audit/benchmark_37_nearest_gene_tss_audit.tsv`.

`derived/chipseq_source_metadata/` contains the small condition metadata,
integration manifest and production summary named as inputs in
`benchmarks/chipseq/2026-09-23/integration_input_v1/provenance.json`. Their
original source is the separate `TFM_multiomics_pipeline_chipseq_generalization`
project under `snapshots/chipseq/`. The manifest refers to external peak and
signal files; copying it does not make those heavy artifacts available or make
ChIP processing runnable within this repository alone.

SHA-256 of imported files (also recorded in the relevant frozen provenance):

| File | SHA-256 |
| --- | --- |
| `curated/benchmark_loci.tsv` | `ea5a30f96853b288791e12ca84ed7c438987ed7cc9605d81fe503cc52438bdf0` |
| `curated/benchmark_loci.bed` | `fabcb8dbc2dfdb79c9a3c6d3e518b7c7a25a63f1137cfc1d7e3028b770261d07` |
| `curated/benchmark_gold_binary.tsv` | `7274efa875783d5545d952d74b155d6bf7c61eb65c6a18d04f74210e062c0203` |
| `curated/benchmark_gold_binary.bed` | `8e3d3999cab1a83deeef0c6996b6ebd79bf3c9d91aa36a2a30c3ec7d97844093` |
| `curated/benchmark_observations.tsv` | `115575d7b014ba18ab97f34f3bf9902d80ab5673356aeb44f27b8743241c62dc` |
| `derived/benchmark_37_nearest_gene_tss_audit.tsv` | `6fbf88c612b02d76a46deb7dd6e9913f89276e12e55e91e8cfb782df05f52cdf` |
| `derived/chipseq_source_metadata/chipseq_conditions.tsv` | `12a6cb7c77a9ec9e1d6f560d24e6fa1e7ac5b1ecfa5192d668134cca6a7a447c` |
| `derived/chipseq_source_metadata/chipseq_integration_manifest.tsv` | `fd36f996d618ec623203bc7107f03116b2c68486010ac07ece89874e4a8e2a65` |
| `derived/chipseq_source_metadata/chipseq_production_summary.tsv` | `898cbf4b20e01c7574a82825881f1d73ddfb95dcd6c054a3528213222b04670c` |

Historical provenance paths are retained unchanged in frozen JSON. The local
copies are portable reference inputs, not a rewritten execution record.
