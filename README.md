# Building training and benchmarking data pipelines for predictive genome engineering in CHO cells

This repository assembles auditable genomic context for integration-site
research in Chinese hamster ovary (CHO) cells. It combines **RNA-seq**
transcriptional context, **ATAC-seq** chromatin accessibility and **ChIP-seq**
histone-mark context at a curated 37-locus benchmark. The main reusable products
are frozen tables with QC, provenance and checksums. No other omics modality is
represented as integrated here.

## Scientific motivation and scope

Genome engineering studies need consistent reference coordinates, traceable
integration-site evidence and comparable molecular context. This project builds
the data infrastructure for that work. It does not establish a validated
classifier or demonstrate generalization to unseen integration sites. The
benchmark has only four negative loci, so current outcome comparisons are
exploratory and descriptive.

```text
Public sequencing data
        |
  RNA / ATAC / ChIP
        |
37 canonical benchmark loci
        |
Multiomic master: 37 x 588
        |
Feature governance: 317 biological candidates
        |
-150 deterministic redundancies
        |
Extended: 167 predictors / Core: 67 predictors
```

## What this repository provides

- A repository-local copy of the curated 37-locus benchmark, its observation
  evidence, the historical 31-locus binary subset, and small derived audit
  inputs ([benchmark provenance](benchmark/README.md)).
- Frozen locus-level RNA, ATAC and ChIP contexts, a 37 × 588 multiomic master,
  a feature manifest, and Core/Extended predictor views.
- Snakemake rules for metadata, reference, RNA-seq and ATAC-seq processing,
  plus standalone validated scripts for final locus and feature integration.
- Historical execution evidence and descriptive downstream figures and tables.

## Key validated outputs

| Product | Location | Shape or purpose |
| --- | --- | --- |
| Multiomic master | [master context](benchmarks/multiomics/2026-09-23/master_context_v1/benchmark_37_multiomic_master_context_v1.tsv) | 37 loci × 588 columns |
| Feature manifest | [governance table](benchmarks/multiomics/2026-09-23/feature_manifest_v1/multiomic_feature_manifest_v1.tsv) | inclusion and exclusion policy |
| Extended view | [extended table](benchmarks/multiomics/2026-09-23/feature_views_v1/benchmark_37_multiomic_extended_feature_view_v1.tsv) | 37 × 169: ID, label, 167 predictors |
| Core view | [core table](benchmarks/multiomics/2026-09-23/feature_views_v1/benchmark_37_multiomic_core_feature_view_v1.tsv) | 37 × 69: ID, label, 67 predictors |
| RNA context | [RNA table](benchmarks/rnaseq/2026-09-23/rna_locus_context_v1/benchmark_37_rna_context_v1.tsv) | gene-level transcriptional context |
| ATAC context | [ATAC table](benchmarks/atacseq/2026-09-23/atac_locus_context_v1/benchmark_37_atac_context_v1.tsv) | accessibility and signal context |
| ChIP context | [ChIP table](benchmarks/chipseq/2026-09-23/locus_context_v1/benchmark_37_chipseq_context_v1.tsv) | six histone-mark blocks |

Each frozen product directory has a QC table, provenance JSON and
`SHA256SUMS.txt`. Read these with the data: column names alone do not capture
the input restrictions or interpretation limits.

## Quick start

### Path A: downstream reuse

Clone the repository and inspect the frozen tables directly. No FASTQ download
is needed. For example, from the repository root:

```bash
python - <<'PY'
import csv
from collections import Counter

path = 'benchmarks/multiomics/2026-09-23/feature_views_v1/benchmark_37_multiomic_core_feature_view_v1.tsv'
with open(path, newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    rows = list(reader)
    print('loci:', len(rows), 'columns:', len(reader.fieldnames))
    print('roles:', dict(Counter(row['benchmark_role'] for row in rows)))
PY
```

`canonical_locus_id` identifies rows and `benchmark_role` is a label only.
Neither is a predictor. Keep `support_only` separate from `negative`; do not
recode it as a negative or use it as evidence of failure. Review the
[feature specification](benchmarks/multiomics/2026-09-23/feature_views_v1/multiomic_feature_view_spec_v1.tsv)
before choosing a view. The master is an audit-rich context table, not a ready
predictor matrix.

### Path B: workflow reproduction

The root [environment](environment.yml) pins the Snakemake/Python environment.
With Conda installed, the starting commands are:

```bash
conda env create -f environment.yml
conda run -n cho-multiomics snakemake -n --cores 1 --default-resources tmpdir=/tmp
```

For actual rule execution, enable the pinned rule environments in
`workflow/envs/` with `--use-conda`; the root environment alone does not contain
all aligners, QC tools and track-processing packages. Use the same flag when
dry-running a target whose rule environments you intend to inspect. The
standalone integration scripts also have input and package requirements that
are not assembled into a single final Snakemake target.

This dry-run may require access to metadata or inputs absent from a fresh clone.
Do not infer full reproducibility from a successful dry-run. The root
[Snakefile](Snakefile) includes metadata, reference, RNA-seq and ATAC-seq rules.
Its default `all` target covers metadata, selected FASTQ/QC, reference and
RNA-seq validation outputs; ATAC production has a separate `atacseq_production`
target. ChIP production is not included in this Snakefile. The final 37-locus
contexts, master and feature views were produced using validated standalone
scripts under [workflow/scripts](workflow/scripts), with frozen execution
products under `benchmarks/`. There is currently **no single-command
raw-to-final multiomic reconstruction** in this repository.

## Repository structure and data sources

- `benchmark/`: curated benchmark and small source audits.
- `benchmarks/`: frozen scientific products and historical validation evidence.
- `config/`: metadata curation decisions, run selection and workflow settings.
- `workflow/rules/`, `workflow/scripts/`, `workflow/envs/`: Snakemake rules,
  standalone producers and component environments.
- `docs/`: [ATAC methods/status](docs/ATAC_WORKFLOW.md) and
  [HPC execution notes](docs/HPC_EXECUTION.md).
- `statistical_analysis/`: [sample catalogue analysis](statistical_analysis/README.md).

The raw sequencing records were retrieved from public archives using the
metadata and FASTQ acquisition scripts. The archived metadata snapshots under
`snapshots/` document selected curation states. FASTQs, BAMs, indexes, BigWigs
and other large production inputs are not distributed in Git.

## Reference genome and metadata curation

All integrated locus contexts use *Cricetulus griseus* CriGri-PICRH-1.0,
RefSeq **GCF_003668045.3**. The mitochondrial accession is **NC_007936.1**.
Reference retrieval and annotation paths are configured in
[config/config.yaml](config/config.yaml). Do not combine locus coordinates from
other assemblies without explicit remapping and an audit.

The [metadata rules](workflow/rules/metadata.smk) query public run metadata,
classify assay and CHO relevance, apply curated study/run decisions and build
FASTQ manifests. Selection is auditable in `config/` and historical `snapshots/`.
Catalogue presence alone does not establish eligibility for molecular analysis.
Network access and archive availability are required for fresh retrieval.

## Modality workflows

**RNA-seq.** [RNA rules](workflow/rules/rnaseq.smk) cover alignment, QC,
strandedness and gene counting. The frozen downstream RNA products add
normalization and study-level aggregation before nearest-gene TSS association
at the 37 loci. Gene expression is a gene-level proxy, not a direct local
transcription measurement at an integration coordinate. Some normalization and
locus-context steps were standalone, not part of the default Snakemake target.

**ATAC-seq.** [ATAC rules](workflow/rules/atacseq.smk) have a validated SINGLE
bulk CHO production path with filtering, Tn5 transformation, TSS QC, MACS3
peaks, FRiP, BigWig signal and IDR reproducibility. The eligible bulk runs are
SRR12774931 and SRR12774932. PAIRED processing has technical tests only; no
PAIRED biological validation is claimed. See the [ATAC workflow notes](docs/ATAC_WORKFLOW.md)
and [production capsule](benchmarks/atacseq/2026-09-19/production_10318071_capsule/README.md).
ATAC production is a separate target and its heavy inputs are external.

**ChIP-seq.** The frozen integration products summarize six histone marks:
H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me3 and H3K9me3. Production
alignment/peak processing was performed in a separate ChIP project; this
repository contains [small source metadata](benchmark/README.md), integration
scripts and frozen context outputs, but no complete local ChIP production DAG.
The source manifest names external peak and signal artifacts that are not in Git.

## Canonical integration-site benchmark

[benchmark/curated/benchmark_loci.tsv](benchmark/curated/benchmark_loci.tsv)
defines **37 loci: 27 positive, four negative and six `support_only`**.
[Observation-level evidence](benchmark/curated/benchmark_observations.tsv) and
the [locus audit](benchmark/derived/benchmark_37_nearest_gene_tss_audit.tsv)
document curation, coordinate semantics and TSS context. The Gold31 files are
the earlier 27-positive/4-negative binary subset. The legacy
[integration rule](workflow/rules/integration.smk) and
[configuration](config/integration.yaml) consume Gold31 and preserve original
VERA paths; the root Snakefile does not include that rule. Gold31 must not be
presented as the final 37-locus integration.

## Multiomic integration and feature governance

Validated standalone producers created the 37-locus RNA, ATAC and ChIP
contexts, then the 37 × 588 [master](benchmarks/multiomics/2026-09-23/master_context_v1/).
The [manifest](benchmarks/multiomics/2026-09-23/feature_manifest_v1/) identifies
**317 biological candidate features**. Deterministic, outcome-blind redundancy
removal discarded **150**, leaving **167 Extended predictors**. The Core view
retains **67**. Feature inclusion/reduction did not use `benchmark_role` or
benchmark outcomes. Coordinates, evidence metadata and high-cardinality
identifiers are not predictors. The 96 ChIP percentile columns defined relative
to these 37 benchmark loci are composition-dependent and are not portable
baseline predictors.

## Downstream exploratory analysis

The [downstream exploration](benchmarks/multiomics/2026-09-24/downstream_exploration_v1/README.md)
contains descriptive within-study sensitivity and a multiomic heatmap. The
small, uneven benchmark, especially its four negatives, cannot support robust
supervised performance or generalization claims. `support_only` remains a
distinct class throughout.

## Reproducibility, computation and availability

Every frozen 2026-09-23 context, master and feature product includes QC,
provenance and a SHA-256 manifest. The JSON records original execution paths,
including historical absolute cluster paths; those paths are part of the audit
record and have not been rewritten. Verify checksums from each product directory
with `sha256sum -c SHA256SUMS.txt`. The local benchmark copies and their source
hashes are documented in [benchmark/README.md](benchmark/README.md).

Downstream table reuse needs only ordinary CSV/TSV tools and modest storage.
Upstream sequencing reproduction needs Conda, Snakemake, public archive access,
the reference genome, substantial scratch/disk capacity and compute suitable
for alignment, peak calling and signal processing. The
[SLURM examples](docs/HPC_EXECUTION.md) target VERA/C3SE; other sites must adapt
accounts, partitions and scratch paths. External ChIP and ATAC production
artifacts and raw sequencing data are not bundled. Local benchmark inputs make
the frozen context easier to audit, but do not supply all heavy inputs required
to rerun its producers.

### Tests and environment limits

Create the root environment with `conda env create -f environment.yml` before
running Python tests. From the repository root:

```bash
conda run -n cho-multiomics python -m unittest discover -s statistical_analysis/tests -v
conda run -n cho-multiomics python -m unittest discover -s workflow/tests -v
```

The statistical-analysis suite additionally needs Plotly and Kaleido; see its
[environment instructions](statistical_analysis/README.md#environment-and-pdf-export).
The full workflow suite also exercises packages from component environments
(including `pysam` and `pyBigWig`), Snakemake DAG checks and an HPC inventory
test that queries `squeue`. A local shell lacking those packages or SLURM may
run only a subset. Test failures caused by missing tools should be distinguished
from selection or scientific-output regressions; no SLURM job submission is
needed for the documented tests.

## Citation and archived release

The exact archived **v1.0.0** release has the
[version-specific DOI](https://doi.org/10.5281/zenodo.22957205).
The [concept DOI](https://doi.org/10.5281/zenodo.22957204) resolves to the
software record across versions. Cite the version-specific DOI when referring
to v1.0.0; [CITATION.cff](CITATION.cff) contains its citation metadata.

## Known limitations and release status

- The four negative loci limit outcome-based inference; no robust supervised
  predictive performance is claimed.
- Some upstream input data and ChIP production workflows live outside this
  repository. The final multiomic integration is a validated sequence of
  standalone scripts and frozen products, not a unified raw-to-final DAG.
- Historical VERA paths in provenance and site-specific scripts describe real
  executions. They are not portable defaults.
- Capability status: **downstream reuse: YES; workflow reproduction: PARTIAL;
  full raw-to-final reproduction: NO, currently**.
- Original repository source code is under the [MIT License](LICENSE).
- The v1.0.0 citation and archive identifiers are recorded above and in
  [CITATION.cff](CITATION.cff). Third-party data retain their source terms.

The v1.0.0 archive records this research repository as released. External-data
packaging and the partial upstream workflow boundary remain as described above.
