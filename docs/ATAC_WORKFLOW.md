# ATAC-seq production module

## Scope and eligibility

The tracked catalogue contains 14 ATAC runs (17 FASTQs), but catalogue presence
is not production eligibility. `atacseq_production_manifest` joins the complete
FASTQ manifest with `config/samples.tsv` and recomputes the shared
`selection_policy.py`/`selection_gate.py` decision. Its file-level output is
`results/atacseq/metadata/production_manifest.tsv`; it preserves accessions,
original and processing layouts, FASTQ structure/role and the selection reason.
It fails closed for absent, ambiguous or inconsistent metadata. The validation
FASTQ manifest is not consulted by production ATAC.

The current eligible bulk CHO cohort is `SRR12774931` and `SRR12774932`.
`SRR12774934` is a non-CHO historical technical pilot. `SRR29929613`,
`SRR29929615` and `SRR29929617` are CHO single-cell/multiome libraries. These
and the remaining non-CHO tissues are classified in the production manifest but
cannot be forced through a downstream path: every run-specific rule depends on
the eligibility gate. `atacseq_production` is the explicit aggregate target; it
is intentionally absent from the global `all` target.

## Reference and validated core

The reference is *Cricetulus griseus* CriGri-PICRH-1.0, RefSeq
GCF_003668045.3, with mitochondrial accession NC_007936.1. Bowtie2 2.5.5,
samtools 1.24, Picard 3.5.0, pysam 0.24.0 and numpy 2.5.2 are pinned. Picard
marks duplicates with `REMOVE_DUPLICATES=false`; filtering subsequently excludes
them. Picard and paired-name sorting use the Snakemake temporary resource.

SINGLE filtering uses MAPQ 30 and excludes flag mask 3844. PAIRED filtering
requires mask 3, excludes 3852, and retains only exactly two records with one R1,
one R2 and reciprocal mate coordinates on the same reference. Secondary,
supplementary, duplicate, QC-fail, unmapped/mate-unmapped, mitochondrial and
low-MAPQ records are removed. No nuclear reference, no surviving record, an
orphan or an inconsistent pair is fatal.

Status: real bulk SINGLE CHO alignment/filtering/TSS behavior is **VALIDATED**.
PAIRED filtering is **TECHNICALLY TESTED** only; no bulk CHO PAIRED library is
currently eligible or biologically validated.

## Tn5 and TSS

The versioned biological Tn5 correction is forward `+4 bp`, reverse `-5 bp`.
The same insertion calculation is imported by TSS and the transformation script.
For the current SINGLE path, two deterministic outputs are explicit:

- `insertions.bed`: canonical 1-bp Tn5 insertion events;
- `shifted_read_intervals.bed`: original read intervals with only their 5-prime
  end shifted, matching the validated pilot input used by MACS3 and FRiP.

Coordinates and references are validated and input/output/drop counts plus
SHA-256 hashes are recorded. Both interval products are deterministically
coordinate-sorted with an external merge sort, including the pooled tagAlign;
this remains bounded-memory for production-sized inputs. PAIRED downstream
transformation is closed.

TSS enrichment uses ±2000 bp, 10-bp bins, 100-bp background on each side,
MAPQ 30 and the Tn5 shifts above. Boundary TSS are reported and omitted because
a complete window is required. Historical resource counts are 56,804 transcript
features, 42,215 unique TSS and 42,148 usable TSS. Bulk regression evidence is
9.693333 (`SRR12774931`) and 9.405718 (`SRR12774932`). These are observations,
not universal quality thresholds.

## Peaks, FRiP and tracks

MACS3 3.0.4 receives the shifted read intervals as `-f BED`, with `--nomodel`,
`--keep-dup all`, `--call-summits`, `-B`, `--SPMR`, `-p 0.01` and a semantic
`smooth_window_bp: 150`, from which `--extsize 150 --shift -75` is derived.
The latter parameters are signal extension/smoothing and are not Tn5 correction.
`--keep-dup all` is intentional because duplicate exclusion already happened.
Genome size is the sum of actual non-mitochondrial FAI reference lengths, recorded
as a nuclear-reference-span approximation rather than a mappability-derived
effective genome size. Pilot peak counts (not rerun here) are 110,311 and 215,037.

FRiP preserves the pilot denominator: usable shifted read intervals. Its numerator
is intervals overlapping at least one per-run narrowPeak. It records the exact
files and hashes; no threshold is imposed. Pilot values are 0.28102314 and
0.30694866.

BigWig is generated from the MACS3 SPMR treatment pileup (150-bp extension), not
literal 1-bp insertion coverage. bedGraph coordinates are clipped to the nuclear
reference (`max(0,start)`, `min(end,length)`); empty intervals are dropped and
unknown chromosomes are fatal. Metrics record input, kept, clipped, dropped,
unknown and clipped bases. The BigWig header may be a proper subset of the
reference, but it must equal exactly the chromosomes carrying retained signal in
the clipped bedGraph; outside-reference chromosomes, missing signal chromosomes
and length mismatches are fatal. The file receives a streaming SHA-256. Pilot
clipping evidence is 107/62/0 and 117/74/0 (clipped/dropped/unknown) for the two
runs.

## Replicates and reproducibility

`config/atacseq_replicate_groups.tsv` is the auditable source of grouping. The
current neutral group `PRJNA667472_CHO-K1_bulk_ATAC` contains the two eligible
comparable runs; it makes no biological-versus-technical replicate claim. A group
must have at least two unique, compatible, eligible SINGLE members. The DAG makes
a pooled tagAlign, pooled MACS3 peaks and IDR 2.0.4.2 results with seed 0. Both
IDR ≤0.10 and ≤0.05 outputs are retained. A final group provenance JSON
records membership, parameters and streaming hashes of the pooled, peak and IDR
artifacts. Pilot counts are 50,054 and 42,576.

## Provenance, environments and limitations

Per-run provenance records selection, source FASTQ identity, reference hash,
parameters, representations, nuclear span and hashes of persistent outputs.
Large-file SHA-256 values are computed in bounded-memory chunks. Paths are
project-relative. Direct production dependencies are pinned in separate
alignment/Picard/QC, downstream MACS3/bedtools, IDR and track environments.

MACS3, FRiP, IDR and BigWig are **PILOT-VALIDATED BUT NEWLY INTEGRATED**. The new
Snakemake downstream implementation is **NOT YET PRODUCTION-VALIDATED** until an
authorized limited SLURM run verifies it. No production run, biological threshold
or PAIRED downstream interpretation is implied by unit tests or dry-runs.

Dry-run examples (never use the global target for ATAC):

```bash
snakemake -n --cores 1 --default-resources tmpdir=/tmp atacseq_production
snakemake -n --cores 1 --default-resources tmpdir=/tmp \
  results/atacseq/provenance/SRR12774931/provenance.json \
  results/atacseq/provenance/SRR12774932/provenance.json
snakemake -n --cores 1 --default-resources tmpdir=/tmp atacseq_reproducibility
```
