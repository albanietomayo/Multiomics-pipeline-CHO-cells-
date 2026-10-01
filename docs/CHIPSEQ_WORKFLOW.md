# ChIP-seq reproduction workflow

The proposed v1.0.1 patch publishes a portable adaptation of the validated
ChIP-seq processing. v1.0.0 remains the frozen scientific release. Existing
benchmark, RNA, ATAC, ChIP integration, master, Core and Extended products are
unchanged. No full ChIP-seq production rerun is claimed for this patch.

Clean-clone validation re-executed **one of the 18 PASS analyses** from FASTQ:
SRR20770297 IP versus SRR20770287 Input (PRJNA865478, fixed 147 bp). For this
analysis, peak count, FRiP and peak/summit coordinates matched history, while
strict signal identity failed and a two-read Input discrepancy remains
unresolved. Substitution of this single analysis's reproduced products into
the frozen downstream aggregation preserved its locus SPMR and coordinate
features, yielded a byte-identical final 37 x 524 ChIP context table, and left
the ChIP values supplied to the 37 x 588 master matrix unchanged. Ancillary
condition/study Input QC fields differ and are excluded from those final
products. **The other 17 PASS analyses were not re-executed from FASTQ.**

## Entry point and targets

Run commands from the repository root with the orchestration dependencies in
`environment.yml` installed. Snakemake requires a supported Conda installation.
The top-level Snakefile and its RNA/ATAC behavior are unchanged.

```bash
snakemake --snakefile workflow/chipseq/Snakefile --use-conda --cores 1 -n chipseq_all
```

This dry-run resolves 22 physical runs, four shared Inputs and 18 analyses:
eight narrow and ten broad. It does not download reads, create rule environments,
submit jobs, require existing BAMs, or generate execution provenance.

| Target | Purpose |
| --- | --- |
| `chipseq_metadata_all` | Retrieve/replay ENA metadata, rebuild eligibility and pairing, and compare with frozen plans |
| `chipseq_plan_validate` | Verify immutable inputs, authoritative pairings and historical parameter contracts |
| `chipseq_reference_all` | Acquire the public nuclear/mitochondrial reference and build its Bowtie2 index |
| `chipseq_download_all` | Download and size/MD5-verify the 22 physical FASTQ files |
| `chipseq_all` | Process shared controls and IPs; generate peaks, SPMR, FRiP, QC and manifests |

The production target consumes the committed authoritative cohort. Metadata
reconstruction is a separate target; current archive metadata cannot silently
replace biological pairings. Differences fail the comparison. The supported
cohort is single-end Illumina; unsupported or ambiguous inputs fail closed.

Execution locations belong in `config/chipseq.yaml`. Scientific settings remain
in the frozen JSON/TSV files. Default outputs are under
`results/chipseq/reproduction/`, separate from historical scientific products.
Snakemake execution-location overrides are forwarded to the adapter. Keep
`config/chipseq_metadata.yaml` consistent when relocating metadata outputs.

Full execution requires substantial compute, network access, storage and a
supported Linux toolchain. It is separate from validation of this patch:

```bash
snakemake chipseq_all --snakefile workflow/chipseq/Snakefile --use-conda --cores 16 \
  --resources chipseq_slots=2
```

The portable adapter limits execution to at most the configured one or two
workers and locks shared controls independently of the scheduler. Intermediate
IP files live in a private scratch directory. Configure `scratch_root`, or use
TMPDIR/the operating system temporary directory. Only job-owned intermediates
are removed. Shared control BAMs and downloaded raw reads remain available;
there is no automatic control cleanup. No SLURM account, partition, module,
branch name, personal directory or queue submission is required.

## Validated scientific behavior

The cohort contains 18 IP runs and four Inputs, with H3K4me3, H3K27ac, H3K9me3,
H3K4me1, H3K27me3 and H3K36me3. Pairings are inherited exactly from
`snapshots/chipseq/incremental_planning_validation_001/chipseq_analysis_plan.tsv`.

Processing follows the validated lean worker: ENA validation; fastp `-Q -L -G`
and FastQC; Bowtie2 local/very-sensitive-local/seed=0 with run-specific RG/SM/LB
and PL:ILLUMINA; coordinate sorting/CSI indexing/full alignment QC; Picard
MarkDuplicates with optical detection disabled; and the original pysam nuclear
filter and full-read verification. Filtering requires MAPQ >=30, rejects MAPQ
255, and excludes flags 3844.

The reference is CriGri-PICRH-1.0 / GCF_003668045.3 plus mitochondrial
NC_007936.1. The mitochondrial contig participates in mapping and is excluded
from filtered BAMs and the 2,366,634,374 bp nuclear span. Reference sequence
hashes and FAI identity must match frozen evidence. Newly built gzip containers,
indexes and inventories record their actual hashes; historical packaging hashes
are retained separately.

PRJNA865478 uses fixed 147 bp fragments. PRJEB9291 requires independent
PhantomPeakQualTools estimates, with no fixed fallback. Newly executed estimates
are recorded and can be compared with historical values; they are not replaced
with historical numbers. MACS3 uses q=0.01, broad cutoff=0.1, keep-dup all,
scale-to small, `--nomodel`, and `-B --SPMR`. Signal is lossless deterministic
gzip bedGraph, not BigWig. FRiP is IP reads overlapping IP peaks divided by IP
filtered reads; Input overlap and IP/Input ratio remain descriptive metrics.
No blacklist filtering or ChIP replicate concordance/IDR was performed.

## Environments and historical evidence

Table A9 describes five modular environments: alignment, duplicates, filtering,
PhantomPeakQualTools and peak calling. Their original YAML specifications are
published without changing package versions. Final validated lean production
used the combined `workflow/envs/chipseq_production.yaml`; the portable jobs use
that same specification. The five modular YAMLs are not the single environment
used by all final production jobs.

The provenance directory records the original per-analysis source commit,
source hashes and four historical explicit Conda package sets. The portable
implementation is an adaptation for reproducibility and is not claimed to be
the byte-identical source file set for every historical run. Historical package
sets differ in some transitive packages despite identical Table A9 tool versions.
A future solve is not claimed to recreate all historical transitive packages.

## Sanitized metadata replay

The committed `snapshots/chipseq/sanitized_metadata_001/` snapshot contains only
ENA experiment/sample accessions, source query URLs, retrieval dates, historical
source hashes and the five scientific label fields consumed by the reviewed
planning functions. It contains no submitter/contact fields or email addresses.
The isolated metadata target replays JSON without fetching or writing raw XML.
The original inventory's XML paths and hashes describe historical evidence,
not XML bodies distributed with this repository. Reconstructed full-record
analysis/pairing and physical-run plans must match the frozen plans exactly.

Raw ENA XML and `ena_xml_snapshot.tar.gz` must never be included in Git, GitHub
release assets or Zenodo. Only accessions, query/date records, required scientific
metadata, inventory, validation summary, frozen plans and the historical archive
checksum are public. The historical archive SHA256 remains:

```text
5b01f03bacf45548c004c70c28e46222e6313a0f457a3549f2165f174c0c54cf
```

See [validation](CHIPSEQ_V101_VALIDATION.md),
[code-availability note](../ERRATA_THESIS.md) and
[production provenance](../provenance/chipseq/production_2026-09-23/README.md).
