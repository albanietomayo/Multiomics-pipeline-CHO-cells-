# Thesis errata: ChIP-seq code availability

The original v1.0.0 public release inadvertently omitted the upstream ChIP-seq
processing implementation used during the project, while retaining the
downstream ChIP-seq integration code, frozen outputs and provenance.

Version v1.0.1 adds the missing processing implementation and associated
software environments.

## Statements inaccurate relative to the public v1.0.0 release

**Table A9.** The table lists five stage-specific ChIP-seq environment files
that were absent from the public v1.0.0 repository:

- `workflow/envs/chipseq_alignment.yaml`
- `workflow/envs/chipseq_duplicates.yaml`
- `workflow/envs/chipseq_filtering.yaml`
- `workflow/envs/chipseq_phantompeakqualtools.yaml`
- `workflow/envs/chipseq_peak_calling.yaml`

v1.0.1 restores these exact names with the project's validated software
versions. Final validated lean production historically used the combined
`workflow/envs/chipseq_production.yaml` environment; the restored stage files
do not change that historical fact.

**Section 2.5 / workflow architecture.** The description of a ChIP-seq module
within the reproducible Snakemake architecture overstated what was actually
available in the public v1.0.0 repository. That release contained downstream
ChIP-seq integration material and frozen products/provenance, but omitted the
complete upstream public processing implementation. v1.0.1 restores it through
the isolated portable ChIP-seq workflow.

**Conclusion 1.** The description of infrastructure going "from raw data to
harmonized products" was too broad when interpreted as a claim about public
v1.0.0 ChIP-seq code availability: the raw-to-product processing implementation
used in production was not included. This is a code-availability discrepancy,
not a correction to frozen ChIP-seq scientific results.

**Existing correct disclosure.** The v1.0.0 README already stated that upstream
ChIP alignment/peak processing lived in a separate project and that the public
repository lacked a complete local ChIP production DAG. This disclosure was
consistent with the release and does not require correction. The discrepancies
above do not imply that every thesis statement about ChIP-seq was inaccurate.

## Validation scope

As validation of the public reproduction interface, **one of the 18 PASS
ChIP-seq analyses** was re-executed from FASTQ in a clean clone: SRR20770297 IP
versus SRR20770287 Input, PRJNA865478, fixed 147 bp fragment size. Peak count,
FRiP and peak/summit coordinates matched historical production. Small localized
numerical signal differences and a two-read Input discrepancy remained, with
their cause unresolved.

Substituting the reproduced products of this one analysis into the frozen
downstream aggregation preserved its locus SPMR and peak-coordinate features
and produced a byte-identical final 37 x 524 ChIP context table. The ChIP values
entering the 37 x 588 multi-omic master matrix were unchanged. Ancillary Input
QC fields differed in condition/study tables and are excluded from those final
products. **The remaining 17 PASS analyses were not re-executed from FASTQ.**

No scientific results, benchmark labels or frozen multi-omic products were
modified. Strict reproduction identity was not obtained for the single
validation analysis; no full ChIP-seq dataset reproduction is claimed.
