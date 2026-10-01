# ChIP-seq v1.0.1 implementation validation

Validation date: 2026-09-30. This records an unstaged development patch, not a
release or a replacement scientific analysis. The initially completed checks
did not execute ChIP-seq processing. Subsequent clean-clone validation
re-executed **one of the 18 PASS analyses** from FASTQ: SRR20770297 IP versus
SRR20770287 Input, PRJNA865478, fixed 147 bp. Job 10413311 completed with exit
code 0:0. **The other 17 PASS analyses were not re-executed from FASTQ.**
Separate evidence is recorded under `provenance/chipseq/reproduction_validation/`.

For this single analysis, peak count, FRiP and peak/summit coordinates match
historical production. Strict reproduction is FAIL: the Input retains two
additional reads, peak numeric content differs, and both decoded signal hashes
differ. The cause remains unresolved. Treatment signal differs at six genomic
bases; none changes this analysis's 592 locus SPMR feature values or its
peak-coordinate-derived locus features. Substituting only this analysis's
reproduced products into the frozen downstream aggregation gives a
byte-identical final 37 x 524 ChIP context table and unchanged ChIP values
entering the 37 x 588 multi-omic master matrix. Ancillary condition/study Input
QC fields differ and are excluded from those final products. This is not
validation that all 18 analyses were reproduced from FASTQ.

## Repository safety

The starting branch was `multiomics-integration`, with HEAD
`826e45882b3196975a7a4a881326a2df1fc28b73`. Implementation is on
`chipseq-v1.0.1-reproducibility`, with the same HEAD and an empty staging index.
The v1.0.0 annotated tag object remains
`c140cf4777505b744403324bdede078dd3e368e7`; its peeled commit remains
`08e4bfd09e1954b72f62ae15e4c4c32cab162afd`.

All 345 previously tracked benchmark/integration files match their pre-patch
SHA256 values. The modified pre-existing tracked files are README.md, the
ChIP-specific raw-XML exclusion in .gitignore, and release-specific CITATION.cff
metadata.
The top-level Snakefile, existing downloader, RNA/ATAC code and
all remaining existing tracked files remain byte-identical. Twelve pre-existing
untracked status entries were recorded and left untouched. Nothing was staged,
committed, pushed, tagged or released.

CITATION.cff is now being prepared for v1.0.1 in both citation blocks. Its
release-specific changes are separate from the scientific validation; the
v1.0.0 tag and its original citation file remain immutable. Both release dates
must use the planned tag-creation day and be updated before tagging if that
day changes.
The planned tag date has not yet been supplied: both `date-released` fields
are intentionally pending rather than carrying the v1.0.0 date into v1.0.1.
Publication metadata is not ready until both dates are set to the same tag day,
no earlier than 2026-09-30. The concept DOI was verified directly against the
Zenodo record API; see `zenodo_concept_doi_verification.json`.

## Observed validation

| Check | Result |
| --- | --- |
| New Python syntax | 28 files parse and compile in memory |
| YAML parsing | Eight files parse |
| JSON parsing | 19 files parse |
| TSV structure | Eleven tables have unique columns and consistent row widths |
| Frozen scientific inputs | SHA256 checks pass |
| Frozen cohort/parameters | 18 analyses, 22 physical runs, four Inputs, eight narrow, ten broad; tests pass |
| Direct Table A9/combined environment pins | Match the original YAMLs and all four explicit package records |
| ChIP-seq tests | 113 passed; 40 subtests passed |
| Existing and new regression suite | 164 passed, 50 subtests passed; one pre-existing ATAC failure (see below) |
| Clean ChIP-seq dry-run | 49 jobs resolve successfully, including 22 downloads and 18 analyses |
| Saved clean-clone dry-run | `provenance/chipseq/reproduction_validation/clean_clone_dry_run.log`: 22 physical runs, 18 analyses, eight narrow, ten broad; no production jobs executed |
| Historical metadata replay | 232 sanitized JSON records replayed without XML; reconstructed plans match frozen plans |
| Security/portability scan | No credentials/email addresses in introduced material; no site paths/identifiers required by executable/configuration files |
| Existing frozen public outputs | All 345 SHA256 comparisons pass |
| Existing tracked files | README.md and .gitignore have ChIP-related changes; CITATION.cff has release-specific changes only |

The test harness used Python 3.14.7, Snakemake 9.25.2, pytest 9.1.1, PyYAML
6.0.3 and pysam 0.24.0 in a temporary environment. Conda 26.1.1 was used for
the dry-run. This harness is distinct from the frozen production environment
specification, which pins Python 3.10. The initial checks did not create
scientific Conda environments. The isolated
real reproduction uses the selected analysis's exact historical explicit package
record; see its dedicated evidence for completion and scientific comparisons.
The harness emitted the standard free-threaded Python extension/GIL warning.

Commands run include:

```bash
python -m pytest -q -p no:cacheprovider \
  workflow/tests/test_chipseq*.py \
  workflow/tests/test_prepare_chipseq_processing_runs.py
python -m pytest -q -p no:cacheprovider workflow/tests statistical_analysis/tests
snakemake --snakefile workflow/chipseq/Snakefile --use-conda --cores 1 -n chipseq_all
```

The dry-run used a fresh temporary reproduction tree containing source,
configuration and compact evidence, with no BAMs or generated results. It did
not download FASTQ, submit jobs, create Conda environments or fabricate
production provenance. The current metadata replay uses the bundled minimal scientific JSON snapshot,
without an external archive, raw XML or metadata network retrieval. Whole-record comparisons of both the
analysis and physical-run plans passed.

## Pre-existing regression failure

`workflow/tests/test_atac_production.py::ATACProductionContracts::test_real_production_snakemake_dag_direct_downstream_gate`
fails because the negative eligibility case expects a nonzero dry-run exit
code, while the checkpoint-based DAG returns zero before the checkpoint runs.
The identical test failure was independently reproduced against a temporary
source tree reconstructed from the untouched starting HEAD. ATAC code and its
tests have not been changed to suppress the failure.

Release readiness requires no new regressions relative to the starting HEAD,
not an unrelated ATAC fix. A complete untouched baseline run has 51 passing
tests, 10 passing subtests and this one failure. The patch run has 164 passing
tests, 50 passing subtests and the identical failure. All 113 ChIP-specific
tests pass; no new regressions were observed. The full suite is not described
as entirely passing.

## Portability and retained historical fields

The reviewed original helpers retain historical method/schema names, including
`slurm_job` and production reservation tokens. The adapter supplies the local
execution label `local_reproduction`; these are not scheduler requirements or
authentication tokens. The frozen scientific JSONs retain relative legacy
latest-submission strings to preserve their original hashes. The portable
entry point does not read those pointer files or invoke legacy preparation/
submission functions. Historical provenance retains original execution paths
as evidence, never as runtime defaults.

Large XML/source/run archives, raw FASTQ, BAMs, reference sequences/indexes,
signal tracks and generated production outputs are excluded from this patch.
Intended release-asset members and historical hashes are documented separately;
the assets themselves are not claimed to be published.
