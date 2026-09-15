# ChIP-seq automation and experimental eligibility

Two explicit Snakemake entry points are maintained:

- `workflow/rules/chipseq_metadata.smk`: catalog plan, XML retrieval, labels,
  conditions, candidate controls and validation of the curated protocol review.
- `workflow/rules/chipseq_preprocessing.smk`: eligibility of the selected pilot,
  pilot manifest, download, preprocessing and QC via `common.smk`.

Both include `workflow/rules/chipseq_eligibility.smk`. That rule calls
`workflow/scripts/validate_chipseq_eligibility.py`. The preprocessing selection
depends on its report, and the SLURM worker explicitly runs the eligibility
target before preparing the manifest or downloading FASTQ.

The scientific review is `config/chipseq_experimental_eligibility.json`.
It is curated input, not an automatically inferred conclusion. The validator
checks its provenance, full catalog coverage, role/study consistency and pilot
approval. New or changed catalog/condition evidence requires a new review.
`excluded`, `review_required`, missing and unknown decisions cannot enter the
pilot. A successful eligibility check does not approve a control association,
establish independent biological replicates, or certify read quality.

## Preparing a submission

From the `chipseq-metadata-audit` checkout:

```bash
bash workflow/slurm/submit_chipseq_preprocessing.sh --check
```

This creates an isolated source copy, checks eligibility and pilot selection,
and records hashes and source provenance without submitting a SLURM job.
It does not validate the Snakemake DAG or run QC.

For an intentional new processing attempt:

```bash
bash workflow/slurm/submit_chipseq_preprocessing.sh --submit
```

Every `--submit` creates a new attempt; it does not resume an earlier attempt
or prevent duplicate concurrent submissions. Do not resubmit simply to install
or check this integration. Existing submissions retain their original sources.

The worker prepares Conda on the compute node, executes the eligibility target,
generates/checks the manifest, performs a dry-run, then executes preprocessing
and verifies saved output copies using the existing worker's workflow.

## Snakemake verification on VERA

Use the project's existing Snakemake environment. From the project root:

```bash
snakemake --snakefile workflow/rules/chipseq_metadata.smk chipseq_metadata_all --cores 1 -n
snakemake --snakefile workflow/rules/chipseq_preprocessing.smk chipseq_preprocessing_all --cores 1 -n
```

For an actual local execution of the eligibility check only:

```bash
snakemake --snakefile workflow/rules/chipseq_preprocessing.smk chipseq_experimental_eligibility --cores 1 -p
```

No FASTQ target is requested by that last command. Checkpoint-dependent parts
of the full dry-run may remain deferred until their metadata are available.

## Source control and scope

Review and commit the JSON, the three new workflow files, this document, and
the changes to both entry points and `chipseq_preprocessing.sbatch` together
with the earlier uncommitted pilot files. The installer does not stage, commit,
push, submit, or modify archived evidence.

This integration covers the audited metadata and preprocessing pilot. It is not
a full ChIP-seq alignment/peak-calling workflow or an automatic scientific
review of new studies. Keep control decisions separate and preserve ambiguous
inputs without merging or relabelling them by inference.
