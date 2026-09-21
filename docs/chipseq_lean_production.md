# Lean ChIP-seq production mode

This mode executes exactly one eligible analysis at a time and inherits its
IP/Input pair from
`snapshots/chipseq/incremental_planning_validation_001/chipseq_analysis_plan.tsv`.
`config/chipseq_peak_calling_plan.tsv` supplies the reviewed peak mode and
fragment policy. Any disagreement between these layers fails before execution.

The worker uses the existing download validator, fastp settings, Bowtie2
settings, alignment QC, Picard duplicate marking, nuclear filtering,
PhantomPeakQualTools parsing, MACS3 contract, FRiP, peak QC, and provenance
implementations. It reads the shared Bowtie2 index in place. A global `flock`
under the production root serializes jobs.

For each physical run, node-local storage follows this validated deletion
order:

1. Download and ENA size/MD5-validate raw FASTQ in `$TMPDIR`.
2. Run fastp and FastQC; validate processed gzip and positive fastp read count,
   then remove raw FASTQ.
3. Create, quickcheck, index and fully QC the sorted alignment BAM; then remove
   processed FASTQ.
4. Create and quickcheck the Picard duplicate-marked BAM; then remove the raw
   alignment BAM/CSI.
5. Create, index and fully verify the filtered BAM; then remove the
   duplicate-marked BAM.

A control is processed once and atomically persisted under
`results/chipseq/production/shared_controls/<accession>/` as only its filtered
BAM/CSI and compact QC/provenance. An IP filtered BAM stays in `$TMPDIR` and is
never copied to persistent storage. MACS3 bedGraphs are validated by the
existing QC and then stored as lossless deterministic gzip files. The final
artifact directory and completion marker become visible atomically.

The worker records scratch used/free bytes at job start, environment creation,
after every download/alignment/duplicate/filter stage, after peak calling and
after persistence. A ten-second sampler supplies maximum observed usage.

The check-only first-pair command is:

```bash
python workflow/slurm/submit_chipseq_production.py \
  --analysis-id SRR20770297 \
  --control-run SRR20770287 \
  --shared-reference-root "$HOME/TFM_multiomics_pipeline_benchmark/resources/reference/CriGri-PICRH-1.0/chipseq" \
  --development-dirty-check
```

`--submit` is the only path to `sbatch`, and dirty submissions are forbidden.
After all analyses for a control are complete, verify cleanup eligibility with:

```bash
python workflow/scripts/chipseq_production.py control-cleanup \
  --production-root results/chipseq/production \
  --control SRR20770287
```

This check verifies every expected analysis manifest and artifact. It prints
`SAFE_TO_REMOVE_CONTROL_BAM` but changes nothing. An operator may explicitly
add `--remove` to remove only the validated control BAM/CSI; compact control QC
and the cleanup marker remain. The production worker never removes a control
automatically, including after failures.
