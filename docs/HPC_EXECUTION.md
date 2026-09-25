# VERA/C3SE execution examples

`run_rnaseq_smoke.sh` and the scripts in `workflow/slurm/` are preserved as
site-specific execution and diagnostic examples. Their `#SBATCH` account and
partition, module names, scratch layout, project paths, output paths and job
resources reflect VERA/C3SE runs. They are not portable entry points and should
not be submitted unchanged on another cluster. Adapt those settings, input
locations and the Conda/module setup to the target site, then dry-run the
relevant Snakemake target before submitting work. These scripts are historical
examples; editing them in place could obscure the recorded execution context.

See [the main README](../README.md#path-b-workflow-reproduction) for the
current workflow boundary. The top-level `Snakefile` includes metadata,
reference, RNA-seq and ATAC-seq rules; its default `all` target does not
include ATAC production or final 37-locus multiomic integration. ChIP-seq
production was performed in a separate project. Large FASTQs, alignments,
tracks and external production artifacts require substantial cluster storage.
