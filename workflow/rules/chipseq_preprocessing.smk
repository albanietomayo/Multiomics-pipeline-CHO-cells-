configfile: "config/config.yaml"

# ------------------------------------------------------------
# ChIP-seq preprocessing consumes the dynamic processing cohort
# produced by workflow/rules/chipseq_metadata.smk.
#
# Metadata/planning is deliberately executed as a separate
# Snakemake stage by the SLURM entrypoint.  This keeps the
# generic FASTQ/QC/preprocessing engine independent of the
# ChIP-specific metadata-planning workflow.
# ------------------------------------------------------------

CHIPSEQ_PROCESSING_RUNS = (
    "results/metadata/chipseq/ena_snapshot_001/"
    "analysis_plan/chipseq_processing_runs.tsv"
)

CHIPSEQ_PROCESSING_SUMMARY = (
    "results/metadata/chipseq/ena_snapshot_001/"
    "analysis_plan/processing_run_summary.json"
)

# Legacy QC summary uses these paths relative to its own project copy.
for section, key, expected in (
    (
        "fastq",
        "validation_manifest",
        "results/metadata/validation_fastq_manifest.tsv",
    ),
    (
        "fastq",
        "validation_output_dir",
        "data/raw/fastq",
    ),
    (
        "preprocessing",
        "fastp_dir",
        "results/preprocessing/fastp",
    ),
    (
        "preprocessing",
        "reports_dir",
        "results/preprocessing/fastp/reports",
    ),
    (
        "qc",
        "raw_fastqc_dir",
        "results/qc/raw/fastqc",
    ),
    (
        "qc",
        "post_fastqc_dir",
        "results/qc/post/fastqc",
    ),
):
    if config[section][key] != expected:
        raise ValueError(
            f"QC summary path mismatch: {section}/{key}"
        )

config["fastq"]["validation_runs"] = CHIPSEQ_PROCESSING_RUNS

include: "common.smk"


rule chipseq_preprocessing_all:
    input:
        CHIPSEQ_PROCESSING_RUNS,
        CHIPSEQ_PROCESSING_SUMMARY,
        config["fastq"]["validation_download_report"],
        config["preprocessing"]["qc_by_fastq"],
        config["preprocessing"]["qc_by_job"],
        (
            f"{config['qc']['raw_multiqc_dir']}/"
            "multiqc_report.html"
        ),
        (
            f"{config['qc']['post_multiqc_dir']}/"
            "multiqc_report.html"
        )

    default_target:
        True
