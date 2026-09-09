# ============================================================
# CHO multi-omics pipeline
# Main Snakemake entry point
# ============================================================

configfile: "config/config.yaml"


include: "workflow/rules/metadata.smk"
include: "workflow/rules/reference.smk"
include: "workflow/rules/common.smk"
include: "workflow/rules/atacseq.smk"
include: "workflow/rules/rnaseq.smk"


rule all:
    input:
        "config/samples.tsv",
        "results/metadata/ena_metadata_after_curation.tsv",
        "results/metadata/reports/curation_comparison.tsv",
        "results/metadata/reports/unresolved_after_curation.tsv",
        config["fastq"]["manifest"],
        config["fastq"]["validation_manifest"],
        get_validation_fastq_outputs,
        config["fastq"]["validation_download_report"],
        get_raw_fastqc_outputs,
        f"{config['qc']['raw_multiqc_dir']}/multiqc_report.html",
        get_fastp_outputs,
        get_post_fastqc_outputs,
        f"{config['qc']['post_multiqc_dir']}/multiqc_report.html",
        config["preprocessing"]["qc_by_fastq"],
        config["preprocessing"]["qc_by_job"],

        # Reference genome
        config["reference"]["fasta"],
        config["reference"]["gff3"],
        config["reference"]["gtf"],
        config["reference"]["sequence_report"],
        config["reference"]["metadata"],
        config["reference"]["sha256"],

        # RNA-seq validation outputs
        config["rnaseq"]["counts_matrix"],
        config["rnaseq"]["alignment_qc_summary"],
        lambda wildcards: get_rnaseq_alignment_qc_outputs(),
        config["rnaseq"]["strandedness_summary"],
        config["rnaseq"]["strandedness_classification"],
        config["rnaseq"]["featurecounts_qc_summary"],

    default_target: True
