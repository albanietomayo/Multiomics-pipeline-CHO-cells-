# ============================================================
# Reference genome resources
# ============================================================

rule reference_genome:
    input:
        script="workflow/scripts/fetch_reference_genome.py"

    output:
        fasta=config["reference"]["fasta"],
        gff3=config["reference"]["gff3"],
        gtf=config["reference"]["gtf"],
        sequence_report=config["reference"]["sequence_report"],
        metadata=config["reference"]["metadata"],
        sha256=config["reference"]["sha256"]

    conda:
        "../envs/reference.yaml"

    shell:
        """
        python {input.script}
        """
