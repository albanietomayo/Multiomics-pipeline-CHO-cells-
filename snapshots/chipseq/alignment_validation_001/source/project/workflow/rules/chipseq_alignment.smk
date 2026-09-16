import json
import re

configfile: "config/config.yaml"

with open("config/chipseq_alignment.json") as handle:
    CA = json.load(handle)
with open("config/chipseq_alignment_inputs.json") as handle:
    ALIGN_INPUTS = json.load(handle)
RUNS = [item["run_accession"] for item in ALIGN_INPUTS["runs"]]
if len(RUNS) != 2 or len(set(RUNS)) != 2 or any(not re.fullmatch(r"[DES]RR[0-9]+", r) for r in RUNS):
    raise ValueError("Expected exactly two distinct pilot runs")
if CA["schema_version"] != 1 or CA["bowtie2_mode"] != "local" or CA["bowtie2_preset"] != "very-sensitive-local":
    raise ValueError("Unsupported alignment configuration; review before changing mode")
REF = config["reference"]
RDIR = REF["dir"] + "/chipseq"
FASTA = RDIR + "/genome_plus_mt.fa"
PREFIX = RDIR + "/bowtie2_index/genome_plus_mt"
INDEX = [PREFIX + suffix for suffix in (".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")]
OUT = "results/chipseq/alignment"
CODE = "workflow/scripts/chipseq_alignment_support.py"
CHIP_ELIGIBILITY_SAMPLES = "config/samples.tsv"
CHIP_ELIGIBILITY_CONDITIONS = "snapshots/chipseq/control_validation_001/chipseq_conditions.tsv"
CHIP_ELIGIBILITY_PILOT = "config/chipseq_pilot.json"
CHIP_ELIGIBILITY_REPORT = OUT + "/pilot_eligibility.json"

include: "chipseq_eligibility.smk"

rule chipseq_alignment_all:
    input:
        expand(OUT + "/{run}/raw.sorted.bam", run=RUNS),
        expand(OUT + "/{run}/raw.sorted.bam.csi", run=RUNS),
        expand(OUT + "/{run}/flagstat.txt", run=RUNS),
        expand(OUT + "/{run}/idxstats.tsv", run=RUNS),
        expand(OUT + "/{run}/samtools_stats.txt", run=RUNS),
        expand(OUT + "/{run}/alignment_qc.json", run=RUNS),
        OUT + "/alignment_qc.tsv",
        FASTA + ".fai",
        RDIR + "/reference_provenance.json"
    default_target: True

rule chipseq_nuclear_reference:
    input:
        script="workflow/scripts/fetch_reference_genome.py",
        configuration="config/config.yaml",
        eligible=CHIP_ELIGIBILITY_REPORT
    output:
        fasta=REF["fasta"], gff3=REF["gff3"], gtf=REF["gtf"],
        sequence_report=REF["sequence_report"], metadata=REF["metadata"], sha256=REF["sha256"]
    conda: "../envs/reference.yaml"
    log: OUT + "/logs/nuclear_reference.log"
    shell: "python {input.script:q} > {log:q} 2>&1"

rule chipseq_mapping_reference:
    input:
        nuclear=REF["fasta"], hashes=REF["sha256"], metadata=REF["metadata"],
        code=CODE, configuration="config/chipseq_alignment.json"
    output:
        fasta=FASTA, mt=RDIR + "/mitochondrial.fa", provenance=RDIR + "/reference_provenance.json"
    conda: "../envs/chipseq_alignment.yaml"
    log: OUT + "/logs/mapping_reference.log"
    shell: "python {input.code:q} reference > {log:q} 2>&1"

rule chipseq_bowtie2_index:
    input: fasta=FASTA, provenance=RDIR + "/reference_provenance.json"
    output: INDEX
    params: prefix=PREFIX, directory=RDIR + "/bowtie2_index"
    threads: CA["threads"]
    conda: "../envs/chipseq_alignment.yaml"
    log: OUT + "/logs/bowtie2_build.log"
    shell:
        "mkdir -p {params.directory:q}; "
        "bowtie2-build --large-index --threads {threads} {input.fasta:q} {params.prefix:q} > {log:q} 2>&1"

rule chipseq_fasta_index:
    input: FASTA
    output: FASTA + ".fai"
    conda: "../envs/chipseq_alignment.yaml"
    shell: "samtools faidx {input:q}"

rule chipseq_align_pilot_run:
    input:
        fastq="inputs/{run}.fastq.gz", index=INDEX, eligible=CHIP_ELIGIBILITY_REPORT,
        verified="inputs/verified.json", configuration="config/chipseq_alignment.json"
    output: bam=OUT + "/{run}/raw.sorted.bam"
    wildcard_constraints: run="|".join(RUNS)
    params:
        prefix=PREFIX, seed=CA["bowtie2_seed"], directory=OUT + "/{run}",
        bt_threads=lambda wildcards, threads: max(1, threads - 2)
    threads: CA["threads"]
    conda: "../envs/chipseq_alignment.yaml"
    log: bt=OUT + "/{run}/bowtie2.log", sort=OUT + "/{run}/sort.log"
    shell:
        r"""
        mkdir -p {params.directory:q}
        bowtie2 --local --very-sensitive-local --seed {params.seed} \
            -p {params.bt_threads} -x {params.prefix:q} -U {input.fastq:q} \
            --rg-id {wildcards.run:q} --rg SM:{wildcards.run} --rg LB:{wildcards.run} --rg PL:ILLUMINA \
            2> {log.bt:q} | \
            samtools sort -@ 1 -m 768M -T {params.directory:q}/sort_tmp \
                -o {output.bam:q} - 2> {log.sort:q}
        samtools quickcheck -v {output.bam:q}
        """

rule chipseq_bam_index:
    input: OUT + "/{run}/raw.sorted.bam"
    output: OUT + "/{run}/raw.sorted.bam.csi"
    conda: "../envs/chipseq_alignment.yaml"
    shell: "samtools index -c {input:q} {output:q}"

rule chipseq_alignment_metrics:
    input:
        bam=OUT + "/{run}/raw.sorted.bam", index=OUT + "/{run}/raw.sorted.bam.csi",
        code=CODE, plan="config/chipseq_alignment_inputs.json", reference=RDIR + "/reference_provenance.json",
        configuration="config/chipseq_alignment.json"
    output:
        flagstat=OUT + "/{run}/flagstat.txt", idxstats=OUT + "/{run}/idxstats.tsv",
        stats=OUT + "/{run}/samtools_stats.txt", qc=OUT + "/{run}/alignment_qc.json"
    conda: "../envs/chipseq_alignment.yaml"
    log: OUT + "/{run}/qc.log"
    shell:
        "samtools flagstat {input.bam:q} > {output.flagstat:q} 2> {log:q}; "
        "samtools idxstats {input.bam:q} > {output.idxstats:q} 2>> {log:q}; "
        "samtools stats {input.bam:q} > {output.stats:q} 2>> {log:q}; "
        "python {input.code:q} qc --run {wildcards.run:q} --bam {input.bam:q} "
        "--output {output.qc:q} 2>> {log:q}"

rule chipseq_alignment_summary:
    input: reports=expand(OUT + "/{run}/alignment_qc.json", run=RUNS), code=CODE
    output: OUT + "/alignment_qc.tsv"
    conda: "../envs/chipseq_alignment.yaml"
    shell: "python {input.code:q} summarize --output {output:q} {input.reports:q}"
