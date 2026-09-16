import json
import re

with open("config/chipseq_filtering.json") as handle:
    CF = json.load(handle)
with open("config/chipseq_filtering_inputs.json") as handle:
    FILTER_INPUTS = json.load(handle)
RUNS = [item["run_accession"] for item in FILTER_INPUTS["runs"]]
if len(RUNS) != 2 or len(set(RUNS)) != 2 or any(not re.fullmatch(r"[DES]RR[0-9]+", run) for run in RUNS):
    raise ValueError("Expected two distinct pilot runs")
OUT = "results/chipseq/filtering"
CODE = "workflow/scripts/chipseq_filtering_support.py"
PLAN = "config/chipseq_filtering_inputs.json"
CHIP_ELIGIBILITY_SAMPLES = "config/samples.tsv"
CHIP_ELIGIBILITY_CONDITIONS = "snapshots/chipseq/control_validation_001/chipseq_conditions.tsv"
CHIP_ELIGIBILITY_PILOT = "config/chipseq_pilot.json"
CHIP_ELIGIBILITY_REPORT = OUT + "/pilot_eligibility.json"

include: "chipseq_eligibility.smk"

rule chipseq_filtering_all:
    input:
        expand(OUT + "/{run}/filtered.bam", run=RUNS),
        expand(OUT + "/{run}/filtered.bam.csi", run=RUNS),
        expand(OUT + "/{run}/filtered_validation.json", run=RUNS),
        expand(OUT + "/{run}/flagstat.txt", run=RUNS),
        expand(OUT + "/{run}/idxstats.tsv", run=RUNS),
        expand(OUT + "/{run}/samtools_stats.txt", run=RUNS),
        expand(OUT + "/{run}/picard_metrics.txt", run=RUNS),
        expand(OUT + "/{run}/filtering_flow.tsv", run=RUNS),
        OUT + "/filtering_qc.tsv"
    default_target: True

rule chipseq_filter_stage:
    input:
        upstream=lambda w: FILTER_INPUTS["source_root"] + "/" + next(
            r["source_relative"] for r in FILTER_INPUTS["runs"] if r["run_accession"] == w.run),
        plan=PLAN, configuration="config/chipseq_filtering.json", code=CODE,
        reference="inputs/reference/genome_plus_mt.fa.fai", eligible=CHIP_ELIGIBILITY_REPORT
    output:
        bam=temp("inputs/{run}/raw.sorted.bam"),
        verified=OUT + "/{run}/input_validation.json"
    wildcard_constraints: run="|".join(RUNS)
    conda: "../envs/chipseq_filtering.yaml"
    log: OUT + "/{run}/stage.log"
    shell:
        "python {input.code:q} stage --run {wildcards.run:q} --output {output.bam:q} "
        "--report {output.verified:q} > {log:q} 2>&1"

rule chipseq_mark_duplicates:
    input:
        bam="inputs/{run}/raw.sorted.bam", verified=OUT + "/{run}/input_validation.json",
        configuration="config/chipseq_filtering.json"
    output:
        bam=temp("work/chipseq_filtering/{run}/dupmarked.sorted.bam"),
        metrics=OUT + "/{run}/picard_metrics.txt"
    params:
        temporary="work/chipseq_filtering/{run}/picard_tmp", heap=CF["picard_heap_mb"]
    threads: CF["threads"]
    resources: mem_mb=18000
    conda: "../envs/chipseq_duplicates.yaml"
    log: OUT + "/{run}/picard.log"
    shell:
        r"""
        mkdir -p {params.temporary:q}
        picard -Xmx{params.heap}m -XX:ActiveProcessorCount={threads} MarkDuplicates \
            INPUT={input.bam:q} OUTPUT={output.bam:q} METRICS_FILE={output.metrics:q} \
            TMP_DIR={params.temporary:q} REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=false \
            READ_NAME_REGEX=null DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES \
            CREATE_INDEX=false VALIDATION_STRINGENCY=STRICT \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=256 > {log:q} 2>&1
        """

rule chipseq_filter_bam:
    input:
        bam="work/chipseq_filtering/{run}/dupmarked.sorted.bam",
        metrics=OUT + "/{run}/picard_metrics.txt", code=CODE,
        configuration="config/chipseq_filtering.json", plan=PLAN
    output:
        bam=OUT + "/{run}/filtered.bam", report=OUT + "/{run}/filtering_qc.json",
        flow=OUT + "/{run}/filtering_flow.tsv"
    resources: mem_mb=2000
    conda: "../envs/chipseq_filtering.yaml"
    log: OUT + "/{run}/filtering.log"
    shell:
        "python {input.code:q} filter --run {wildcards.run:q} --bam {input.bam:q} "
        "--metrics {input.metrics:q} --output {output.bam:q} --report {output.report:q} "
        "--flow {output.flow:q} > {log:q} 2>&1"

rule chipseq_filtered_index:
    input: OUT + "/{run}/filtered.bam"
    output: OUT + "/{run}/filtered.bam.csi"
    conda: "../envs/chipseq_filtering.yaml"
    shell: "samtools index -c {input:q} {output:q}"

rule chipseq_filtered_qc:
    input:
        bam=OUT + "/{run}/filtered.bam", index=OUT + "/{run}/filtered.bam.csi",
        report=OUT + "/{run}/filtering_qc.json", code=CODE, configuration="config/chipseq_filtering.json"
    output:
        verified=OUT + "/{run}/filtered_validation.json", flagstat=OUT + "/{run}/flagstat.txt",
        idxstats=OUT + "/{run}/idxstats.tsv", stats=OUT + "/{run}/samtools_stats.txt"
    conda: "../envs/chipseq_filtering.yaml"
    log: OUT + "/{run}/filtered_qc.log"
    shell:
        "samtools quickcheck -v {input.bam:q} 2> {log:q}; "
        "samtools flagstat {input.bam:q} > {output.flagstat:q} 2>> {log:q}; "
        "samtools idxstats {input.bam:q} > {output.idxstats:q} 2>> {log:q}; "
        "samtools stats {input.bam:q} > {output.stats:q} 2>> {log:q}; "
        "python {input.code:q} verify --run {wildcards.run:q} --bam {input.bam:q} "
        "--report {input.report:q} --output {output.verified:q} >> {log:q} 2>&1"

rule chipseq_filtering_summary:
    input:
        reports=expand(OUT + "/{run}/filtering_qc.json", run=RUNS),
        verified=expand(OUT + "/{run}/filtered_validation.json", run=RUNS), code=CODE
    output: OUT + "/filtering_qc.tsv"
    conda: "../envs/chipseq_filtering.yaml"
    shell: "python {input.code:q} summarize --output {output:q} {input.reports:q}"
