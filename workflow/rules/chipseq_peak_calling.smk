import json
import re

with open("config/chipseq_peak_calling.json") as handle:
    PC = json.load(handle)

with open("config/chipseq_peak_calling_inputs.json") as handle:
    PEAK_INPUTS = json.load(handle)

OUT = "results/chipseq/peak_calling"
CODE = "workflow/scripts/chipseq_peak_calling_support.py"
PLAN = "config/chipseq_peak_calling_inputs.json"

IP = PC["treatment_run"]
CONTROL = PC["control_run"]
RUNS = [CONTROL, IP]
NAME = PC["analysis_name"]

if len(set(RUNS)) != 2:
    raise ValueError("Peak calling requires distinct IP and Input runs")

if any(not re.fullmatch(r"[DES]RR[0-9]+", run) for run in RUNS):
    raise ValueError("Unexpected run accession")

if PC["format"] != "BAM":
    raise ValueError("Pilot peak calling requires SINGLE-END BAM mode")

if PC["peak_mode"] != "narrow":
    raise ValueError("Pilot peak calling currently supports narrow peaks only")

if PC["keep_dup"] != "all":
    raise ValueError("Upstream-deduplicated BAMs require --keep-dup all")

if PC["spmr"] and not PC["store_bdg"]:
    raise ValueError("--SPMR requires bedGraph output")

CHIP_ELIGIBILITY_SAMPLES = "config/samples.tsv"
CHIP_ELIGIBILITY_CONDITIONS = "snapshots/chipseq/control_validation_001/chipseq_conditions.tsv"
CHIP_ELIGIBILITY_PILOT = "config/chipseq_pilot.json"
CHIP_ELIGIBILITY_REPORT = OUT + "/pilot_eligibility.json"

include: "chipseq_eligibility.smk"

EXTRA_FLAGS = []

if PC["store_bdg"]:
    EXTRA_FLAGS.append("-B")

if PC["spmr"]:
    EXTRA_FLAGS.append("--SPMR")

if PC["call_summits"]:
    EXTRA_FLAGS.append("--call-summits")

if PC["cutoff_analysis"]:
    EXTRA_FLAGS.append("--cutoff-analysis")

EXTRA_FLAGS = " ".join(EXTRA_FLAGS)

if PC["model"] == "auto":
    MODEL_FLAGS = ""
elif PC["model"] == "fixed":
    if int(PC.get("extsize", 0)) <= 0:
        raise ValueError("Fixed MACS3 model requires positive extsize")
    MODEL_FLAGS = "--nomodel --extsize {}".format(int(PC["extsize"]))
else:
    raise ValueError("Unsupported MACS3 model mode")


rule chipseq_peak_calling_all:
    input:
        OUT + "/peaks/" + NAME + "_peaks.narrowPeak",
        OUT + "/peaks/" + NAME + "_peaks.xls",
        OUT + "/peaks/" + NAME + "_summits.bed",
        OUT + "/peaks/" + NAME + "_treat_pileup.bdg",
        OUT + "/peaks/" + NAME + "_control_lambda.bdg",
        OUT + "/peak_qc.json",
        OUT + "/peak_qc.tsv",
        OUT + "/peak_calling_parameters.json",
        OUT + "/input_provenance.json"
    default_target: True


rule chipseq_peak_stage:
    input:
        plan=PLAN,
        configuration="config/chipseq_peak_calling.json",
        code=CODE,
        eligible=CHIP_ELIGIBILITY_REPORT
    output:
        bam=temp("inputs/{run}/filtered.bam"),
        index=temp("inputs/{run}/filtered.bam.csi"),
        report=OUT + "/inputs/{run}_validation.json"
    wildcard_constraints:
        run="|".join(RUNS)
    conda:
        "../envs/chipseq_peak_calling.yaml"
    shell:
        """
        python {input.code:q} stage \
            --run {wildcards.run:q} \
            --plan {input.plan:q} \
            --bam {output.bam:q} \
            --index {output.index:q} \
            --report {output.report:q}
        """


rule chipseq_call_peaks:
    input:
        ip="inputs/" + IP + "/filtered.bam",
        ip_index="inputs/" + IP + "/filtered.bam.csi",
        control="inputs/" + CONTROL + "/filtered.bam",
        control_index="inputs/" + CONTROL + "/filtered.bam.csi",
        ip_validation=OUT + "/inputs/" + IP + "_validation.json",
        control_validation=OUT + "/inputs/" + CONTROL + "_validation.json",
        reference="inputs/reference/genome_plus_mt.fa.fai",
        configuration="config/chipseq_peak_calling.json"
    output:
        peaks=OUT + "/peaks/" + NAME + "_peaks.narrowPeak",
        xls=OUT + "/peaks/" + NAME + "_peaks.xls",
        summits=OUT + "/peaks/" + NAME + "_summits.bed",
        treat_bdg=OUT + "/peaks/" + NAME + "_treat_pileup.bdg",
        control_bdg=OUT + "/peaks/" + NAME + "_control_lambda.bdg",
        version=OUT + "/macs3_version.txt"
    params:
        outdir=OUT + "/peaks",
        name=NAME,
        gsize=PC["effective_genome_size"],
        qvalue=PC["qvalue"],
        keep_dup=PC["keep_dup"],
        scale_to=PC["scale_to"],
        model_flags=MODEL_FLAGS,
        extra=EXTRA_FLAGS
    resources:
        mem_mb=24000
    conda:
        "../envs/chipseq_peak_calling.yaml"
    log:
        OUT + "/macs3.log"
    shell:
        r"""
        mkdir -p {params.outdir:q}
        macs3 --version > {output.version:q}

        macs3 callpeak \
            -t {input.ip:q} \
            -c {input.control:q} \
            -f BAM \
            -g {params.gsize} \
            -n {params.name:q} \
            --outdir {params.outdir:q} \
            -q {params.qvalue} \
            --keep-dup {params.keep_dup:q} \
            --scale-to {params.scale_to:q} \
            {params.model_flags} \
            {params.extra} \
            > {log:q} 2>&1
        """


rule chipseq_peak_frip:
    input:
        peaks=OUT + "/peaks/" + NAME + "_peaks.narrowPeak",
        ip="inputs/" + IP + "/filtered.bam",
        ip_index="inputs/" + IP + "/filtered.bam.csi",
        control="inputs/" + CONTROL + "/filtered.bam",
        control_index="inputs/" + CONTROL + "/filtered.bam.csi"
    output:
        ip_total=OUT + "/frip/ip_total_reads.txt",
        control_total=OUT + "/frip/input_total_reads.txt",
        ip_overlap=OUT + "/frip/ip_reads_in_peaks.txt",
        control_overlap=OUT + "/frip/input_reads_overlapping_ip_peaks.txt"
    resources:
        mem_mb=4000
    conda:
        "../envs/chipseq_peak_calling.yaml"
    shell:
        r"""
        mkdir -p {OUT}/frip

        samtools view -c {input.ip:q} > {output.ip_total:q}
        samtools view -c {input.control:q} > {output.control_total:q}

        if [[ -s {input.peaks:q} ]]; then
            bedtools intersect \
                -abam {input.ip:q} \
                -b {input.peaks:q} \
                -u -bed \
                | wc -l \
                | tr -d ' ' \
                > {output.ip_overlap:q}

            bedtools intersect \
                -abam {input.control:q} \
                -b {input.peaks:q} \
                -u -bed \
                | wc -l \
                | tr -d ' ' \
                > {output.control_overlap:q}
        else
            printf '0\n' > {output.ip_overlap:q}
            printf '0\n' > {output.control_overlap:q}
        fi
        """


rule chipseq_peak_qc:
    input:
        peaks=OUT + "/peaks/" + NAME + "_peaks.narrowPeak",
        xls=OUT + "/peaks/" + NAME + "_peaks.xls",
        summits=OUT + "/peaks/" + NAME + "_summits.bed",
        treat_bdg=OUT + "/peaks/" + NAME + "_treat_pileup.bdg",
        control_bdg=OUT + "/peaks/" + NAME + "_control_lambda.bdg",
        log=OUT + "/macs3.log",
        version=OUT + "/macs3_version.txt",
        ip_total=OUT + "/frip/ip_total_reads.txt",
        control_total=OUT + "/frip/input_total_reads.txt",
        ip_overlap=OUT + "/frip/ip_reads_in_peaks.txt",
        control_overlap=OUT + "/frip/input_reads_overlapping_ip_peaks.txt",
        configuration="config/chipseq_peak_calling.json",
        plan=PLAN,
        code=CODE
    output:
        qc_json=OUT + "/peak_qc.json",
        qc_tsv=OUT + "/peak_qc.tsv",
        parameters=OUT + "/peak_calling_parameters.json",
        provenance=OUT + "/input_provenance.json"
    conda:
        "../envs/chipseq_peak_calling.yaml"
    shell:
        """
        python {input.code:q} summarize \
            --peaks {input.peaks:q} \
            --summits {input.summits:q} \
            --xls {input.xls:q} \
            --treat-bdg {input.treat_bdg:q} \
            --control-bdg {input.control_bdg:q} \
            --log {input.log:q} \
            --version {input.version:q} \
            --ip-total {input.ip_total:q} \
            --control-total {input.control_total:q} \
            --ip-overlap {input.ip_overlap:q} \
            --control-overlap {input.control_overlap:q} \
            --qc-json {output.qc_json:q} \
            --qc-tsv {output.qc_tsv:q} \
            --parameters-json {output.parameters:q} \
            --provenance-json {output.provenance:q}
        """
