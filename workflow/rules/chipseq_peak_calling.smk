import json

RUNTIME_PATH = "config/chipseq_peak_calling_runtime.json"
CODE = "workflow/scripts/chipseq_peak_calling_dynamic.py"
OUT = "results/chipseq/peak_calling"

with open(RUNTIME_PATH) as handle:
    RUNTIME = json.load(handle)

if (RUNTIME.get("schema_version") != 1
        or RUNTIME.get("scope") != "dynamic_chipseq_peak_calling_runtime"
        or RUNTIME.get("full_bam_and_csi_hashes_verified") is not True):
    raise ValueError("Unverified or unsupported dynamic peak-calling runtime")

ROWS = RUNTIME.get("analyses", [])
if len(ROWS) != 18 or RUNTIME.get("analysis_count") != 18:
    raise ValueError("Dynamic peak-calling runtime must contain exactly 18 analyses")
BY_ANALYSIS = {row["analysis_id"]: row for row in ROWS}
if len(BY_ANALYSIS) != len(ROWS):
    raise ValueError("Duplicate analysis_id in dynamic peak-calling runtime")

NARROW = sorted(row["analysis_id"] for row in ROWS if row["peak_mode"] == "narrow")
BROAD = sorted(row["analysis_id"] for row in ROWS if row["peak_mode"] == "broad")
EXTERNAL = sorted(row["analysis_id"] for row in ROWS
                  if row["fragment_size_policy"] == "phantompeakqualtools")
FIXED = sorted(row["analysis_id"] for row in ROWS
               if row["fragment_size_policy"] == "fixed")
if (len(NARROW), len(BROAD), len(EXTERNAL), len(FIXED)) != (8, 10, 12, 6):
    raise ValueError("Unexpected peak mode or fragment-policy counts")
if any(BY_ANALYSIS[name]["study_accession"] != "PRJEB9291" for name in EXTERNAL):
    raise ValueError("External fragment estimation leaked outside PRJEB9291")
if any(BY_ANALYSIS[name]["study_accession"] != "PRJNA865478" for name in FIXED):
    raise ValueError("Fixed fragment policy leaked outside PRJNA865478")


def run_item(analysis, role):
    row = BY_ANALYSIS[analysis]
    run = row["ip_run_accession"] if role == "ip" else row["control_run_accession"]
    return RUNTIME["runs"][run]


def primary_peak(analysis):
    suffix = "narrowPeak" if BY_ANALYSIS[analysis]["peak_mode"] == "narrow" else "broadPeak"
    return f"{OUT}/{analysis}/peaks/{analysis}_peaks.{suffix}"


def secondary_peak(analysis):
    suffix = "summits.bed" if BY_ANALYSIS[analysis]["peak_mode"] == "narrow" else "peaks.gappedPeak"
    return f"{OUT}/{analysis}/peaks/{analysis}_{suffix}"


rule chipseq_peak_calling_all:
    input:
        [primary_peak(name) for name in sorted(BY_ANALYSIS)],
        [secondary_peak(name) for name in sorted(BY_ANALYSIS)],
        expand(OUT + "/{analysis}/peak_qc.json", analysis=sorted(BY_ANALYSIS)),
        expand(OUT + "/{analysis}/parameters.json", analysis=sorted(BY_ANALYSIS)),
        expand(OUT + "/{analysis}/provenance.json", analysis=sorted(BY_ANALYSIS)),
        OUT + "/cohort_summary.json"
    default_target: True


rule chipseq_phantompeakqualtools:
    input:
        bam=lambda w: run_item(w.analysis, "ip")["bam_path"],
        csi=lambda w: run_item(w.analysis, "ip")["csi_path"],
        runtime=RUNTIME_PATH
    output:
        table=OUT + "/{analysis}/fragment/phantompeakqualtools.tsv",
        plot=OUT + "/{analysis}/fragment/cross_correlation.pdf"
    wildcard_constraints:
        analysis="|".join(EXTERNAL)
    threads: 4
    conda: "../envs/chipseq_phantompeakqualtools.yaml"
    log: OUT + "/{analysis}/fragment/phantompeakqualtools.log"
    shell:
        r"""
        mkdir -p $(dirname {output.table:q})
        CHIP_PPQT_TABLE_TMP={output.table:q}.tmp.$$
        CHIP_PPQT_PLOT_TMP={output.plot:q}.tmp.$$
        trap 'rm -f "$CHIP_PPQT_TABLE_TMP" "$CHIP_PPQT_PLOT_TMP"' EXIT
        run_spp.R -c={input.bam:q} -p={threads} -savp="$CHIP_PPQT_PLOT_TMP" -out="$CHIP_PPQT_TABLE_TMP" > {log:q} 2>&1
        test -s "$CHIP_PPQT_TABLE_TMP"
        test -s "$CHIP_PPQT_PLOT_TMP"
        mv -f "$CHIP_PPQT_TABLE_TMP" {output.table:q}
        mv -f "$CHIP_PPQT_PLOT_TMP" {output.plot:q}
        trap - EXIT
        """


rule chipseq_parse_phantompeakqualtools:
    input:
        table=OUT + "/{analysis}/fragment/phantompeakqualtools.tsv",
        runtime=RUNTIME_PATH,
        code=CODE
    output: OUT + "/{analysis}/fragment/fragment_estimate.json"
    wildcard_constraints: analysis="|".join(EXTERNAL)
    conda: "../envs/chipseq_peak_calling.yaml"
    shell:
        "python {input.code:q} parse-phantom --table {input.table:q} --runtime {input.runtime:q} "
        "--analysis {wildcards.analysis:q} --output {output:q}"


rule chipseq_peak_parameters:
    input:
        runtime=RUNTIME_PATH,
        code=CODE,
        fragment=lambda w: (OUT + f"/{w.analysis}/fragment/fragment_estimate.json"
                            if w.analysis in EXTERNAL else [])
    output: OUT + "/{analysis}/parameters.json"
    wildcard_constraints: analysis="|".join(sorted(BY_ANALYSIS))
    conda: "../envs/chipseq_peak_calling.yaml"
    params:
        fragment_arg=lambda w: ("--fragment " + OUT + f"/{w.analysis}/fragment/fragment_estimate.json"
                                if w.analysis in EXTERNAL else "")
    shell:
        "python {input.code:q} analysis-parameters --runtime {input.runtime:q} "
        "--analysis {wildcards.analysis:q} {params.fragment_arg} --output {output:q}"


rule chipseq_call_narrow_peaks:
    input:
        runtime=RUNTIME_PATH,
        parameters=OUT + "/{analysis}/parameters.json",
        fragment=lambda w: (OUT + f"/{w.analysis}/fragment/fragment_estimate.json"
                            if w.analysis in EXTERNAL else []),
        code=CODE
    output:
        peaks=OUT + "/{analysis}/peaks/{analysis}_peaks.narrowPeak",
        summits=OUT + "/{analysis}/peaks/{analysis}_summits.bed",
        xls=OUT + "/{analysis}/peaks/{analysis}_peaks.xls",
        treat_bdg=OUT + "/{analysis}/peaks/{analysis}_treat_pileup.bdg",
        control_bdg=OUT + "/{analysis}/peaks/{analysis}_control_lambda.bdg",
        version=OUT + "/{analysis}/macs3_version.txt"
    wildcard_constraints: analysis="|".join(NARROW)
    resources: mem_mb=24000
    conda: "../envs/chipseq_peak_calling.yaml"
    log: OUT + "/{analysis}/macs3.log"
    params:
        outdir=OUT + "/{analysis}/peaks",
        fragment_arg=lambda w: ("--fragment " + OUT + f"/{w.analysis}/fragment/fragment_estimate.json"
                                if w.analysis in EXTERNAL else "")
    shell:
        "python {input.code:q} call-macs3 --runtime {input.runtime:q} --parameters {input.parameters:q} "
        "--analysis {wildcards.analysis:q} --output-dir {params.outdir:q} "
        "--log {log:q} --version {output.version:q} {params.fragment_arg}"


rule chipseq_call_broad_peaks:
    input:
        runtime=RUNTIME_PATH,
        parameters=OUT + "/{analysis}/parameters.json",
        fragment=lambda w: (OUT + f"/{w.analysis}/fragment/fragment_estimate.json"
                            if w.analysis in EXTERNAL else []),
        code=CODE
    output:
        peaks=OUT + "/{analysis}/peaks/{analysis}_peaks.broadPeak",
        gapped=OUT + "/{analysis}/peaks/{analysis}_peaks.gappedPeak",
        xls=OUT + "/{analysis}/peaks/{analysis}_peaks.xls",
        treat_bdg=OUT + "/{analysis}/peaks/{analysis}_treat_pileup.bdg",
        control_bdg=OUT + "/{analysis}/peaks/{analysis}_control_lambda.bdg",
        version=OUT + "/{analysis}/macs3_version.txt"
    wildcard_constraints: analysis="|".join(BROAD)
    resources: mem_mb=24000
    conda: "../envs/chipseq_peak_calling.yaml"
    log: OUT + "/{analysis}/macs3.log"
    params:
        outdir=OUT + "/{analysis}/peaks",
        fragment_arg=lambda w: ("--fragment " + OUT + f"/{w.analysis}/fragment/fragment_estimate.json"
                                if w.analysis in EXTERNAL else "")
    shell:
        "python {input.code:q} call-macs3 --runtime {input.runtime:q} --parameters {input.parameters:q} "
        "--analysis {wildcards.analysis:q} --output-dir {params.outdir:q} "
        "--log {log:q} --version {output.version:q} {params.fragment_arg}"


rule chipseq_peak_frip:
    input:
        peaks=lambda w: primary_peak(w.analysis),
        ip=lambda w: run_item(w.analysis, "ip")["bam_path"],
        ip_index=lambda w: run_item(w.analysis, "ip")["csi_path"],
        control=lambda w: run_item(w.analysis, "input")["bam_path"],
        control_index=lambda w: run_item(w.analysis, "input")["csi_path"]
    output:
        ip_total=OUT + "/{analysis}/frip/ip_total_reads.txt",
        input_total=OUT + "/{analysis}/frip/input_total_reads.txt",
        ip_overlap=OUT + "/{analysis}/frip/ip_reads_in_peaks.txt",
        input_overlap=OUT + "/{analysis}/frip/input_reads_overlapping_ip_peaks.txt"
    resources: mem_mb=4000
    conda: "../envs/chipseq_peak_calling.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.ip_total:q})
        samtools view -c {input.ip:q} > {output.ip_total:q}
        samtools view -c {input.control:q} > {output.input_total:q}
        if [[ -s {input.peaks:q} ]]; then
            bedtools intersect -abam {input.ip:q} -b {input.peaks:q} -u -bed | wc -l | tr -d ' ' > {output.ip_overlap:q}
            bedtools intersect -abam {input.control:q} -b {input.peaks:q} -u -bed | wc -l | tr -d ' ' > {output.input_overlap:q}
        else
            printf '0\n' > {output.ip_overlap:q}
            printf '0\n' > {output.input_overlap:q}
        fi
        """


rule chipseq_peak_qc:
    input:
        runtime=RUNTIME_PATH,
        parameters=OUT + "/{analysis}/parameters.json",
        primary=lambda w: primary_peak(w.analysis),
        secondary=lambda w: secondary_peak(w.analysis),
        xls=OUT + "/{analysis}/peaks/{analysis}_peaks.xls",
        treat_bdg=OUT + "/{analysis}/peaks/{analysis}_treat_pileup.bdg",
        control_bdg=OUT + "/{analysis}/peaks/{analysis}_control_lambda.bdg",
        log=OUT + "/{analysis}/macs3.log",
        version=OUT + "/{analysis}/macs3_version.txt",
        ip_total=OUT + "/{analysis}/frip/ip_total_reads.txt",
        input_total=OUT + "/{analysis}/frip/input_total_reads.txt",
        ip_overlap=OUT + "/{analysis}/frip/ip_reads_in_peaks.txt",
        input_overlap=OUT + "/{analysis}/frip/input_reads_overlapping_ip_peaks.txt",
        fragment=lambda w: (OUT + f"/{w.analysis}/fragment/fragment_estimate.json"
                            if w.analysis in EXTERNAL else []),
        code=CODE
    output:
        qc=OUT + "/{analysis}/peak_qc.json",
        provenance=OUT + "/{analysis}/provenance.json"
    conda: "../envs/chipseq_peak_calling.yaml"
    params:
        fragment_arg=lambda w: ("--fragment " + OUT + f"/{w.analysis}/fragment/fragment_estimate.json"
                                if w.analysis in EXTERNAL else "")
    shell:
        "python {input.code:q} summarize-analysis --runtime {input.runtime:q} --parameters {input.parameters:q} "
        "--analysis {wildcards.analysis:q} --primary-peak {input.primary:q} --secondary-peak {input.secondary:q} "
        "--xls {input.xls:q} --treat-bdg {input.treat_bdg:q} --control-bdg {input.control_bdg:q} "
        "--log {input.log:q} --version {input.version:q} --ip-total {input.ip_total:q} "
        "--input-total {input.input_total:q} --ip-overlap {input.ip_overlap:q} "
        "--input-overlap {input.input_overlap:q} --qc {output.qc:q} --provenance {output.provenance:q} "
        "{params.fragment_arg}"


rule chipseq_peak_cohort_summary:
    input:
        records=expand(OUT + "/{analysis}/peak_qc.json", analysis=sorted(BY_ANALYSIS)),
        code=CODE
    output: OUT + "/cohort_summary.json"
    conda: "../envs/chipseq_peak_calling.yaml"
    shell: "python {input.code:q} summarize-cohort --output {output:q} {input.records:q}"
