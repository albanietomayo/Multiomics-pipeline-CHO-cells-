import csv
import importlib.util
import shlex
from pathlib import Path

PORTABLE_CODE = "workflow/scripts/chipseq_portable_runtime.py"
_spec = importlib.util.spec_from_file_location("chipseq_portable_dag", PORTABLE_CODE)
CHIP_PORTABLE = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(CHIP_PORTABLE)
CHIP_SETTINGS = CHIP_PORTABLE.configuration(overrides=config["chipseq"])
CHIP_PORTABLE.validate_contract()
CHIP_OUT = config["chipseq"]["output_root"]
CHIP_REF = config["chipseq"]["reference_root"]
CHIP_ROWS, _ = CHIP_PORTABLE.PRODUCTION.validated_plans()
CHIP_BY_ANALYSIS = {row["analysis_id"]: row for row in CHIP_ROWS}
CHIP_CONTROLS = sorted({row["control_run_accession"] for row in CHIP_ROWS})
CHIP_PHYSICAL = sorted({row[k] for row in CHIP_ROWS for k in ("ip_run_accession", "control_run_accession")})
CHIP_MANIFEST = CHIP_OUT + "/fastq_manifest.tsv"
CHIP_PLAN_VALID = CHIP_OUT + "/plan_validation.json"
CHIP_ENV = "../envs/chipseq_production.yaml"
CHIP_EXECUTION_ARGS = (
    "--output-root " + shlex.quote(CHIP_OUT)
    + " --reference-root " + shlex.quote(CHIP_REF)
    + " --max-parallel-analyses " + str(CHIP_SETTINGS["max_parallel_analyses"])
    + ((" --scratch-root " + shlex.quote(CHIP_SETTINGS["scratch_root"]))
       if CHIP_SETTINGS.get("scratch_root") else "")
)
CHIP_SOURCES = [PORTABLE_CODE] + ["workflow/scripts/" + n + ".py" for n in (
    "chipseq_production", "chipseq_alignment_support", "chipseq_filtering_support",
    "chipseq_peak_calling_dynamic", "chipseq_benchmark", "chipseq_download_fastq")]
CHIP_FROZEN = [config["chipseq"][k] for k in (
    "samples", "analysis_plan", "processing_plan", "peak_plan", "peak_plan_summary",
    "peak_policy", "alignment_policy", "filtering_policy")]


def chip_primary(name):
    suffix = "narrowPeak" if CHIP_BY_ANALYSIS[name]["peak_mode"] == "narrow" else "broadPeak"
    return CHIP_OUT + f"/analyses/{name}/artifacts/peaks/{name}_peaks.{suffix}"


rule chipseq_all:
    input:
        expand(CHIP_OUT + "/analyses/{analysis}/completed.json", analysis=sorted(CHIP_BY_ANALYSIS)),
        [chip_primary(n) for n in sorted(CHIP_BY_ANALYSIS)],
        CHIP_OUT + "/cohort_summary.json"
    default_target: True

rule chipseq_plan_validate:
    input:
        code=CHIP_SOURCES,
        frozen=CHIP_FROZEN,
        hashes="provenance/chipseq/production_2026-09-23/frozen_inputs.sha256",
        contract="provenance/chipseq/production_2026-09-23/frozen_contract.json"
    output: CHIP_PLAN_VALID
    conda: CHIP_ENV
    shell: "python {PORTABLE_CODE:q} validate-plan --output {output:q}"

rule chipseq_fastq_manifest:
    input: validated=CHIP_PLAN_VALID, sources=CHIP_SOURCES, frozen=CHIP_FROZEN
    output: CHIP_MANIFEST
    conda: CHIP_ENV
    shell: "python {PORTABLE_CODE:q} fastq-manifest --output {output:q}"

rule chipseq_download_all:
    input: expand(CHIP_OUT + "/raw/{run}.fastq.gz", run=CHIP_PHYSICAL)

rule chipseq_download_fastq:
    input:
        manifest=CHIP_MANIFEST,
        code="workflow/scripts/chipseq_download_fastq.py"
    output: CHIP_OUT + "/raw/{run}.fastq.gz"
    wildcard_constraints: run="|".join(CHIP_PHYSICAL)
    conda: CHIP_ENV
    params: filename=lambda w: w.run + ".fastq.gz"
    shell:
        "python {input.code:q} --manifest {input.manifest:q} --run-accession {wildcards.run:q} "
        "--filename {params.filename:q} --output {output:q}"

rule chipseq_reference_all:
    input: CHIP_REF + "/reference_inventory.json"

rule chipseq_reference_build:
    input:
        validated=CHIP_PLAN_VALID,
        code=CHIP_SOURCES,
        acquisition="workflow/scripts/fetch_reference_genome.py",
        configuration="config/config.yaml"
    output: CHIP_REF + "/reference_inventory.json"
    threads: 8
    resources: mem_mb=24000
    conda: CHIP_ENV
    params: root=CHIP_REF
    shell: "python {PORTABLE_CODE:q} reference --output {params.root:q}"

rule chipseq_shared_control:
    input:
        raw=CHIP_OUT + "/raw/{control}.fastq.gz",
        reference=CHIP_REF + "/reference_inventory.json",
        validated=CHIP_PLAN_VALID,
        code=CHIP_SOURCES, frozen=CHIP_FROZEN
    output:
        manifest=CHIP_OUT + "/shared_controls/{control}/control_manifest.json",
        bam=CHIP_OUT + "/shared_controls/{control}/filtered.bam",
        csi=CHIP_OUT + "/shared_controls/{control}/filtered.bam.csi"
    wildcard_constraints: control="|".join(CHIP_CONTROLS)
    threads: 8
    resources: mem_mb=32000, chipseq_slots=1
    conda: CHIP_ENV
    shell: "python {PORTABLE_CODE:q} {CHIP_EXECUTION_ARGS} control --control {wildcards.control:q} --raw {input.raw:q}"

rule chipseq_narrow_analysis:
    input:
        raw=lambda w: CHIP_OUT + "/raw/" + CHIP_BY_ANALYSIS[w.analysis]["ip_run_accession"] + ".fastq.gz",
        control=lambda w: CHIP_OUT + "/shared_controls/" + CHIP_BY_ANALYSIS[w.analysis]["control_run_accession"] + "/control_manifest.json",
        bam=lambda w: CHIP_OUT + "/shared_controls/" + CHIP_BY_ANALYSIS[w.analysis]["control_run_accession"] + "/filtered.bam",
        csi=lambda w: CHIP_OUT + "/shared_controls/" + CHIP_BY_ANALYSIS[w.analysis]["control_run_accession"] + "/filtered.bam.csi",
        reference=CHIP_REF + "/reference_inventory.json",
        validated=CHIP_PLAN_VALID,
        code=CHIP_SOURCES, frozen=CHIP_FROZEN
    output:
        completed=CHIP_OUT + "/analyses/{analysis}/completed.json",
        primary=CHIP_OUT + "/analyses/{analysis}/artifacts/peaks/{analysis}_peaks.narrowPeak",
        secondary=CHIP_OUT + "/analyses/{analysis}/artifacts/peaks/{analysis}_summits.bed",
        qc=CHIP_OUT + "/analyses/{analysis}/artifacts/peak_qc.json",
        parameters=CHIP_OUT + "/analyses/{analysis}/artifacts/parameters.json",
        signal=CHIP_OUT + "/analyses/{analysis}/artifacts/peaks/{analysis}_treat_pileup.bdg.gz",
        control_signal=CHIP_OUT + "/analyses/{analysis}/artifacts/peaks/{analysis}_control_lambda.bdg.gz"
    wildcard_constraints:
        analysis="|".join(sorted(n for n, r in CHIP_BY_ANALYSIS.items() if r["peak_mode"] == "narrow"))
    threads: 8
    resources: mem_mb=32000, chipseq_slots=1
    conda: CHIP_ENV
    shell: "python {PORTABLE_CODE:q} {CHIP_EXECUTION_ARGS} analysis --analysis {wildcards.analysis:q} --raw {input.raw:q}"

rule chipseq_broad_analysis:
    input:
        raw=lambda w: CHIP_OUT + "/raw/" + CHIP_BY_ANALYSIS[w.analysis]["ip_run_accession"] + ".fastq.gz",
        control=lambda w: CHIP_OUT + "/shared_controls/" + CHIP_BY_ANALYSIS[w.analysis]["control_run_accession"] + "/control_manifest.json",
        bam=lambda w: CHIP_OUT + "/shared_controls/" + CHIP_BY_ANALYSIS[w.analysis]["control_run_accession"] + "/filtered.bam",
        csi=lambda w: CHIP_OUT + "/shared_controls/" + CHIP_BY_ANALYSIS[w.analysis]["control_run_accession"] + "/filtered.bam.csi",
        reference=CHIP_REF + "/reference_inventory.json",
        validated=CHIP_PLAN_VALID,
        code=CHIP_SOURCES, frozen=CHIP_FROZEN
    output:
        completed=CHIP_OUT + "/analyses/{analysis}/completed.json",
        primary=CHIP_OUT + "/analyses/{analysis}/artifacts/peaks/{analysis}_peaks.broadPeak",
        secondary=CHIP_OUT + "/analyses/{analysis}/artifacts/peaks/{analysis}_peaks.gappedPeak",
        qc=CHIP_OUT + "/analyses/{analysis}/artifacts/peak_qc.json",
        parameters=CHIP_OUT + "/analyses/{analysis}/artifacts/parameters.json",
        signal=CHIP_OUT + "/analyses/{analysis}/artifacts/peaks/{analysis}_treat_pileup.bdg.gz",
        control_signal=CHIP_OUT + "/analyses/{analysis}/artifacts/peaks/{analysis}_control_lambda.bdg.gz"
    wildcard_constraints:
        analysis="|".join(sorted(n for n, r in CHIP_BY_ANALYSIS.items() if r["peak_mode"] == "broad"))
    threads: 8
    resources: mem_mb=32000, chipseq_slots=1
    conda: CHIP_ENV
    shell: "python {PORTABLE_CODE:q} {CHIP_EXECUTION_ARGS} analysis --analysis {wildcards.analysis:q} --raw {input.raw:q}"

rule chipseq_cohort_summary:
    input:
        completed=expand(CHIP_OUT + "/analyses/{analysis}/completed.json", analysis=sorted(CHIP_BY_ANALYSIS)),
        code=CHIP_SOURCES
    output: CHIP_OUT + "/cohort_summary.json"
    conda: CHIP_ENV
    shell: "python {PORTABLE_CODE:q} {CHIP_EXECUTION_ARGS} summary --output {output:q}"
