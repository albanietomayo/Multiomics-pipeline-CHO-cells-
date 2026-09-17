import json
import re
from pathlib import PurePosixPath

with open(
    "config/chipseq_filtering.json"
) as handle:
    CF = json.load(handle)

with open(
    "config/chipseq_filtering_inputs.json"
) as handle:
    FILTER_INPUTS = json.load(handle)


def _safe_relative(value):
    path = PurePosixPath(value)

    return (
        bool(path.parts)
        and not path.is_absolute()
        and ".." not in path.parts
    )


if FILTER_INPUTS.get("schema_version") != 1:
    raise ValueError(
        "Unsupported filtering input-plan schema"
    )

if FILTER_INPUTS.get("library_layout") != "SINGLE":
    raise ValueError(
        "Dynamic ChIP-seq filtering currently supports "
        "SINGLE only"
    )

if FILTER_INPUTS.get("instrument_platform") != "ILLUMINA":
    raise ValueError(
        "Dynamic ChIP-seq filtering currently supports "
        "ILLUMINA only"
    )

SOURCE_ROOT = FILTER_INPUTS.get(
    "source_root"
)

if (
    not isinstance(SOURCE_ROOT, str)
    or not SOURCE_ROOT.strip()
):
    raise ValueError(
        "Filtering input plan has no source_root"
    )

RUN_ROWS = FILTER_INPUTS.get(
    "runs"
)

if (
    not isinstance(RUN_ROWS, list)
    or not RUN_ROWS
):
    raise ValueError(
        "Filtering input plan contains no runs"
    )

RUNS = []
ROLE_COUNTS = {
    "ip": 0,
    "input": 0,
}

for item in RUN_ROWS:

    if not isinstance(item, dict):
        raise ValueError(
            "Filtering run entry must be an object"
        )

    run = item.get(
        "run_accession"
    )

    if (
        not isinstance(run, str)
        or not re.fullmatch(
            r"[DES]RR[0-9]+",
            run,
        )
    ):
        raise ValueError(
            f"Invalid filtering run accession: {run!r}"
        )

    if run in RUNS:
        raise ValueError(
            f"Duplicate filtering run: {run}"
        )

    role = item.get(
        "role"
    )

    if role not in ROLE_COUNTS:
        raise ValueError(
            f"Invalid filtering role for {run}: {role!r}"
        )

    relative = item.get(
        "source_relative"
    )

    expected_relative = (
        f"outputs/{run}/raw.sorted.bam"
    )

    if (
        relative != expected_relative
        or not _safe_relative(relative)
    ):
        raise ValueError(
            f"Unexpected upstream BAM path for {run}: "
            f"{relative!r}"
        )

    digest = item.get(
        "sha256"
    )

    if (
        not isinstance(digest, str)
        or not re.fullmatch(
            r"[0-9a-f]{64}",
            digest,
        )
    ):
        raise ValueError(
            f"Invalid upstream BAM SHA-256 for {run}"
        )

    try:
        bam_bytes = int(
            item["bytes"]
        )

        input_reads = int(
            item["input_reads"]
        )

        mapped_reads = int(
            item["mapped_reads"]
        )

        nuclear_mapq = int(
            item[
                "nuclear_mapq_ge_threshold"
            ]
        )

    except (
        KeyError,
        TypeError,
        ValueError,
    ) as exc:

        raise ValueError(
            f"Invalid filtering counts for {run}"
        ) from exc

    if bam_bytes <= 0:
        raise ValueError(
            f"Empty upstream BAM for {run}"
        )

    if input_reads <= 0:
        raise ValueError(
            f"Non-positive upstream read count for {run}"
        )

    if not (
        0
        <= mapped_reads
        <= input_reads
    ):
        raise ValueError(
            f"Invalid mapped-read count for {run}"
        )

    if not (
        0
        <= nuclear_mapq
        <= mapped_reads
    ):
        raise ValueError(
            f"Invalid nuclear MAPQ count for {run}"
        )

    RUNS.append(
        run
    )

    ROLE_COUNTS[
        role
    ] += 1

if ROLE_COUNTS["ip"] < 1:
    raise ValueError(
        "Filtering cohort contains no IP runs"
    )

if ROLE_COUNTS["input"] < 1:
    raise ValueError(
        "Filtering cohort contains no Input runs"
    )


OUT = "results/chipseq/filtering"
CODE = "workflow/scripts/chipseq_filtering_support.py"
PLAN = "config/chipseq_filtering_inputs.json"

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
        reference="inputs/reference/genome_plus_mt.fa.fai"
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
