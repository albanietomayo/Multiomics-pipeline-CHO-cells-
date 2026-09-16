# ============================================================
# RNA-seq processing
# ============================================================

def load_rnaseq_validation_manifest():
    """
    Load RNA-seq records from the validation FASTQ manifest.

    The manifest retains the structured FASTQ roles required to distinguish
    paired, single, and unpaired sequencing inputs.
    """

    manifest = _load_validation_manifest_for_preprocessing()

    rnaseq_manifest = manifest[
        manifest["omics"] == "RNA-seq"
    ].copy()

    rnaseq_manifest = rnaseq_manifest[
        rnaseq_manifest["run_accession"].map(lambda acc: selection_decision(acc, omics="RNA-seq")[0])
    ].copy()

    if rnaseq_manifest.empty:
        raise ValueError(
            "No RNA-seq records found in the validation FASTQ manifest."
        )

    return rnaseq_manifest


def get_rnaseq_run_structure(run_accession):
    """
    Determine and validate the FASTQ structure of one RNA-seq run.
    """

    manifest = load_rnaseq_validation_manifest()

    group = manifest[
        manifest["run_accession"] == run_accession
    ]

    if group.empty:
        raise ValueError(
            f"RNA-seq run {run_accession} was not found in the validation manifest."
        )

    structures = group["fastq_structure"].dropna().unique().tolist()

    if len(structures) != 1:
        raise ValueError(
            f"Expected exactly one fastq_structure for {run_accession}, "
            f"found {structures}."
        )

    structure = structures[0]
    roles = set(group["fastq_role"])

    expected_roles = {
        "PAIRED": {"R1", "R2"},
        "SINGLE": {"SINGLE"},
        "UNPAIRED_FROM_PAIRED": {"UNPAIRED"},
        "PAIRED_PLUS_UNPAIRED": {"R1", "R2", "UNPAIRED"},
    }

    if structure not in expected_roles:
        raise ValueError(
            f"Unsupported RNA-seq FASTQ structure for {run_accession}: "
            f"{structure}"
        )

    if roles != expected_roles[structure]:
        raise ValueError(
            f"Inconsistent FASTQ roles for {run_accession}. "
            f"Structure {structure} expects "
            f"{sorted(expected_roles[structure])}, "
            f"but found {sorted(roles)}."
        )

    return structure


def get_rnaseq_validation_runs():
    """
    Return unique RNA-seq run accessions after validating their FASTQ structure.
    """

    manifest = load_rnaseq_validation_manifest()

    runs = sorted(
        manifest["run_accession"]
        .dropna()
        .unique()
        .tolist()
    )

    for run_accession in runs:
        get_rnaseq_run_structure(run_accession)

    return runs

def get_rnaseq_mate_provenance_record(run_accession):
    """
    Read the completed mate-provenance checkpoint for one RNA-seq run.
    """

    provenance_path = (
        checkpoints.rnaseq_inspect_unpaired_mate_provenance
        .get(run_accession=run_accession)
        .output.provenance
    )

    with open(
        provenance_path,
        "r",
        encoding="utf-8",
    ) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        values = handle.readline().rstrip("\n").split("\t")

    if not values or len(header) != len(values):
        raise ValueError(
            f"Invalid mate-provenance TSV for {run_accession}: "
            f"{provenance_path}"
        )

    row = dict(zip(header, values))

    if row.get("run_accession") != run_accession:
        raise ValueError(
            f"Mate-provenance run mismatch: expected "
            f"{run_accession}, found {row.get('run_accession')}"
        )

    return row


def get_rnaseq_alignment_units(run_accession):
    """
    Return effective alignment units for one RNA-seq run.

    Runs with mixed orphan R1/R2 reads are dynamically expanded into
    UNPAIRED_R1 and UNPAIRED_R2 after mate-provenance inspection.
    """

    structure = get_rnaseq_run_structure(run_accession)

    static_units = {
        "PAIRED": ["PAIRED"],
        "SINGLE": ["SINGLE"],
    }

    if structure in static_units:
        return static_units[structure]

    if structure not in {
        "UNPAIRED_FROM_PAIRED",
        "PAIRED_PLUS_UNPAIRED",
    }:
        raise ValueError(
            f"Unsupported RNA-seq FASTQ structure for "
            f"{run_accession}: {structure}"
        )

    mate_row = get_rnaseq_mate_provenance_record(
        run_accession
    )

    provenance = mate_row.get(
        "mate_provenance",
        "",
    )

    if provenance == "MIXED_R1_R2":
        orphan_units = [
            "UNPAIRED_R1",
            "UNPAIRED_R2",
        ]
    else:
        orphan_units = ["UNPAIRED"]

    if structure == "PAIRED_PLUS_UNPAIRED":
        return ["PAIRED"] + orphan_units

    return orphan_units



def get_rnaseq_preprocessed_fastq(run_accession, role):
    """
    Resolve one preprocessed RNA-seq FASTQ from its structured FASTQ role.
    """

    manifest = load_rnaseq_validation_manifest()

    match = manifest[
        (manifest["run_accession"] == run_accession)
        & (manifest["fastq_role"] == role)
    ]

    if len(match) != 1:
        raise ValueError(
            f"Expected exactly one {role} FASTQ for {run_accession}, "
            f"found {len(match)}."
        )

    fastp_dir = config["preprocessing"]["fastp_dir"]

    if role in {"R1", "R2"}:
        return (
            f"{fastp_dir}/{run_accession}/PAIRED/"
            f"{run_accession}_{role}.fastq.gz"
        )

    if role in {"SINGLE", "UNPAIRED"}:
        filename = match.iloc[0]["filename"]
        basename = filename.removesuffix(".fastq.gz")

        return (
            f"{fastp_dir}/{run_accession}/{role}/"
            f"{basename}.fastq.gz"
        )

    raise ValueError(
        f"Unsupported RNA-seq FASTQ role for {run_accession}: {role}"
    )


def get_rnaseq_alignment_reads(wildcards):
    """
    Return the preprocessed FASTQ input(s) for one STAR alignment unit.
    """

    run_accession = wildcards.run_accession
    alignment_unit = wildcards.alignment_unit

    valid_units = get_rnaseq_alignment_units(run_accession)

    if alignment_unit not in valid_units:
        raise ValueError(
            f"Alignment unit {alignment_unit} is not valid for "
            f"{run_accession}. Expected one of {valid_units}."
        )

    if alignment_unit == "PAIRED":
        return [
            get_rnaseq_preprocessed_fastq(run_accession, "R1"),
            get_rnaseq_preprocessed_fastq(run_accession, "R2"),
        ]

    if alignment_unit == "SINGLE":
        return [
            get_rnaseq_preprocessed_fastq(run_accession, "SINGLE")
        ]

    if alignment_unit == "UNPAIRED":
        return [
            get_rnaseq_preprocessed_fastq(run_accession, "UNPAIRED")
        ]

    if alignment_unit == "UNPAIRED_R1":
        return [
            (
                config["rnaseq"]["mate_split_dir"]
                + f"/{run_accession}/UNPAIRED_R1.fastq.gz"
            )
        ]

    if alignment_unit == "UNPAIRED_R2":
        return [
            (
                config["rnaseq"]["mate_split_dir"]
                + f"/{run_accession}/UNPAIRED_R2.fastq.gz"
            )
        ]

    raise ValueError(
        f"Unsupported alignment unit: {alignment_unit}"
    )


def get_rnaseq_unpaired_runs():
    """
    Return validation RNA-seq runs containing an UNPAIRED FASTQ.
    """

    runs = []

    for run_accession in get_rnaseq_validation_runs():
        structure = get_rnaseq_run_structure(run_accession)

        if structure in {
            "UNPAIRED_FROM_PAIRED",
            "PAIRED_PLUS_UNPAIRED",
        }:
            runs.append(run_accession)

    return runs


def get_rnaseq_mate_provenance_outputs():
    """
    Return mate-provenance TSV outputs for RNA-seq UNPAIRED FASTQs.
    """

    provenance_dir = config["rnaseq"]["mate_provenance_dir"]

    return [
        (
            f"{provenance_dir}/{run_accession}/"
            "mate_provenance.tsv"
        )
        for run_accession in get_rnaseq_unpaired_runs()
    ]


def get_rnaseq_effective_alignment_units_outputs():
    """
    Return effective-alignment-unit audit TSVs for all
    validation RNA-seq runs.
    """

    alignment_units_dir = config["rnaseq"]["alignment_units_dir"]

    return [
        (
            f"{alignment_units_dir}/{run_accession}/"
            "effective_alignment_units.tsv"
        )
        for run_accession in get_rnaseq_validation_runs()
    ]


def get_rnaseq_run_mate_provenance_input(wildcards):
    """
    Return mate-provenance input for runs containing UNPAIRED reads.
    """

    structure = get_rnaseq_run_structure(
        wildcards.run_accession
    )

    if structure in {
        "UNPAIRED_FROM_PAIRED",
        "PAIRED_PLUS_UNPAIRED",
    }:
        return (
            config["rnaseq"]["mate_provenance_dir"]
            + f"/{wildcards.run_accession}/mate_provenance.tsv"
        )

    return []


def get_rnaseq_run_mate_provenance_arg(wildcards):
    """
    Build the optional command-line argument for the resolver.
    """

    provenance = get_rnaseq_run_mate_provenance_input(
        wildcards
    )

    if not provenance:
        return ""

    return f"--mate-provenance {provenance}"


def get_rnaseq_star_validation_outputs():
    """
    Return the expected STAR BAM outputs for all validation RNA-seq runs.
    """

    align_dir = config["rnaseq"]["star_align_dir"]
    outputs = []

    for run_accession in get_rnaseq_validation_runs():
        for alignment_unit in get_rnaseq_alignment_units(run_accession):
            outputs.append(
                f"{align_dir}/{run_accession}/{alignment_unit}/"
                "Aligned.sortedByCoord.out.bam"
            )

    return outputs


def get_rnaseq_alignment_qc_outputs():
    """
    Return the expected alignment-QC outputs for all validation RNA-seq units.
    """

    qc_dir = config["rnaseq"]["alignment_qc_dir"]
    outputs = []

    for run_accession in get_rnaseq_validation_runs():
        for alignment_unit in get_rnaseq_alignment_units(run_accession):
            outputs.extend(
                [
                    (
                        f"{qc_dir}/{run_accession}/{alignment_unit}/"
                        "flagstat.txt"
                    ),
                    (
                        f"{qc_dir}/{run_accession}/{alignment_unit}/"
                        "quickcheck.ok"
                    ),
                    (
                        f"{qc_dir}/{run_accession}/{alignment_unit}/"
                        "star_metrics.tsv"
                    ),
                ]
            )

    return outputs


def get_rnaseq_star_metrics_outputs():
    """
    Return STAR metrics TSV files for all validation RNA-seq alignment units.
    """

    qc_dir = config["rnaseq"]["alignment_qc_dir"]
    outputs = []

    for run_accession in get_rnaseq_validation_runs():
        for alignment_unit in get_rnaseq_alignment_units(run_accession):
            outputs.append(
                f"{qc_dir}/{run_accession}/{alignment_unit}/"
                "star_metrics.tsv"
            )

    return outputs


rule rnaseq_annotation_bed12:
    input:
        gtf=config["reference"]["gtf"]

    output:
        bed=config["rnaseq"]["annotation_bed12"]

    conda:
        "../envs/rnaseq.yaml"

    shell:
        """
        mkdir -p $(dirname {output.bed})

        tmp_genepred=$(mktemp)

        gtfToGenePred \
            -ignoreGroupsWithoutExons \
            {input.gtf} \
            "$tmp_genepred"

        genePredToBed \
            "$tmp_genepred" \
            {output.bed}

        rm -f "$tmp_genepred"
        """

def get_rnaseq_strandedness_outputs():
    """
    Return RSeQC strandedness reports for all validation RNA-seq units.
    """

    strandedness_dir = config["rnaseq"]["strandedness_dir"]
    outputs = []

    for run_accession in get_rnaseq_validation_runs():
        for alignment_unit in get_rnaseq_alignment_units(run_accession):
            outputs.append(
                f"{strandedness_dir}/{run_accession}/{alignment_unit}/"
                "infer_experiment.txt"
            )

    return outputs

def get_rnaseq_strandedness_metrics_outputs():
    """
    Return structured strandedness metrics for all validation RNA-seq units.
    """

    strandedness_dir = config["rnaseq"]["strandedness_dir"]
    outputs = []

    for run_accession in get_rnaseq_validation_runs():
        for alignment_unit in get_rnaseq_alignment_units(run_accession):
            outputs.append(
                f"{strandedness_dir}/{run_accession}/{alignment_unit}/"
                "strandedness_metrics.tsv"
            )

    return outputs


def get_rnaseq_run_strandedness_classification_inputs(wildcards):
    """
    Return unit-level strandedness classifications for one RNA-seq run.
    """

    run_accession = wildcards.run_accession
    strandedness_dir = config["rnaseq"]["strandedness_dir"]

    return [
        (
            f"{strandedness_dir}/{run_accession}/{alignment_unit}/"
            "strandedness_classification.tsv"
        )
        for alignment_unit in get_rnaseq_alignment_units(run_accession)
    ]


def get_rnaseq_featurecounts_outputs():
    """
    Return featureCounts outputs for all validation RNA-seq alignment units.
    """

    counts_dir = config["rnaseq"]["counts_dir"]
    outputs = []

    for run_accession in get_rnaseq_validation_runs():
        for alignment_unit in get_rnaseq_alignment_units(run_accession):
            outputs.extend(
                [
                    (
                        f"{counts_dir}/{run_accession}/{alignment_unit}/"
                        "featurecounts.txt"
                    ),
                    (
                        f"{counts_dir}/{run_accession}/{alignment_unit}/"
                        "featurecounts.txt.summary"
                    ),
                    (
                        f"{counts_dir}/{run_accession}/{alignment_unit}/"
                        "featurecounts_strandedness.txt"
                    ),
                ]
            )

    return outputs


def get_rnaseq_run_featurecounts_inputs(wildcards):
    """
    Return featureCounts files that belong to one RNA-seq run.
    """

    run_accession = wildcards.run_accession
    counts_dir = config["rnaseq"]["counts_dir"]

    units = get_rnaseq_alignment_units(run_accession)

    return [
        (
            f"{counts_dir}/{run_accession}/{alignment_unit}/"
            "featurecounts.txt"
        )
        for alignment_unit in units
    ]

def get_rnaseq_run_counts_outputs():
    """
    Return final gene-count files for all validation RNA-seq runs.
    """

    run_counts_dir = config["rnaseq"]["run_counts_dir"]

    return [
        f"{run_counts_dir}/{run_accession}/gene_counts.tsv"
        for run_accession in get_rnaseq_validation_runs()
    ]

def get_rnaseq_featurecounts_qc_outputs():
    """
    Return structured featureCounts QC files for all validation RNA-seq units.
    """

    qc_dir = config["rnaseq"]["featurecounts_qc_dir"]
    outputs = []

    for run_accession in get_rnaseq_validation_runs():
        for alignment_unit in get_rnaseq_alignment_units(run_accession):
            outputs.append(
                f"{qc_dir}/{run_accession}/{alignment_unit}/"
                "featurecounts_metrics.tsv"
            )

    return outputs

checkpoint rnaseq_inspect_unpaired_mate_provenance:
    input:
        fastq=lambda wildcards: get_rnaseq_preprocessed_fastq(
            wildcards.run_accession,
            "UNPAIRED",
        ),
        script="workflow/scripts/inspect_unpaired_mate_provenance.py",
    output:
        provenance=(
            config["rnaseq"]["mate_provenance_dir"]
            + "/{run_accession}/mate_provenance.tsv"
        ),
    conda:
        "../envs/reporting.yaml"
    shell:
        """
        python {input.script} \
            --fastq {input.fastq} \
            --run-accession {wildcards.run_accession} \
            --output {output.provenance}
        """


rule rnaseq_mate_provenance_validation:
    input:
        lambda wildcards: get_rnaseq_mate_provenance_outputs()


rule rnaseq_split_mixed_unpaired:
    input:
        fastq=lambda wildcards: get_rnaseq_preprocessed_fastq(
            wildcards.run_accession,
            "UNPAIRED",
        ),
        provenance=lambda wildcards: (
            config["rnaseq"]["mate_provenance_dir"]
            + f"/{wildcards.run_accession}/mate_provenance.tsv"
        ),
        script="workflow/scripts/split_mixed_unpaired_fastq.py",
    output:
        r1=(
            config["rnaseq"]["mate_split_dir"]
            + "/{run_accession}/UNPAIRED_R1.fastq.gz"
        ),
        r2=(
            config["rnaseq"]["mate_split_dir"]
            + "/{run_accession}/UNPAIRED_R2.fastq.gz"
        ),
        summary=(
            config["rnaseq"]["mate_split_dir"]
            + "/{run_accession}/split_summary.tsv"
        ),
    conda:
        "../envs/reporting.yaml"
    shell:
        """
        python {input.script} \
            --input-fastq {input.fastq} \
            --run-accession {wildcards.run_accession} \
            --r1-output {output.r1} \
            --r2-output {output.r2} \
            --summary {output.summary}
        """


rule rnaseq_effective_alignment_units:
    input:
        mate_provenance=get_rnaseq_run_mate_provenance_input
    output:
        units=(
            config["rnaseq"]["alignment_units_dir"]
            + "/{run_accession}/effective_alignment_units.tsv"
        )
    run:
        import csv
        from pathlib import Path

        run_accession = wildcards.run_accession
        structure = get_rnaseq_run_structure(run_accession)

        if structure == "PAIRED":
            provenance = "NOT_APPLICABLE"
            alignment_units = ["PAIRED"]

        elif structure == "SINGLE":
            provenance = "NOT_APPLICABLE"
            alignment_units = ["SINGLE"]

        elif structure in {
            "UNPAIRED_FROM_PAIRED",
            "PAIRED_PLUS_UNPAIRED",
        }:
            provenance_path = str(input.mate_provenance)

            with open(
                provenance_path,
                "r",
                encoding="utf-8",
            ) as handle:
                rows = list(
                    csv.DictReader(
                        handle,
                        delimiter="\t",
                    )
                )

            if len(rows) != 1:
                raise ValueError(
                    f"Expected one mate-provenance record for "
                    f"{run_accession}, found {len(rows)}."
                )

            row = rows[0]

            if row["run_accession"] != run_accession:
                raise ValueError(
                    f"Mate-provenance run mismatch for "
                    f"{run_accession}."
                )

            provenance = row["mate_provenance"]

            if provenance == "MIXED_R1_R2":
                orphan_units = [
                    "UNPAIRED_R1",
                    "UNPAIRED_R2",
                ]
            else:
                orphan_units = ["UNPAIRED"]

            if structure == "PAIRED_PLUS_UNPAIRED":
                alignment_units = ["PAIRED"] + orphan_units
            else:
                alignment_units = orphan_units

        else:
            raise ValueError(
                f"Unsupported RNA-seq structure for "
                f"{run_accession}: {structure}"
            )

        output_path = Path(output.units)

        output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with output_path.open(
            "w",
            encoding="utf-8",
        ) as handle:
            handle.write(
                "run_accession\t"
                "fastq_structure\t"
                "mate_provenance\t"
                "alignment_unit\n"
            )

            for alignment_unit in alignment_units:
                handle.write(
                    f"{run_accession}\t"
                    f"{structure}\t"
                    f"{provenance}\t"
                    f"{alignment_unit}\n"
                )


rule star_index:
    input:
        fasta=config["reference"]["fasta"],
        gtf=config["reference"]["gtf"]

    output:
        index=directory(config["rnaseq"]["star_index"])

    params:
        sjdb_overhang=config["rnaseq"]["star_sjdb_overhang"]

    threads: 8

    resources:
        mem_mb=40000

    conda:
        "../envs/rnaseq.yaml"

    shell:
        """
        mkdir -p {output.index}

        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output.index} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.sjdb_overhang}
        """

rule star_align_unit:
    input:
        index=rules.star_index.output.index,
        reads=get_rnaseq_alignment_reads

    output:
        bam=(
            config["rnaseq"]["star_align_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "Aligned.sortedByCoord.out.bam"
        ),
        bai=(
            config["rnaseq"]["star_align_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "Aligned.sortedByCoord.out.bam.bai"
        ),
        log=(
            config["rnaseq"]["star_align_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "Log.final.out"
        ),
        sj=(
            config["rnaseq"]["star_align_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "SJ.out.tab"
        )

    params:
        out_dir=lambda wildcards: (
            f"{config['rnaseq']['star_align_dir']}/"
            f"{wildcards.run_accession}/"
            f"{wildcards.alignment_unit}"
        ),
        rg_id=lambda wildcards: (
            f"{wildcards.run_accession}.{wildcards.alignment_unit}"
        ),
        rg_sample=lambda wildcards: wildcards.run_accession

    threads: 8

    resources:
        mem_mb=40000

    conda:
        "../envs/rnaseq.yaml"

    wildcard_constraints:
        alignment_unit="PAIRED|SINGLE|UNPAIRED|UNPAIRED_R1|UNPAIRED_R2"

    shell:
        """
        mkdir -p {params.out_dir}

        STAR \
            --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.reads} \
            --readFilesCommand zcat \
            --outSAMattrRGline ID:{params.rg_id} SM:{params.rg_sample} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.out_dir}/

        samtools index \
            -@ {threads} \
            {output.bam} \
            {output.bai}
        """

rule rnaseq_star_validation:
    input:
        lambda wildcards: get_rnaseq_star_validation_outputs()


rule rnaseq_alignment_qc:
    input:
        bam=(
            config["rnaseq"]["star_align_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "Aligned.sortedByCoord.out.bam"
        ),
        bai=(
            config["rnaseq"]["star_align_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "Aligned.sortedByCoord.out.bam.bai"
        )

    output:
        flagstat=(
            config["rnaseq"]["alignment_qc_dir"]
            + "/{run_accession}/{alignment_unit}/flagstat.txt"
        ),
        quickcheck=(
            config["rnaseq"]["alignment_qc_dir"]
            + "/{run_accession}/{alignment_unit}/quickcheck.ok"
        )

    threads: 2

    conda:
        "../envs/rnaseq.yaml"

    shell:
        """
        mkdir -p $(dirname {output.flagstat})

        samtools quickcheck -v {input.bam}

        touch {output.quickcheck}

        samtools flagstat \
            -@ {threads} \
            {input.bam} \
            > {output.flagstat}
        """

rule rnaseq_alignment_qc_validation:
    input:
        lambda wildcards: get_rnaseq_alignment_qc_outputs()

rule summarize_star_log:
    input:
        log=(
            config["rnaseq"]["star_align_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "Log.final.out"
        ),
        script="workflow/scripts/summarize_star_log.py"

    output:
        metrics=(
            config["rnaseq"]["alignment_qc_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "star_metrics.tsv"
        )

    conda:
        "../envs/reporting.yaml"

    wildcard_constraints:
        alignment_unit="PAIRED|SINGLE|UNPAIRED|UNPAIRED_R1|UNPAIRED_R2"

    shell:
        """
        python {input.script} \
            {input.log} \
            {wildcards.run_accession} \
            {wildcards.alignment_unit} \
            {output.metrics}
        """
rule combine_star_metrics:
    input:
        metrics=lambda wildcards: get_rnaseq_star_metrics_outputs(),
        script="workflow/scripts/combine_star_metrics.py"

    output:
        summary=config["rnaseq"]["alignment_qc_summary"]

    conda:
        "../envs/reporting.yaml"

    shell:
        """
        python {input.script} \
            {output.summary} \
            {input.metrics}
        """

rule infer_rnaseq_strandedness:
    input:
        bam=(
            config["rnaseq"]["star_align_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "Aligned.sortedByCoord.out.bam"
        ),
        bai=(
            config["rnaseq"]["star_align_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "Aligned.sortedByCoord.out.bam.bai"
        ),
        bed=config["rnaseq"]["annotation_bed12"]

    output:
        report=(
            config["rnaseq"]["strandedness_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "infer_experiment.txt"
        )

    params:
        sample_size=config["rnaseq"]["rseqc_sample_size"],
        mapq=config["rnaseq"]["rseqc_mapq"]

    conda:
        "../envs/rnaseq.yaml"

    wildcard_constraints:
        alignment_unit="PAIRED|SINGLE|UNPAIRED|UNPAIRED_R1|UNPAIRED_R2"

    shell:
        """
        mkdir -p $(dirname {output.report})

        infer_experiment.py \
            -i {input.bam} \
            -r {input.bed} \
            -s {params.sample_size} \
            -q {params.mapq} \
            > {output.report}
        """

rule rnaseq_strandedness_validation:
    input:
        lambda wildcards: get_rnaseq_strandedness_outputs()

rule summarize_rseqc_strandedness:
    input:
        report=(
            config["rnaseq"]["strandedness_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "infer_experiment.txt"
        ),
        script="workflow/scripts/summarize_rseqc_strandedness.py"

    output:
        metrics=(
            config["rnaseq"]["strandedness_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "strandedness_metrics.tsv"
        )

    conda:
        "../envs/reporting.yaml"

    wildcard_constraints:
        alignment_unit="PAIRED|SINGLE|UNPAIRED|UNPAIRED_R1|UNPAIRED_R2"

    shell:
        """
        python {input.script} \
            {input.report} \
            {wildcards.run_accession} \
            {wildcards.alignment_unit} \
            {output.metrics}
        """

rule classify_rseqc_strandedness_unit:
    input:
        metrics=(
            config["rnaseq"]["strandedness_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "strandedness_metrics.tsv"
        ),
        script="workflow/scripts/classify_rseqc_strandedness.py"

    output:
        classified=(
            config["rnaseq"]["strandedness_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "strandedness_classification.tsv"
        )

    params:
        stranded_threshold=config["rnaseq"]["stranded_threshold"],
        unstranded_threshold=config["rnaseq"]["unstranded_threshold"]

    conda:
        "../envs/reporting.yaml"

    wildcard_constraints:
        alignment_unit="PAIRED|SINGLE|UNPAIRED|UNPAIRED_R1|UNPAIRED_R2"

    shell:
        """
        python {input.script} \
            {input.metrics} \
            {params.stranded_threshold} \
            {params.unstranded_threshold} \
            {output.classified}
        """

rule combine_rseqc_strandedness:
    input:
        metrics=lambda wildcards: get_rnaseq_strandedness_metrics_outputs(),
        script="workflow/scripts/combine_star_metrics.py"

    output:
        summary=config["rnaseq"]["strandedness_summary"]

    conda:
        "../envs/reporting.yaml"

    shell:
        """
        python {input.script} \
            {output.summary} \
            {input.metrics}
        """

rule classify_rseqc_strandedness:
    input:
        summary=config["rnaseq"]["strandedness_summary"],
        script="workflow/scripts/classify_rseqc_strandedness.py"

    output:
        classified=config["rnaseq"]["strandedness_classification"]

    params:
        stranded_threshold=config["rnaseq"]["stranded_threshold"],
        unstranded_threshold=config["rnaseq"]["unstranded_threshold"]

    conda:
        "../envs/reporting.yaml"

    shell:
        """
        python {input.script} \
            {input.summary} \
            {params.stranded_threshold} \
            {params.unstranded_threshold} \
            {output.classified}
        """

rule resolve_rseqc_strandedness_run:
    input:
        classifications=get_rnaseq_run_strandedness_classification_inputs,
        mate_provenance=get_rnaseq_run_mate_provenance_input,
        script="workflow/scripts/resolve_rseqc_strandedness_run.py"
    output:
        resolved=(
            config["rnaseq"]["strandedness_dir"]
            + "/{run_accession}/resolved_strandedness.tsv"
        )
    params:
        mate_provenance_arg=get_rnaseq_run_mate_provenance_arg
    conda:
        "../envs/reporting.yaml"
    shell:
        """
        python {input.script} \
            {wildcards.run_accession} \
            {output.resolved} \
            {params.mate_provenance_arg} \
            {input.classifications}
        """



rule featurecounts_primary:
    input:
        bam=(
            config["rnaseq"]["star_align_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "Aligned.sortedByCoord.out.bam"
        ),
        gtf=config["reference"]["gtf"],
        classification=(
            config["rnaseq"]["strandedness_dir"]
            + "/{run_accession}/resolved_strandedness.tsv"
        ),
        script="workflow/scripts/run_featurecounts.py",
        recovery="workflow/scripts/featurecounts_recovery.py"

    output:
        attempt=directory(config["rnaseq"]["counts_dir"] + "/{run_accession}/{alignment_unit}/.primary")
    threads: 4
    conda:
        "../envs/rnaseq.yaml"
    wildcard_constraints:
        alignment_unit="PAIRED|SINGLE|UNPAIRED|UNPAIRED_R1|UNPAIRED_R2"
    shell:
        """
        python {input.recovery:q} primary {output.attempt:q} {input.script:q} \
            {input.classification:q} {wildcards.run_accession:q} {wildcards.alignment_unit:q} \
            {input.gtf:q} {input.bam:q} {threads} unused unused unused
        """

rule featurecounts_unit:
    input:
        bam=(
            config["rnaseq"]["star_align_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "Aligned.sortedByCoord.out.bam"
        ),
        gtf=config["reference"]["gtf"],
        classification=(
            config["rnaseq"]["strandedness_dir"]
            + "/{run_accession}/resolved_strandedness.tsv"
        ),
        script="workflow/scripts/run_featurecounts.py",
        recovery="workflow/scripts/featurecounts_recovery.py",
        attempt=config["rnaseq"]["counts_dir"] + "/{run_accession}/{alignment_unit}/.primary"

    output:
        counts=(
            config["rnaseq"]["counts_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "featurecounts.txt"
        ),
        summary=(
            config["rnaseq"]["counts_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "featurecounts.txt.summary"
        ),
        strandedness=(
            config["rnaseq"]["counts_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "featurecounts_strandedness.txt"
        ),
        provenance=config["rnaseq"]["counts_dir"] + "/{run_accession}/{alignment_unit}/featurecounts_provenance.json"


    threads: 1

    conda:
        "../envs/featurecounts_fallback.yaml"

    wildcard_constraints:
        alignment_unit="PAIRED|SINGLE|UNPAIRED|UNPAIRED_R1|UNPAIRED_R2"

    shell:
        """
        python {input.recovery:q} finalize {input.attempt:q} {input.script:q} \
            {input.classification:q} {wildcards.run_accession:q} {wildcards.alignment_unit:q} \
            {input.gtf:q} {input.bam:q} {threads} \
            {output.counts:q} {output.strandedness:q} {output.provenance:q}
        """

rule rnaseq_featurecounts_validation:
    input:
        lambda wildcards: get_rnaseq_featurecounts_outputs()

rule combine_featurecounts_run:
    input:
        counts=get_rnaseq_run_featurecounts_inputs,
        script="workflow/scripts/combine_featurecounts_units.py"

    output:
        combined=(
            config["rnaseq"]["run_counts_dir"]
            + "/{run_accession}/gene_counts.tsv"
        )

    conda:
        "../envs/reporting.yaml"

    shell:
        """
        python {input.script} \
            {wildcards.run_accession} \
            {output.combined} \
            {input.counts}
        """

rule rnaseq_run_counts_validation:
    input:
        lambda wildcards: get_rnaseq_run_counts_outputs()


rule build_rnaseq_count_matrix:
    input:
        counts=lambda wildcards: get_rnaseq_run_counts_outputs(),
        script="workflow/scripts/build_rnaseq_count_matrix.py",
        selection=list(SELECTION_FILES)

    output:
        matrix=config["rnaseq"]["counts_matrix"],
        selection=config["rnaseq"]["counts_matrix"] + ".selection.tsv"

    conda:
        "../envs/reporting.yaml"

    shell:
        """
        python {input.script:q} \
            {output.matrix:q} \
            {input.counts:q}
        """

rule rnaseq_validation:
    input:
        counts_matrix=config["rnaseq"]["counts_matrix"],
        alignment_summary=config["rnaseq"]["alignment_qc_summary"],
        alignment_qc=lambda wildcards: get_rnaseq_alignment_qc_outputs(),
        strandedness_summary=config["rnaseq"]["strandedness_summary"],
        strandedness_classification=config["rnaseq"]["strandedness_classification"],
        mate_provenance=lambda wildcards: get_rnaseq_mate_provenance_outputs(),
        effective_alignment_units=lambda wildcards: get_rnaseq_effective_alignment_units_outputs(),
        featurecounts_qc_summary=config["rnaseq"]["featurecounts_qc_summary"]

rule summarize_featurecounts_unit:
    input:
        summary=(
            config["rnaseq"]["counts_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "featurecounts.txt.summary"
        ),
        script="workflow/scripts/summarize_featurecounts.py"

    output:
        metrics=(
            config["rnaseq"]["featurecounts_qc_dir"]
            + "/{run_accession}/{alignment_unit}/"
            + "featurecounts_metrics.tsv"
        )

    conda:
        "../envs/reporting.yaml"

    wildcard_constraints:
        alignment_unit="PAIRED|SINGLE|UNPAIRED|UNPAIRED_R1|UNPAIRED_R2"

    shell:
        """
        python {input.script} \
            {input.summary} \
            {wildcards.run_accession} \
            {wildcards.alignment_unit} \
            {output.metrics}
        """


rule combine_featurecounts_qc:
    input:
        metrics=lambda wildcards: get_rnaseq_featurecounts_qc_outputs(),
        script="workflow/scripts/combine_star_metrics.py"

    output:
        summary=config["rnaseq"]["featurecounts_qc_summary"]

    conda:
        "../envs/reporting.yaml"

    shell:
        """
        python {input.script} \
            {output.summary} \
            {input.metrics}
        """

# Explicit production targets: reconcile saved counts, then build the eligible matrix.
# Existing default validation targets remain separate from the production matrix.
def production_runs_root(wildcards):
    root = config.get("rnaseq_production_runs_root", "")
    if not root:
        raise ValueError("Set --config rnaseq_production_runs_root=/absolute/path/to/production/runs")
    return root

checkpoint rnaseq_production_inventory:
    input:
        runs_root=production_runs_root,
        selection=list(SELECTION_FILES),
        script="workflow/scripts/inventory_selection_results.py"
    output:
        counts="results/rnaseq/production_selection/eligible_counts.list",
        summary="results/rnaseq/production_selection/summary.json",
        runs="results/rnaseq/production_selection/run_status.tsv",
        missing="results/rnaseq/production_selection/missing_or_unverified.tsv",
        files="results/rnaseq/production_selection/excluded_files.tsv"
    shell:
        """
        python {input.script:q} --runs-root {input.runs_root:q} \\
          --outdir results/rnaseq/production_selection --refresh
        """

def production_count_files(wildcards):
    checked = checkpoints.rnaseq_production_inventory.get()
    with open(checked.output.counts) as handle:
        return [line.strip() for line in handle if line.strip()]

rule rnaseq_production_count_matrix:
    input:
        counts=production_count_files,
        counts_list=lambda wildcards: checkpoints.rnaseq_production_inventory.get().output.counts,
        selection=list(SELECTION_FILES),
        script="workflow/scripts/build_rnaseq_count_matrix.py"
    output:
        matrix="results/rnaseq/production_matrix/gene_counts.tsv",
        selection="results/rnaseq/production_matrix/gene_counts.tsv.selection.tsv",
        provenance="results/rnaseq/production_matrix/gene_counts.tsv.provenance.json"
    shell:
        """
        python {input.script:q} {output.matrix:q} --counts-list {input.counts_list:q}
        """
