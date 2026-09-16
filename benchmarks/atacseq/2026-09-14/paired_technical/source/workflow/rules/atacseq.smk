# ============================================================
# ATAC-seq
# ============================================================

rule atacseq_reference:
    input:
        script="workflow/scripts/build_atac_reference.py",
        nuclear_fasta=config["reference"]["fasta"]

    output:
        mitochondrial_fasta=config["atacseq"]["mitochondrial_fasta"],
        mapping_fasta=config["atacseq"]["mapping_fasta"],
        metadata=config["atacseq"]["reference_metadata"],
        sha256=config["atacseq"]["reference_sha256"]

    params:
        mitochondrial_accession=config["atacseq"]["mitochondrial_accession"]

    conda:
        "../envs/reference.yaml"

    shell:
        """
        python {input.script} \
          --nuclear-fasta {input.nuclear_fasta} \
          --mitochondrial-accession {params.mitochondrial_accession} \
          --mitochondrial-fasta {output.mitochondrial_fasta} \
          --mapping-fasta {output.mapping_fasta} \
          --metadata {output.metadata} \
          --sha256 {output.sha256}
        """



rule atacseq_bowtie2_index:
    input:
        fasta=config["atacseq"]["mapping_fasta"]

    output:
        index_dir=directory(config["atacseq"]["bowtie2_index_dir"])

    params:
        prefix=config["atacseq"]["bowtie2_index"]

    threads:
        8

    conda:
        "../envs/atacseq.yaml"

    shell:
        """
        mkdir -p {output.index_dir}
        bowtie2-build --threads {threads} {input.fasta:q} {params.prefix:q}
        """



# ------------------------------------------------------------
# ATAC-seq validation samples
# ------------------------------------------------------------

def load_atacseq_validation_manifest():
    """
    Load ATAC-seq records from the validation FASTQ manifest.
    """
    manifest = _load_validation_manifest_for_preprocessing()

    atacseq_manifest = manifest[
        manifest["omics"] == "ATAC-seq"
    ].copy()

    if atacseq_manifest.empty:
        raise ValueError(
            "No ATAC-seq records found in the validation FASTQ manifest."
        )

    return atacseq_manifest


def get_atacseq_run_structure(run_accession):
    """
    Determine and validate the FASTQ structure of one ATAC-seq run.
    """
    manifest = load_atacseq_validation_manifest()

    group = manifest[
        manifest["run_accession"] == run_accession
    ]

    if group.empty:
        raise ValueError(
            f"ATAC-seq run {run_accession} was not found in the validation manifest."
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
    }

    if structure not in expected_roles:
        raise ValueError(
            f"Unsupported ATAC-seq FASTQ structure for {run_accession}: "
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


def get_atacseq_preprocessed_fastq(run_accession, role):
    """
    Resolve one preprocessed ATAC-seq FASTQ from its structured FASTQ role.
    """
    manifest = load_atacseq_validation_manifest()

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

    if role == "SINGLE":
        filename = match.iloc[0]["filename"]
        basename = filename.removesuffix(".fastq.gz")

        return (
            f"{fastp_dir}/{run_accession}/SINGLE/"
            f"{basename}.fastq.gz"
        )

    raise ValueError(
        f"Unsupported ATAC-seq FASTQ role for {run_accession}: {role}"
    )


def get_atacseq_alignment_reads(wildcards):
    """
    Return the FASTQ input(s) required for one ATAC-seq alignment.
    """
    run_accession = wildcards.run_accession
    structure = get_atacseq_run_structure(run_accession)

    if structure == "PAIRED":
        return [
            get_atacseq_preprocessed_fastq(run_accession, "R1"),
            get_atacseq_preprocessed_fastq(run_accession, "R2"),
        ]

    return [
        get_atacseq_preprocessed_fastq(run_accession, "SINGLE")
    ]


def get_atacseq_bowtie2_read_args(wildcards):
    """
    Build Bowtie2 FASTQ arguments according to the run structure.
    """
    run_accession = wildcards.run_accession
    structure = get_atacseq_run_structure(run_accession)

    if structure == "PAIRED":
        r1 = get_atacseq_preprocessed_fastq(run_accession, "R1")
        r2 = get_atacseq_preprocessed_fastq(run_accession, "R2")
        max_fragment = config["atacseq"]["paired_max_fragment_length"]

        return f"-1 {r1} -2 {r2} -X {max_fragment}"

    single = get_atacseq_preprocessed_fastq(run_accession, "SINGLE")
    return f"-U {single}"


# ------------------------------------------------------------
# Bowtie2 alignment
# ------------------------------------------------------------

def _hpc_atac_output(path):
    return temp(path) if config.get("hpc_validation_cleanup", False) else path


rule atacseq_align_run:
    input:
        index_dir=config["atacseq"]["bowtie2_index_dir"],
        reads=get_atacseq_alignment_reads

    output:
        bam=_hpc_atac_output(config["atacseq"]["align_dir"] + "/{run_accession}/raw.sorted.bam"),
        bai=_hpc_atac_output(config["atacseq"]["align_dir"] + "/{run_accession}/raw.sorted.bam.bai"),
        log=config["atacseq"]["align_dir"] + "/{run_accession}/bowtie2.log"

    params:
        index=config["atacseq"]["bowtie2_index"],
        read_args=get_atacseq_bowtie2_read_args,
        mode="--" + config["atacseq"]["bowtie2_mode"],
        preset="--" + config["atacseq"]["bowtie2_preset"]

    threads:
        8

    resources:
        mem_mb=8000

    conda:
        "../envs/atacseq.yaml"

    shell:
        """
        mkdir -p $(dirname {output.bam:q})

        bowtie2 \
          {params.mode} \
          {params.preset} \
          --threads {threads} \
          --rg-id {wildcards.run_accession} \
          --rg SM:{wildcards.run_accession} \
          -x {params.index:q} \
          {params.read_args} \
          2> {output.log:q} \
        | samtools sort -o {output.bam:q} -

        samtools index {output.bam:q} {output.bai:q}
        """

rule atacseq_mark_duplicates:
    input:
        bam=config["atacseq"]["align_dir"] + "/{run_accession}/raw.sorted.bam"

    output:
        bam=_hpc_atac_output(config["atacseq"]["align_dir"] + "/{run_accession}/dupmarked.sorted.bam"),
        bai=_hpc_atac_output(config["atacseq"]["align_dir"] + "/{run_accession}/dupmarked.sorted.bai"),
        metrics=config["atacseq"]["qc_dir"] + "/{run_accession}/duplication_metrics.txt"

    resources:
        mem_mb=4000

    conda:
        "../envs/atacseq_picard.yaml"

    shell:
        """
        mkdir -p $(dirname {output.bam:q})
        mkdir -p $(dirname {output.metrics:q})

        picard_tmp=$(mktemp -d {resources.tmpdir:q}/picard_markduplicates.XXXXXX)
        trap 'rm -rf -- "$picard_tmp"' EXIT

        printf 'PICARD_TMP_DIR=%s\n' "$picard_tmp"
        df -h "$picard_tmp"
        df -i "$picard_tmp"

        JAVA_TOOL_OPTIONS="${{JAVA_TOOL_OPTIONS:-}} -Djava.io.tmpdir=$picard_tmp" \
        picard -Xmx3g MarkDuplicates \
          I={input.bam:q} \
          O={output.bam:q} \
          M={output.metrics:q} \
          TMP_DIR="$picard_tmp" \
          REMOVE_DUPLICATES=false \
          CREATE_INDEX=true

        test -s {output.bai:q}
        """


rule atacseq_filter_bam:
    input:
        manifest=lambda wildcards: (
            checkpoints.validation_fastq_manifest
            .get()
            .output.manifest
        ),
        bam=config["atacseq"]["align_dir"] + "/{run_accession}/dupmarked.sorted.bam",
        bai=config["atacseq"]["align_dir"] + "/{run_accession}/dupmarked.sorted.bai"

    output:
        bam=config["atacseq"]["filtered_bam_dir"] + "/{run_accession}/filtered.bam",
        bai=config["atacseq"]["filtered_bam_dir"] + "/{run_accession}/filtered.bam.bai"

    params:
        structure=lambda wildcards: get_atacseq_run_structure(
            wildcards.run_accession
        ),
        min_mapq=config["atacseq"]["min_mapq"],
        single_exclude_flags=config["atacseq"]["single_exclude_flags"],
        paired_require_flags=config["atacseq"]["paired_require_flags"],
        paired_exclude_flags=config["atacseq"]["paired_exclude_flags"],
        mitochondrial_accession=config["atacseq"]["mitochondrial_accession"]

    threads:
        4

    resources:
        mem_mb=8000

    conda:
        "../envs/atacseq.yaml"

    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.bam:q})

        NUCLEAR_REFS=$(
            samtools idxstats {input.bam:q} \
            | awk -v mt="{params.mitochondrial_accession}" \
                '$1 != mt && $1 != "*" {{print $1}}'
        )

        if [[ "{params.structure}" == "SINGLE" ]]; then
            samtools view \
              -b \
              -q {params.min_mapq} \
              -F {params.single_exclude_flags} \
              -o {output.bam:q} \
              {input.bam:q} \
              $NUCLEAR_REFS

        elif [[ "{params.structure}" == "PAIRED" ]]; then
            TMP_ROOT="${{TMPDIR:-/tmp}}"
            TMP_DIR=$(
                mktemp -d \
                  "$TMP_ROOT/atacseq_filter_{wildcards.run_accession}.XXXXXX"
            )
            trap 'rm -rf "$TMP_DIR"' EXIT

            samtools view \
              -u \
              -q {params.min_mapq} \
              -f {params.paired_require_flags} \
              -F {params.paired_exclude_flags} \
              {input.bam:q} \
              $NUCLEAR_REFS \
            | samtools sort \
                -n \
                -@ {threads} \
                -m 750M \
                -T "$TMP_DIR/queryname" \
                -o "$TMP_DIR/candidates.name.bam" \
                -

            samtools view "$TMP_DIR/candidates.name.bam" \
            | cut -f1 \
            | uniq -c \
            | awk '$1 == 2 {{print $2}}' \
            > "$TMP_DIR/complete_pair_names.txt"

            test -s "$TMP_DIR/complete_pair_names.txt"

            samtools view \
              -u \
              -N "$TMP_DIR/complete_pair_names.txt" \
              "$TMP_DIR/candidates.name.bam" \
            | samtools sort \
                -@ {threads} \
                -m 750M \
                -T "$TMP_DIR/coordinate" \
                -o {output.bam:q} \
                -

            RECORD_COUNT=$(samtools view -c {output.bam:q})

            if (( RECORD_COUNT == 0 || RECORD_COUNT % 2 != 0 )); then
                printf 'Invalid paired ATAC-seq BAM record count: %s\n' \
                  "$RECORD_COUNT" >&2
                exit 1
            fi

        else
            printf 'Unsupported ATAC-seq structure: %s\n' \
              "{params.structure}" >&2
            exit 1
        fi

        samtools quickcheck -v {output.bam:q}
        samtools index {output.bam:q} {output.bai:q}
        """

rule atacseq_tss_reference:
    input:
        script="workflow/scripts/build_atac_tss_bed.py",
        gtf=config["reference"]["gtf"]

    output:
        bed=config["atacseq"]["tss_bed"],
        summary=config["atacseq"]["tss_summary"]

    conda:
        "../envs/reference.yaml"

    shell:
        """
        python {input.script} \
          --gtf {input.gtf:q} \
          --bed {output.bed:q} \
          --summary {output.summary:q}
        """

rule atacseq_tss_enrichment:
    input:
        bam=config["atacseq"]["filtered_bam_dir"] + "/{run_accession}/filtered.bam",
        bai=config["atacseq"]["filtered_bam_dir"] + "/{run_accession}/filtered.bam.bai",
        tss=config["atacseq"]["tss_bed"],
        script="workflow/scripts/calculate_atac_tss_enrichment.py"

    output:
        profile=config["atacseq"]["tss_profile"],
        summary=config["atacseq"]["tss_enrichment_summary"]

    params:
        window=config["atacseq"]["tss_window_bp"],
        bin_size=config["atacseq"]["tss_bin_size_bp"],
        background=config["atacseq"]["tss_background_bp"],
        forward_shift=config["atacseq"]["tn5_forward_shift"],
        reverse_shift=config["atacseq"]["tn5_reverse_shift"],
        min_mapq=config["atacseq"]["min_mapq"]

    conda:
        "../envs/atacseq_qc.yaml"

    shell:
        """
        python {input.script} \
          --bam {input.bam:q} \
          --tss-bed {input.tss:q} \
          --profile {output.profile:q} \
          --summary {output.summary:q} \
          --window {params.window} \
          --bin-size {params.bin_size} \
          --background {params.background} \
          --forward-shift {params.forward_shift} \
          --reverse-shift {params.reverse_shift} \
          --min-mapq {params.min_mapq}
        """
