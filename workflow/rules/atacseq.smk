"""Production ATAC-seq DAG. ATAC remains opt-in and is not in rule all."""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path("workflow/scripts").resolve()))
from atac_manifest import eligible_run

ATAC = config["atacseq"]
ATAC_MANIFEST = ATAC["production_manifest"]

def _manifest(_wildcards=None):
    return checkpoints.atacseq_production_manifest.get().output.manifest

def _rows(run, single=False):
    return eligible_run(str(_manifest()), run, require_single=single)

# Compatibility helper for the small, explicitly validation-only preprocessing
# subset. Production rules below never call it or validation_fastq_manifest.
def load_atacseq_validation_manifest():
    manifest = _load_validation_manifest_for_preprocessing()
    subset = manifest[manifest["omics"] == "ATAC-seq"].copy()
    return subset[subset["run_accession"].map(
        lambda accession: selection_decision(accession, omics="ATAC-seq")[0]
    )].copy()

def _structure(wildcards):
    return _rows(wildcards.run_accession)[0]["fastq_structure"]

def _gated(pattern, single=False):
    def resolve(wildcards):
        _rows(wildcards.run_accession, single)
        return pattern.format(run_accession=wildcards.run_accession)
    return resolve

def _raw_manifest(wildcards):
    matches=[r for r in _rows(wildcards.run_accession) if r["filename"]==wildcards.filename]
    if len(matches)!=1: raise ValueError("FASTQ is not uniquely production eligible")
    return str(_manifest())

def _preprocessed(run, role):
    rows=_rows(run); match=[r for r in rows if r["fastq_role"]==role]
    if len(match)!=1: raise ValueError(f"Expected one {role} FASTQ for {run}")
    base=match[0]["filename"].removesuffix(".fastq.gz")
    return f"{ATAC['preprocessed_dir']}/{run}/{match[0]['fastq_structure']}/{base}.fastq.gz"

def _alignment_reads(wildcards):
    roles=["SINGLE"] if _structure(wildcards)=="SINGLE" else ["R1","R2"]
    return [_preprocessed(wildcards.run_accession,r) for r in roles]

def _read_args(wildcards):
    run=wildcards.run_accession
    if _structure(wildcards)=="SINGLE": return "-U "+_preprocessed(run,"SINGLE")
    return f"-1 {_preprocessed(run,'R1')} -2 {_preprocessed(run,'R2')} -X {ATAC['paired_max_fragment_length']}"

with open(ATAC["replicate_groups"],newline="",encoding="utf-8") as _handle:
    _group_rows=list(csv.DictReader(_handle,delimiter="\t"))
ATAC_GROUP_MEMBERS={}
for _row in _group_rows:
    ATAC_GROUP_MEMBERS.setdefault(_row["replicate_group"],[]).append(_row["run_accession"])
ATAC_GROUPS=sorted(ATAC_GROUP_MEMBERS)

checkpoint atacseq_production_manifest:
    input:
        full_manifest=config["fastq"]["manifest"],
        samples="config/samples.tsv",
        selection=list(SELECTION_FILES),
        script="workflow/scripts/build_atac_production_manifest.py"
    output:
        manifest=ATAC_MANIFEST,
        provenance=ATAC["production_manifest_provenance"]
    conda: "../envs/fastq_io.yaml"
    shell:
        "python {input.script:q} --manifest {input.full_manifest:q} --output {output.manifest:q} --provenance {output.provenance:q}"

rule atacseq_reference:
    input: script="workflow/scripts/build_atac_reference.py", nuclear_fasta=config["reference"]["fasta"]
    output:
        mitochondrial_fasta=ATAC["mitochondrial_fasta"], mapping_fasta=ATAC["mapping_fasta"],
        metadata=ATAC["reference_metadata"], sha256=ATAC["reference_sha256"]
    params: mitochondrial=ATAC["mitochondrial_accession"]
    conda: "../envs/reference.yaml"
    shell:
        "python {input.script:q} --nuclear-fasta {input.nuclear_fasta:q} --mitochondrial-accession {params.mitochondrial:q} --mitochondrial-fasta {output.mitochondrial_fasta:q} --mapping-fasta {output.mapping_fasta:q} --metadata {output.metadata:q} --sha256 {output.sha256:q}"

rule atacseq_bowtie2_index:
    input: fasta=ATAC["mapping_fasta"]
    output: index_dir=directory(ATAC["bowtie2_index_dir"])
    params: prefix=ATAC["bowtie2_index"]
    threads: 8
    conda: "../envs/atacseq.yaml"
    shell: "mkdir -p {output.index_dir:q} && bowtie2-build --threads {threads} {input.fasta:q} {params.prefix:q}"

rule atacseq_reference_span:
    input: fasta=ATAC["mapping_fasta"], script="workflow/scripts/atac_reference_span.py"
    output: fai=ATAC["mapping_fasta"]+".fai", sizes=ATAC["nuclear_chrom_sizes"], metadata=ATAC["nuclear_span_metadata"]
    params: mitochondrial=ATAC["mitochondrial_accession"]
    conda: "../envs/atacseq.yaml"
    shell:
        "samtools faidx {input.fasta:q} && python {input.script:q} --fai {output.fai:q} --mitochondrial {params.mitochondrial:q} --sizes {output.sizes:q} --metadata {output.metadata:q}"

rule atacseq_download_fastq:
    input: manifest=_raw_manifest, script="workflow/scripts/download_fastq.py"
    output: fastq=ATAC["raw_fastq_dir"]+"/{run_accession}/{filename}"
    conda: "../envs/fastq_io.yaml"
    shell:
        "python {input.script:q} --manifest {input.manifest:q} --run-accession {wildcards.run_accession:q} --filename {wildcards.filename:q} --output {output.fastq:q}"

rule atacseq_fastp_single:
    input:
        manifest=_manifest,
        fastq=lambda w: ATAC["raw_fastq_dir"]+f"/{w.run_accession}/"+_rows(w.run_accession)[0]["filename"]
    output:
        fastq=ATAC["preprocessed_dir"]+"/{run_accession}/SINGLE/{run_accession}.fastq.gz",
        html=ATAC["preprocessed_dir"]+"/{run_accession}/SINGLE/fastp.html",
        json=ATAC["preprocessed_dir"]+"/{run_accession}/SINGLE/fastp.json"
    params: structure=_structure
    threads: 4
    conda: "../envs/preprocessing.yaml"
    shell:
        "test {params.structure:q} = SINGLE && mkdir -p $(dirname {output.fastq:q}) && fastp --thread {threads} --in1 {input.fastq:q} --out1 {output.fastq:q} --html {output.html:q} --json {output.json:q}"

rule atacseq_fastp_paired:
    input:
        manifest=_manifest,
        r1=lambda w: ATAC["raw_fastq_dir"]+f"/{w.run_accession}/{w.run_accession}_1.fastq.gz",
        r2=lambda w: ATAC["raw_fastq_dir"]+f"/{w.run_accession}/{w.run_accession}_2.fastq.gz"
    output:
        r1=ATAC["preprocessed_dir"]+"/{run_accession}/PAIRED/{run_accession}_1.fastq.gz",
        r2=ATAC["preprocessed_dir"]+"/{run_accession}/PAIRED/{run_accession}_2.fastq.gz",
        html=ATAC["preprocessed_dir"]+"/{run_accession}/PAIRED/fastp.html",
        json=ATAC["preprocessed_dir"]+"/{run_accession}/PAIRED/fastp.json"
    params: structure=_structure
    threads: 4
    conda: "../envs/preprocessing.yaml"
    shell:
        "test {params.structure:q} = PAIRED && mkdir -p $(dirname {output.r1:q}) && fastp --thread {threads} --in1 {input.r1:q} --in2 {input.r2:q} --out1 {output.r1:q} --out2 {output.r2:q} --html {output.html:q} --json {output.json:q}"

rule atacseq_align_run:
    input: manifest=_manifest, index_dir=ATAC["bowtie2_index_dir"], reads=_alignment_reads
    output:
        bam=ATAC["align_dir"]+"/{run_accession}/raw.sorted.bam",
        bai=ATAC["align_dir"]+"/{run_accession}/raw.sorted.bam.bai",
        log=ATAC["align_dir"]+"/{run_accession}/bowtie2.log"
    params: index=ATAC["bowtie2_index"], args=_read_args, mode="--"+ATAC["bowtie2_mode"], preset="--"+ATAC["bowtie2_preset"]
    threads: 8
    resources: mem_mb=8000
    conda: "../envs/atacseq.yaml"
    shell:
        "mkdir -p $(dirname {output.bam:q}) && bowtie2 {params.mode} {params.preset} --threads {threads} --rg-id {wildcards.run_accession:q} --rg SM:{wildcards.run_accession:q} -x {params.index:q} {params.args} 2> {output.log:q} | samtools sort -o {output.bam:q} - && samtools index {output.bam:q} {output.bai:q}"

rule atacseq_mark_duplicates:
    input: manifest=_manifest, bam=_gated(ATAC["align_dir"]+"/{run_accession}/raw.sorted.bam")
    output:
        bam=ATAC["align_dir"]+"/{run_accession}/dupmarked.sorted.bam",
        bai=ATAC["align_dir"]+"/{run_accession}/dupmarked.sorted.bai",
        metrics=ATAC["qc_dir"]+"/{run_accession}/duplication_metrics.txt"
    resources: mem_mb=4000
    conda: "../envs/atacseq_picard.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bam:q}) $(dirname {output.metrics:q})
        picard_tmp=$(mktemp -d {resources.tmpdir:q}/picard_markduplicates.XXXXXX)
        trap 'rm -rf -- "$picard_tmp"' EXIT
        JAVA_TOOL_OPTIONS="${{JAVA_TOOL_OPTIONS:-}} -Djava.io.tmpdir=$picard_tmp" picard -Xmx3g MarkDuplicates I={input.bam:q} O={output.bam:q} M={output.metrics:q} TMP_DIR="$picard_tmp" REMOVE_DUPLICATES=false CREATE_INDEX=true
        test -s {output.bai:q}
        """

rule atacseq_filter_bam:
    input: manifest=_manifest, bam=_gated(ATAC["align_dir"]+"/{run_accession}/dupmarked.sorted.bam"), script="workflow/scripts/filter_atac_bam.py"
    output:
        bam=ATAC["filtered_bam_dir"]+"/{run_accession}/filtered.bam",
        bai=ATAC["filtered_bam_dir"]+"/{run_accession}/filtered.bam.bai",
        metrics=ATAC["filtered_bam_dir"]+"/{run_accession}/filter_metrics.json"
    params:
        layout=_structure, min_mapq=ATAC["min_mapq"], single=ATAC["single_exclude_flags"],
        preq=ATAC["paired_require_flags"], pexc=ATAC["paired_exclude_flags"], mitochondrial=ATAC["mitochondrial_accession"]
    resources: mem_mb=8000
    conda: "../envs/atacseq_qc.yaml"
    shell:
        "python {input.script:q} --input {input.bam:q} --output {output.bam:q} --metrics {output.metrics:q} --layout {params.layout:q} --min-mapq {params.min_mapq} --single-exclude {params.single} --paired-require {params.preq} --paired-exclude {params.pexc} --mitochondrial {params.mitochondrial:q} --tmpdir {resources.tmpdir:q}"

rule atacseq_tss_reference:
    input: script="workflow/scripts/build_atac_tss_bed.py", gtf=config["reference"]["gtf"]
    output: bed=ATAC["tss_bed"], summary=ATAC["tss_summary"]
    conda: "../envs/reference.yaml"
    shell: "python {input.script:q} --gtf {input.gtf:q} --bed {output.bed:q} --summary {output.summary:q}"

rule atacseq_tss_enrichment:
    input:
        manifest=_manifest, bam=_gated(ATAC["filtered_bam_dir"]+"/{run_accession}/filtered.bam"),
        bai=ATAC["filtered_bam_dir"]+"/{run_accession}/filtered.bam.bai", tss=ATAC["tss_bed"],
        script="workflow/scripts/calculate_atac_tss_enrichment.py"
    output: profile=ATAC["tss_profile"], summary=ATAC["tss_enrichment_summary"]
    params:
        window=ATAC["tss_window_bp"], bin=ATAC["tss_bin_size_bp"], background=ATAC["tss_background_bp"],
        forward_shift=ATAC["tn5_forward_shift"], reverse_shift=ATAC["tn5_reverse_shift"], mapq=ATAC["min_mapq"]
    conda: "../envs/atacseq_qc.yaml"
    shell:
        "python {input.script:q} --bam {input.bam:q} --tss-bed {input.tss:q} --profile {output.profile:q} --summary {output.summary:q} --window {params.window} --bin-size {params.bin} --background {params.background} --forward-shift {params.forward_shift} --reverse-shift {params.reverse_shift} --min-mapq {params.mapq}"

rule atacseq_tn5_transform:
    input: manifest=_manifest, bam=_gated(ATAC["filtered_bam_dir"]+"/{run_accession}/filtered.bam",True), script="workflow/scripts/tn5_transform.py"
    output:
        insertions=ATAC["tn5_dir"]+"/{run_accession}/insertions.bed",
        tagalign=ATAC["tn5_dir"]+"/{run_accession}/shifted_read_intervals.bed",
        metrics=ATAC["tn5_dir"]+"/{run_accession}/metrics.json"
    params: forward_shift=ATAC["tn5_forward_shift"], reverse_shift=ATAC["tn5_reverse_shift"]
    resources: mem_mb=2000
    conda: "../envs/atacseq_qc.yaml"
    shell: "python {input.script:q} --bam {input.bam:q} --insertions {output.insertions:q} --tagalign {output.tagalign:q} --metrics {output.metrics:q} --forward-shift {params.forward_shift} --reverse-shift {params.reverse_shift} --tmpdir {resources.tmpdir:q}"

if ATAC["macs3_smooth_window_bp"]<=0 or ATAC["macs3_smooth_window_bp"]%2:
    raise ValueError("macs3_smooth_window_bp must be positive and even")
if ATAC["macs3_keep_dup"] != "all":
    raise ValueError("macs3_keep_dup must remain 'all': duplicates are removed upstream")
MACS_SHIFT=-(ATAC["macs3_smooth_window_bp"]//2)

rule atacseq_macs3_peaks:
    input: manifest=_manifest, tagalign=_gated(ATAC["tn5_dir"]+"/{run_accession}/shifted_read_intervals.bed",True), span=ATAC["nuclear_span_metadata"]
    output:
        peaks=ATAC["peaks_dir"]+"/{run_accession}/{run_accession}_peaks.narrowPeak",
        pileup=ATAC["peaks_dir"]+"/{run_accession}/{run_accession}_treat_pileup.bdg",
        summits=ATAC["peaks_dir"]+"/{run_accession}/{run_accession}_summits.bed",
        xls=ATAC["peaks_dir"]+"/{run_accession}/{run_accession}_peaks.xls"
    params:
        outdir=lambda w:ATAC["peaks_dir"]+"/"+w.run_accession, extsize=ATAC["macs3_smooth_window_bp"],
        shift=MACS_SHIFT, pvalue=ATAC["macs3_pvalue"], keepdup=ATAC["macs3_keep_dup"]
    conda: "../envs/atacseq_downstream.yaml"
    shell:
        "mkdir -p {params.outdir:q} && genome_size=$(python -c 'import json,sys;print(json.load(open(sys.argv[1]))[\"nuclear_reference_span\"])' {input.span:q}) && macs3 callpeak -t {input.tagalign:q} -f BED -g \"$genome_size\" -n {wildcards.run_accession:q} --outdir {params.outdir:q} --nomodel --shift {params.shift} --extsize {params.extsize} --keep-dup {params.keepdup:q} --call-summits -B --SPMR -p {params.pvalue}"

rule atacseq_frip:
    input:
        manifest=_manifest, tagalign=_gated(ATAC["tn5_dir"]+"/{run_accession}/shifted_read_intervals.bed",True),
        peaks=ATAC["peaks_dir"]+"/{run_accession}/{run_accession}_peaks.narrowPeak", script="workflow/scripts/calculate_atac_frip.py"
    output: metrics=ATAC["qc_dir"]+"/{run_accession}/frip.json"
    conda: "../envs/atacseq_downstream.yaml"
    shell: "python {input.script:q} --tagalign {input.tagalign:q} --peaks {input.peaks:q} --run {wildcards.run_accession:q} --output {output.metrics:q}"

rule atacseq_bigwig:
    input:
        manifest=_manifest, pileup=_gated(ATAC["peaks_dir"]+"/{run_accession}/{run_accession}_treat_pileup.bdg",True),
        sizes=ATAC["nuclear_chrom_sizes"], clip="workflow/scripts/clip_atac_bedgraph.py", validate="workflow/scripts/validate_atac_bigwig.py"
    output:
        clipped=ATAC["tracks_dir"]+"/{run_accession}/treatment_pileup.clipped.bedGraph",
        metrics=ATAC["tracks_dir"]+"/{run_accession}/clipping_metrics.json",
        bigwig=ATAC["tracks_dir"]+"/{run_accession}/treatment_pileup.SPMR.bw",
        sha256=ATAC["tracks_dir"]+"/{run_accession}/treatment_pileup.SPMR.bw.sha256"
    conda: "../envs/atacseq_tracks.yaml"
    shell:
        "python {input.clip:q} --input {input.pileup:q} --chrom-sizes {input.sizes:q} --output {output.clipped:q} --metrics {output.metrics:q} && bedGraphToBigWig {output.clipped:q} {input.sizes:q} {output.bigwig:q} && python {input.validate:q} --bigwig {output.bigwig:q} --chrom-sizes {input.sizes:q} --clipped-bedgraph {output.clipped:q} --sha256 {output.sha256:q}"

rule atacseq_run_provenance:
    input:
        manifest=_manifest, reference=ATAC["reference_metadata"], reference_hashes=ATAC["reference_sha256"], span=ATAC["nuclear_span_metadata"],
        filtered=_gated(ATAC["filtered_bam_dir"]+"/{run_accession}/filtered.bam",True), tn5=ATAC["tn5_dir"]+"/{run_accession}/metrics.json",
        filter_metrics=ATAC["filtered_bam_dir"]+"/{run_accession}/filter_metrics.json",
        tss_profile=ATAC["tss_profile"], tss_summary=ATAC["tss_enrichment_summary"],
        peaks=ATAC["peaks_dir"]+"/{run_accession}/{run_accession}_peaks.narrowPeak", frip=ATAC["qc_dir"]+"/{run_accession}/frip.json",
        bigwig=ATAC["tracks_dir"]+"/{run_accession}/treatment_pileup.SPMR.bw", script="workflow/scripts/atac_provenance.py"
    output: json=ATAC["provenance_dir"]+"/{run_accession}/provenance.json"
    params:
        values=json.dumps({"min_mapq":ATAC["min_mapq"],"single_exclude_flags":ATAC["single_exclude_flags"],
          "tn5":{"forward":ATAC["tn5_forward_shift"],"reverse":ATAC["tn5_reverse_shift"]},
          "macs3":{"format":"BED","nomodel":True,"smooth_window_bp":ATAC["macs3_smooth_window_bp"],"shift":MACS_SHIFT,
          "extsize":ATAC["macs3_smooth_window_bp"],"keep_dup":ATAC["macs3_keep_dup"],"pvalue":ATAC["macs3_pvalue"]},
          "software":{"bowtie2":"2.5.5","samtools":"1.24","picard":"3.5.0","pysam":"0.24.0","numpy":"2.5.2","macs3":"3.0.4","idr":"2.0.4.2"}})
    shell:
        "python {input.script:q} --run {wildcards.run_accession:q} --manifest {input.manifest:q} --reference-metadata {input.reference:q} --input {input.reference_hashes:q} --input {input.span:q} --input {input.filtered:q} --input {input.filter_metrics:q} --input {input.tn5:q} --input {input.tss_profile:q} --input {input.tss_summary:q} --input {input.peaks:q} --input {input.frip:q} --input {input.bigwig:q} --parameters-json {params.values:q} --output {output.json:q}"

rule atacseq_validate_run:
    input:
        manifest=_manifest, bam=_gated(ATAC["filtered_bam_dir"]+"/{run_accession}/filtered.bam",True),
        bai=ATAC["filtered_bam_dir"]+"/{run_accession}/filtered.bam.bai", profile=ATAC["tss_profile"],
        summary=ATAC["tss_enrichment_summary"], peaks=ATAC["peaks_dir"]+"/{run_accession}/{run_accession}_peaks.narrowPeak",
        frip=ATAC["qc_dir"]+"/{run_accession}/frip.json", bigwig=ATAC["tracks_dir"]+"/{run_accession}/treatment_pileup.SPMR.bw",
        script="workflow/scripts/validate_atac_run.py"
    output: report=ATAC["qc_dir"]+"/{run_accession}/production_validation.json"
    params: mapq=ATAC["min_mapq"], mitochondrial=ATAC["mitochondrial_accession"], window=ATAC["tss_window_bp"], bin=ATAC["tss_bin_size_bp"], background=ATAC["tss_background_bp"]
    conda: "../envs/atacseq_qc.yaml"
    shell:
        "python {input.script:q} --run {wildcards.run_accession:q} --manifest {input.manifest:q} --bam {input.bam:q} --bai {input.bai:q} --profile {input.profile:q} --summary {input.summary:q} --peaks {input.peaks:q} --frip {input.frip:q} --bigwig {input.bigwig:q} --min-mapq {params.mapq} --mitochondrial {params.mitochondrial:q} --window {params.window} --bin-size {params.bin} --background {params.background} --output {output.report:q}"

rule atacseq_validate_replicate_group:
    input: manifest=_manifest, groups=ATAC["replicate_groups"], script="workflow/scripts/validate_atac_replicates.py"
    output: marker=ATAC["reproducibility_dir"]+"/{replicate_group}/eligibility.json"
    shell: "python {input.script:q} --groups {input.groups:q} --manifest {input.manifest:q} --group {wildcards.replicate_group:q} --output {output.marker:q}"

def _group_tags(wildcards):
    runs=ATAC_GROUP_MEMBERS.get(wildcards.replicate_group,[])
    if not runs: raise ValueError("Unknown ATAC replicate group")
    return [ATAC["tn5_dir"]+f"/{r}/shifted_read_intervals.bed" for r in runs]

def _group_peaks(wildcards):
    runs=ATAC_GROUP_MEMBERS.get(wildcards.replicate_group,[])
    if len(runs)!=2: raise ValueError("Current IDR integration requires exactly two members")
    return [ATAC["peaks_dir"]+f"/{r}/{r}_peaks.narrowPeak" for r in runs]

rule atacseq_pool_tagalign:
    input: eligible=ATAC["reproducibility_dir"]+"/{replicate_group}/eligibility.json", tags=_group_tags
    output: pooled=ATAC["reproducibility_dir"]+"/{replicate_group}/pooled.shifted_read_intervals.bed"
    resources: mem_mb=2000
    shell: "LC_ALL=C sort --stable --temporary-directory {resources.tmpdir:q} -k1,1 -k2,2n -k3,3n -k4,4 -k6,6 -k5,5n {input.tags:q} > {output.pooled:q}"

rule atacseq_pooled_macs3:
    input: eligible=ATAC["reproducibility_dir"]+"/{replicate_group}/eligibility.json", tagalign=ATAC["reproducibility_dir"]+"/{replicate_group}/pooled.shifted_read_intervals.bed", span=ATAC["nuclear_span_metadata"]
    output: peaks=ATAC["reproducibility_dir"]+"/{replicate_group}/pooled_peaks.narrowPeak"
    params: outdir=lambda w:ATAC["reproducibility_dir"]+"/"+w.replicate_group, shift=MACS_SHIFT, extsize=ATAC["macs3_smooth_window_bp"], pvalue=ATAC["macs3_pvalue"], keepdup=ATAC["macs3_keep_dup"]
    conda: "../envs/atacseq_downstream.yaml"
    shell:
        "genome_size=$(python -c 'import json,sys;print(json.load(open(sys.argv[1]))[\"nuclear_reference_span\"])' {input.span:q}) && macs3 callpeak -t {input.tagalign:q} -f BED -g \"$genome_size\" -n pooled --outdir {params.outdir:q} --nomodel --shift {params.shift} --extsize {params.extsize} --keep-dup {params.keepdup:q} --call-summits -B --SPMR -p {params.pvalue}"

rule atacseq_idr:
    input: eligible=ATAC["reproducibility_dir"]+"/{replicate_group}/eligibility.json", peaks=_group_peaks, pooled=ATAC["reproducibility_dir"]+"/{replicate_group}/pooled_peaks.narrowPeak"
    output: ranked=ATAC["reproducibility_dir"]+"/{replicate_group}/idr.all.tsv"
    params: seed=ATAC["idr_seed"]
    conda: "../envs/atacseq_idr.yaml"
    shell: "idr --samples {input.peaks:q} --peak-list {input.pooled:q} --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.10 --random-seed {params.seed} --output-file {output.ranked:q}"

rule atacseq_idr_thresholds:
    input: ranked=ATAC["reproducibility_dir"]+"/{replicate_group}/idr.all.tsv", script="workflow/scripts/filter_idr.py"
    output:
        idr010=ATAC["reproducibility_dir"]+"/{replicate_group}/idr_le_0.10.tsv",
        idr005=ATAC["reproducibility_dir"]+"/{replicate_group}/idr_le_0.05.tsv",
        metrics=ATAC["reproducibility_dir"]+"/{replicate_group}/reproducibility_metrics.json"
    shell: "python {input.script:q} --input {input.ranked:q} --idr-010 {output.idr010:q} --idr-005 {output.idr005:q} --metrics {output.metrics:q}"

rule atacseq_reproducibility_provenance:
    input:
        eligibility=ATAC["reproducibility_dir"]+"/{replicate_group}/eligibility.json",
        span=ATAC["nuclear_span_metadata"],
        pooled=ATAC["reproducibility_dir"]+"/{replicate_group}/pooled.shifted_read_intervals.bed",
        pooled_peaks=ATAC["reproducibility_dir"]+"/{replicate_group}/pooled_peaks.narrowPeak",
        replicate_peaks=_group_peaks,
        ranked=ATAC["reproducibility_dir"]+"/{replicate_group}/idr.all.tsv",
        idr010=ATAC["reproducibility_dir"]+"/{replicate_group}/idr_le_0.10.tsv",
        idr005=ATAC["reproducibility_dir"]+"/{replicate_group}/idr_le_0.05.tsv",
        metrics=ATAC["reproducibility_dir"]+"/{replicate_group}/reproducibility_metrics.json",
        script="workflow/scripts/atac_reproducibility_provenance.py"
    output: provenance=ATAC["reproducibility_dir"]+"/{replicate_group}/provenance.json"
    params:
        values=json.dumps({"macs3":{"smooth_window_bp":ATAC["macs3_smooth_window_bp"],"shift":MACS_SHIFT,
          "extsize":ATAC["macs3_smooth_window_bp"],"keep_dup":ATAC["macs3_keep_dup"],"pvalue":ATAC["macs3_pvalue"]},
          "idr":{"version":"2.0.4.2","seed":ATAC["idr_seed"],"thresholds":[0.10,0.05]}})
    shell:
        "python {input.script:q} --group {wildcards.replicate_group:q} --eligibility {input.eligibility:q} --input {input.span:q} --input {input.pooled:q} --input {input.pooled_peaks:q} --input {input.replicate_peaks:q} --input {input.ranked:q} --input {input.idr010:q} --input {input.idr005:q} --input {input.metrics:q} --parameters-json {params.values:q} --output {output.provenance:q}"

def _production_targets(_wildcards):
    from atac_manifest import eligible_runs
    targets=[]
    for run in eligible_runs(str(_manifest()),require_single=True):
        targets += [ATAC["tss_enrichment_summary"].format(run_accession=run), ATAC["qc_dir"]+f"/{run}/frip.json",
          ATAC["tracks_dir"]+f"/{run}/treatment_pileup.SPMR.bw.sha256", ATAC["provenance_dir"]+f"/{run}/provenance.json",
          ATAC["qc_dir"]+f"/{run}/production_validation.json"]
    for group in ATAC_GROUPS:
        targets += [ATAC["reproducibility_dir"]+f"/{group}/idr_le_0.10.tsv", ATAC["reproducibility_dir"]+f"/{group}/idr_le_0.05.tsv",
          ATAC["reproducibility_dir"]+f"/{group}/reproducibility_metrics.json", ATAC["reproducibility_dir"]+f"/{group}/provenance.json"]
    return targets

rule atacseq_production:
    input: _production_targets

rule atacseq_reproducibility:
    input:
        expand(ATAC["reproducibility_dir"]+"/{replicate_group}/idr_le_0.10.tsv",replicate_group=ATAC_GROUPS),
        expand(ATAC["reproducibility_dir"]+"/{replicate_group}/idr_le_0.05.tsv",replicate_group=ATAC_GROUPS),
        expand(ATAC["reproducibility_dir"]+"/{replicate_group}/provenance.json",replicate_group=ATAC_GROUPS)
