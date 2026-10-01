#!/usr/bin/env python3
"""Portable adapters around the frozen, validated ChIP-seq scientific helpers.

No execution occurs on import. Historical evidence is read-only. All generated
runtime manifests and intermediates belong to an explicit reproduction root.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EVIDENCE = ROOT / "provenance/chipseq/production_2026-09-23"


def module(name):
    path = ROOT / "workflow/scripts" / (name + ".py")
    spec = importlib.util.spec_from_file_location(name + "_portable", path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


PRODUCTION = module("chipseq_production")
DYNAMIC = PRODUCTION.DYNAMIC
BENCHMARK = PRODUCTION.BENCHMARK


def sha(path):
    return PRODUCTION.sha256(path)


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    PRODUCTION.atomic_json(path, value)


def rows(path):
    return PRODUCTION.table(path)


def location(value):
    p = Path(value)
    return p.resolve() if p.is_absolute() else (ROOT / p).resolve()


def configuration(path="config/chipseq.yaml", overrides=None):
    import yaml
    value = yaml.safe_load(location(path).read_text())["chipseq"]
    value.update(overrides or {})
    if value["mode"] not in {"frozen", "live_metadata"}:
        raise ValueError("Unsupported metadata mode")
    expected = {
        "samples": "config/samples.tsv", "analysis_plan": str(PRODUCTION.ANALYSIS_PLAN.relative_to(ROOT)),
        "processing_plan": str(PRODUCTION.PROCESSING_PLAN.relative_to(ROOT)),
        "peak_plan": str(PRODUCTION.PEAK_PLAN.relative_to(ROOT)),
        "peak_plan_summary": str(PRODUCTION.PEAK_SUMMARY.relative_to(ROOT)),
        "peak_policy": str(PRODUCTION.PEAK_POLICY.relative_to(ROOT)),
        "alignment_policy": "config/chipseq_alignment.json", "filtering_policy": "config/chipseq_filtering.json",
    }
    if any(value.get(k) != v for k, v in expected.items()):
        raise ValueError("Frozen scientific input locations cannot be redirected")
    for key in ("output_root", "reference_root"):
        p = location(value[key])
        if p == ROOT or ROOT in p.parents and p.relative_to(ROOT).parts[0] not in {"results", "resources", "data"}:
            raise ValueError("Generated output may not overwrite source or frozen evidence")
    if type(value["max_parallel_analyses"]) is not int or value["max_parallel_analyses"] not in (1, 2):
        raise ValueError("Concurrency must be one or two analyses")
    return value


@contextmanager
def production_slot(settings):
    """Limit concurrent workers independently of the scheduler or site."""
    import fcntl
    root = location(settings["output_root"])
    root.mkdir(parents=True, exist_ok=True)
    handle = None
    try:
        while handle is None:
            for number in range(settings["max_parallel_analyses"]):
                candidate = (root / f".portable.slot.{number}.lock").open("a+")
                try:
                    fcntl.flock(candidate, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    candidate.close()
                else:
                    handle = candidate
                    break
            if handle is None:
                time.sleep(0.2)
        yield
    finally:
        if handle is not None:
            fcntl.flock(handle, fcntl.LOCK_UN)
            handle.close()


def check_manifest(path, root=ROOT):
    for line in Path(path).read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        p = Path(name.lstrip("*"))
        if p.is_absolute() or ".." in p.parts:
            raise ValueError("Unsafe checksum manifest member")
        DYNAMIC.verify_hash(Path(root) / p, digest, name)


def validate_contract():
    check_manifest(EVIDENCE / "frozen_inputs.sha256")
    peaks, ready = PRODUCTION.validated_plans()
    contract = read_json(EVIDENCE / "frozen_contract.json")
    ids = {r["analysis_id"] for r in peaks}
    if ids != set(contract["analysis_parameters"]):
        raise ValueError("Historical analysis identity mismatch")
    physical = {r[k] for r in peaks for k in ("ip_run_accession", "control_run_accession")}
    controls = {r["control_run_accession"] for r in peaks}
    if (len(ids), len(physical), len(controls)) != (18, 22, 4):
        raise ValueError("Frozen cohort size mismatch")
    for r in peaks:
        selected = PRODUCTION.select_analysis(r["analysis_id"], r["control_run_accession"])
        if selected["analysis"] != r:
            raise ValueError("Frozen analysis selection mismatch")
        old = contract["analysis_parameters"][r["analysis_id"]]
        if any(old[k] != r[k] for k in ("analysis_id", "ip_run_accession", "control_run_accession", "peak_mode")):
            raise ValueError("Historical pairing/mode mismatch")
        if r["study_accession"] == "PRJNA865478":
            size, _ = DYNAMIC.resolve_fragment_size(r)
            if size != 147:
                raise ValueError("Fixed-fragment mismatch")
        else:
            size = old["fragment_size_bp"]
            if r["fragment_size_policy"] != "phantompeakqualtools" or r["fragment_size_bp"]:
                raise ValueError("External-fragment policy mismatch")
        if float(r["qvalue"]) != old["qvalue"] or r["scale_to"] != old["scale_to"] or r["keep_dup"] != old["keep_dup"]:
            raise ValueError("Historical MACS3 policy mismatch")
        expected_flags = ["--nomodel", "--extsize", str(size), "--keep-dup", "all", "-B", "--SPMR"]
        if r["peak_mode"] == "broad":
            expected_flags += ["--broad", "--broad-cutoff", str(r["broad_cutoff"])]
        if old["flags"] != expected_flags or r["spmr"] != "true" or r["store_bdg"] != "true":
            raise ValueError("Historical MACS3 flags mismatch")
        if r["fragment_size_policy"] == "fixed":
            if DYNAMIC.macs3_contract(r, 2366634374)["flags"] != old["flags"]:
                raise ValueError("Fixed MACS3 contract mismatch")
    return {"schema_version": 1, "scope": "frozen_plan_validation", "analyses": len(peaks),
            "physical_runs": len(physical), "shared_inputs": len(controls),
            "narrow": sum(r["peak_mode"] == "narrow" for r in peaks),
            "broad": sum(r["peak_mode"] == "broad" for r in peaks)}


def compare_plans(analysis, processing):
    for actual, frozen, key in ((analysis, PRODUCTION.ANALYSIS_PLAN, "analysis_id"),
                               (processing, PRODUCTION.PROCESSING_PLAN, "run_accession")):
        a, b = rows(actual), rows(frozen)
        if len(a) != len({r[key] for r in a}) or {r[key]: r for r in a} != {r[key]: r for r in b}:
            raise ValueError("Reconstructed plan differs from the authoritative historical plan")
    return validate_contract()


def metadata_record(plan, outdir, accession):
    validate_contract()
    module("chipseq_metadata_snapshot").materialize(plan, outdir, accession)


def fastq_manifest(output):
    validate_contract()
    peaks, _ = PRODUCTION.validated_plans()
    values = {}
    for r in peaks:
        for role in ("ip", "input"):
            run = r["ip_run_accession" if role == "ip" else "control_run_accession"]
            values[run] = PRODUCTION.sample(run, r["study_accession"])
    import io
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, ["run_accession", "filename", "fastq_ftp", "fastq_md5", "fastq_bytes"], delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(values[k] for k in sorted(values))
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    PRODUCTION.atomic_text(output, buffer.getvalue())


@contextmanager
def cwd(path):
    before = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(before)


def run(command, log=None, output=None, output_binary=None):
    if log:
        Path(log).parent.mkdir(parents=True, exist_ok=True)
    with open(log or os.devnull, "w") as err:
        if output_binary:
            with open(output_binary, "wb") as out:
                subprocess.run([str(x) for x in command], check=True, stdout=out, stderr=err)
        elif output:
            with open(output, "w") as out:
                subprocess.run([str(x) for x in command], check=True, stdout=out, stderr=err)
        else:
            subprocess.run([str(x) for x in command], check=True, stdout=err, stderr=subprocess.STDOUT)


def helper(name, *args):
    run([sys.executable, ROOT / "workflow/scripts" / (name + ".py"), *args])


def bowtie2_command(run_id, fastq, index):
    return ["bowtie2", "--local", "--very-sensitive-local", "--seed", "0", "-p", "6", "-x", str(index),
            "-U", str(fastq), "--rg-id", run_id, "--rg", "SM:" + run_id,
            "--rg", "LB:" + run_id, "--rg", "PL:ILLUMINA"]


def picard_command(bam, output, metrics, temporary):
    return ["picard", "-Xmx12000m", "-XX:ActiveProcessorCount=4", "MarkDuplicates",
            "INPUT=" + str(bam), "OUTPUT=" + str(output), "METRICS_FILE=" + str(metrics),
            "TMP_DIR=" + str(temporary), "REMOVE_DUPLICATES=false", "REMOVE_SEQUENCING_DUPLICATES=false",
            "READ_NAME_REGEX=null", "DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES",
            "CREATE_INDEX=false", "VALIDATION_STRINGENCY=STRICT", "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=256"]


def remove_output_scaffold(path, directories=("",)):
    """Claim only Snakemake's known empty parent directories; never delete data."""
    path = Path(path)
    if not path.exists() and not path.is_symlink():
        return
    members = [path, *path.rglob("*")]
    if any(p.is_symlink() or not p.is_dir() or str(p.relative_to(path)) not in
           {"." if d == "" else d for d in directories} for p in members):
        raise ValueError("Output destination collision: nonempty or unexpected scaffold")
    for p in sorted(members, key=lambda p: len(p.parts), reverse=True):
        p.rmdir()


def assert_versions():
    import re
    expected = {"bowtie2": "2.5.5", "samtools": "1.24", "macs3": "3.0.4", "bedtools": "2.31.1", "fastp": "1.3.6", "fastqc": "0.12.1"}
    records = {}
    for binary, version in expected.items():
        text = subprocess.check_output([binary, "--version"], stderr=subprocess.STDOUT, text=True)
        if not re.search(r"(?<![0-9.])" + re.escape(version) + r"(?![0-9.])", text):
            raise ValueError("Unexpected tool version: " + binary)
        records[binary] = text.strip()
    import pysam
    if pysam.__version__ != "0.24.0":
        raise ValueError("Unexpected pysam version")
    records["pysam"] = pysam.__version__
    # Package-manager provenance verifies Picard and PhantomPeakQualTools, which
    # do not both expose reliable --version interfaces.
    prefix = Path(sys.prefix)
    for name, version in (("picard", "3.5.0"), ("phantompeakqualtools", "1.2.2")):
        matches = list((prefix / "conda-meta").glob(name + "-*.json"))
        packages = [read_json(f) for f in matches]
        if len(packages) != 1 or packages[0].get("version") != version:
            raise ValueError("Unexpected or unverifiable Conda package: " + name)
        records[name] = packages[0]
    return records


def reference_build(dest):
    import yaml
    validate_contract()
    dest = Path(dest).resolve()
    remove_output_scaffold(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".chipseq-reference-", dir=dest.parent) as temp:
        project = Path(temp); (project / "config").mkdir()
        cfg = yaml.safe_load((ROOT / "config/config.yaml").read_text())
        ref = cfg["reference"]
        if any(str(ref.get(k)) != str(BENCHMARK.REFERENCE_IDENTITY[k]) for k in
               ("assembly_name", "refseq_accession", "genbank_accession", "species", "taxid")):
            raise ValueError("Reference acquisition identity differs from frozen production")
        # Keep all identity fields; relocate acquisition products into this job.
        ref["dir"] = "nuclear"
        for key, name in (("fasta", "genome.fa"), ("gff3", "annotation.gff3"), ("gtf", "annotation.gtf"),
                          ("sequence_report", "sequence_report.jsonl"), ("metadata", "reference_metadata.tsv"), ("sha256", "SHA256SUMS.txt")):
            ref[key] = "nuclear/" + name
        (project / "config/config.yaml").write_text(yaml.safe_dump(cfg))
        shutil.copyfile(ROOT / "config/chipseq_alignment.json", project / "config/chipseq_alignment.json")
        with cwd(project):
            helper("fetch_reference_genome")
            helper("chipseq_alignment_support", "reference")
        built = project / "nuclear/chipseq"
        for label, name in (("nuclear", project / "nuclear/genome.fa"), ("mitochondrial", built / "mitochondrial.fa"), ("mapping", built / "genome_plus_mt.fa")):
            DYNAMIC.verify_hash(name, BENCHMARK.REFERENCE_CONTENT_SHA256[label], label + " sequence")
        run(["samtools", "faidx", built / "genome_plus_mt.fa"])
        DYNAMIC.verify_hash(built / "genome_plus_mt.fa.fai", BENCHMARK.REFERENCE_SOURCE_SHA256["genome_plus_mt.fa.fai"], "reference FAI")
        (built / "bowtie2_index").mkdir()
        version = subprocess.check_output(["bowtie2-build", "--version"], text=True)
        if "2.5.5" not in version:
            raise ValueError("Reference requires Bowtie2 2.5.5")
        run(["bowtie2-build", "--large-index", "--threads", "8", built / "genome_plus_mt.fa", built / "bowtie2_index/genome_plus_mt"], built / "bowtie2_build.log")
        inventory = {"schema_version": 1, "scope": "portable_rebuilt_reference", "identity": BENCHMARK.REFERENCE_IDENTITY,
                     "source_content_sha256": BENCHMARK.REFERENCE_CONTENT_SHA256,
                     "files": [{"relative_path": str(f.relative_to(built)), "sha256": sha(f), "bytes": f.stat().st_size}
                               for f in sorted(built.rglob("*")) if f.is_file()]}
        write_json(built / "reference_inventory.json", inventory)
        reference_validate(built)
        built.replace(dest)


def reference_validate(root):
    root = Path(root).resolve()
    metadata = read_json(root / "reference_provenance.json")
    inventory = read_json(root / "reference_inventory.json")
    if inventory.get("identity") != BENCHMARK.REFERENCE_IDENTITY:
        raise ValueError("Reference inventory identity mismatch")
    for label, key in (("nuclear", "nuclear_sha256"), ("mitochondrial", "mitochondrial_sha256"), ("mapping", "mapping_sha256")):
        if metadata.get(key) != BENCHMARK.REFERENCE_CONTENT_SHA256[label]:
            raise ValueError("Reference sequence provenance mismatch")
    expected = {"genome_plus_mt.fa.fai", "reference_provenance.json", *BENCHMARK.INDEX_COMPONENTS}
    names = {r["relative_path"] for r in inventory["files"]}
    if not expected <= names:
        raise ValueError("Reference inventory missing required files")
    for item in inventory["files"]:
        name = Path(item["relative_path"])
        if name.is_absolute() or ".." in name.parts:
            raise ValueError("Unsafe reference inventory path")
        f = root / name
        if f.stat().st_size != item["bytes"]:
            raise ValueError("Reference size mismatch")
        DYNAMIC.verify_hash(f, item["sha256"], "reference component")
    fasta = root / "genome_plus_mt.fa"
    if fasta.exists():
        observed = sha(fasta)
    else:
        digest = hashlib.sha256()
        with gzip.open(root / "genome_plus_mt.fa.gz", "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        observed = digest.hexdigest()
    if observed != BENCHMARK.REFERENCE_CONTENT_SHA256["mapping"]:
        raise ValueError("Reference mapping sequence mismatch")
    DYNAMIC.verify_hash(root / "genome_plus_mt.fa.fai", BENCHMARK.REFERENCE_SOURCE_SHA256["genome_plus_mt.fa.fai"], "FAI")
    return DYNAMIC.validate_reference(root / "genome_plus_mt.fa.fai", metadata, sha(root / "genome_plus_mt.fa.fai"))


def execution_plan(analysis, reference):
    selected = PRODUCTION.select_analysis(analysis)
    reference_validate(reference)
    value = dict(selected)
    value["reference"] = {"root": str(Path(reference).resolve())}
    value["input_sha256"] = {name: sha(path) for name, path in {
        "analysis_plan": PRODUCTION.ANALYSIS_PLAN, "peak_plan": PRODUCTION.PEAK_PLAN,
        "peak_plan_summary": PRODUCTION.PEAK_SUMMARY, "peak_policy": PRODUCTION.PEAK_POLICY,
        "processing_plan": PRODUCTION.PROCESSING_PLAN, "samples": PRODUCTION.SAMPLES,
        "filtering_policy": PRODUCTION.FILTER_CONFIG}.items()}
    value["scope"] = "portable_reproduction_not_historical_execution"
    sources = [ROOT / "workflow/scripts" / (name + ".py") for name in
               ("chipseq_portable_runtime", "chipseq_production", "chipseq_alignment_support",
                "chipseq_filtering_support", "chipseq_peak_calling_dynamic", "chipseq_benchmark",
                "chipseq_download_fastq")]
    sources += [ROOT / "workflow/envs/chipseq_production.yaml"]
    value["implementation_sha256"] = {str(p.relative_to(ROOT)): sha(p) for p in sources}
    return value


def scratch_context(settings, raw):
    parent = settings.get("scratch_root") or os.environ.get("TMPDIR") or tempfile.gettempdir()
    parent = location(parent) if settings.get("scratch_root") else Path(parent).resolve()
    if not parent.is_dir():
        raise ValueError("Scratch directory does not exist")
    needed = max(8 * 1024**3, 4 * Path(raw).stat().st_size + 2 * 1024**3)
    if shutil.disk_usage(parent).free < needed:
        raise ValueError("Insufficient scratch space")
    return tempfile.TemporaryDirectory(prefix="chipseq-reproduction-", dir=parent)


def process_run(plan, role, raw, work):
    run_id = plan["analysis"]["ip_run_accession" if role == "ip" else "control_run_accession"]
    raw, work = Path(raw).resolve(), Path(work).resolve()
    download = module("chipseq_download_fastq")
    sample = plan["samples"][role]
    if not download.validate_existing_file(raw, sample["fastq_bytes"], sample["fastq_md5"]):
        raise ValueError("FASTQ no longer matches ENA evidence")
    for name in ("config", "inputs/reference", "filtered", "alignment", "fastqc", "picard_tmp"):
        (work / name).mkdir(parents=True, exist_ok=True)
    for name in ("chipseq_alignment.json", "chipseq_filtering.json"):
        shutil.copyfile(ROOT / "config" / name, work / "config" / name)
    ref = Path(plan["reference"]["root"])
    shutil.copyfile(ref / "genome_plus_mt.fa.fai", work / "inputs/reference/genome_plus_mt.fa.fai")
    processed = work / "processed.fastq.gz"; filtered = work / "filtered"; align = work / "alignment"
    run(["fastp", "--in1", raw, "--out1", processed, "--json", work / "fastp.json", "--html", work / "fastp.html", "--thread", "4", "-Q", "-L", "-G"], work / "fastp.log")
    with gzip.open(processed, "rb") as handle:
        while handle.read(1024 * 1024):
            pass
    run(["fastqc", "--threads", "1", "--outdir", work / "fastqc", processed], work / "fastqc.log")
    PRODUCTION.alignment_plan(plan, run_id, role, processed, work / "fastp.json", work / "config/chipseq_alignment_inputs.json")
    bam = align / "raw.sorted.bam"
    with (align / "bowtie2.log").open("w") as err, (align / "sort.log").open("w") as sorterr:
        mapper = subprocess.Popen(bowtie2_command(run_id, processed, ref / "bowtie2_index/genome_plus_mt"), stdout=subprocess.PIPE, stderr=err)
        try:
            sorter = subprocess.run(["samtools", "sort", "-@", "1", "-m", "768M", "-T", str(work / "sort_tmp"), "-o", str(bam), "-"], stdin=mapper.stdout, stderr=sorterr)
            mapper.stdout.close()
            status = mapper.wait()
            if sorter.returncode or status:
                raise RuntimeError("Bowtie2/samtools pipeline failed")
        finally:
            if mapper.poll() is None:
                mapper.kill(); mapper.wait()
    run(["samtools", "quickcheck", "-v", bam]); run(["samtools", "index", "-c", bam])
    for command, name in (("flagstat", "flagstat.txt"), ("idxstats", "idxstats.tsv"), ("stats", "samtools_stats.txt")):
        run(["samtools", command, bam], output=align / name)
    with cwd(work):
        helper("chipseq_alignment_support", "qc", "--run", run_id, "--bam", bam, "--output", align / "alignment_qc.json")
    processed.unlink()
    PRODUCTION.filtering_plan(plan, run_id, role, bam, align / "alignment_qc.json", work / "config/chipseq_filtering_inputs.json")
    marked = work / "dupmarked.sorted.bam"
    run(picard_command(bam, marked, filtered / "picard_metrics.txt", work / "picard_tmp"), filtered / "picard.log")
    run(["samtools", "quickcheck", "-v", marked]); bam.unlink(); Path(str(bam) + ".csi").unlink()
    with cwd(work):
        helper("chipseq_filtering_support", "filter", "--run", run_id, "--bam", marked, "--metrics", filtered / "picard_metrics.txt",
               "--output", filtered / "filtered.bam", "--report", filtered / "filtering_qc.json", "--flow", filtered / "filtering_flow.tsv")
        run(["samtools", "quickcheck", "-v", filtered / "filtered.bam"]); run(["samtools", "index", "-c", filtered / "filtered.bam"])
        helper("chipseq_filtering_support", "verify", "--run", run_id, "--bam", filtered / "filtered.bam", "--report", filtered / "filtering_qc.json", "--output", filtered / "filtered_validation.json")
    marked.unlink()
    for f in (align / "alignment_qc.json", work / "fastp.json", work / "fastp.html"):
        shutil.copyfile(f, filtered / f.name)
    shutil.copytree(work / "fastqc", filtered / "fastqc")
    return filtered


def process_control(settings, control, raw):
    validate_contract(); versions = assert_versions()
    peaks, _ = PRODUCTION.validated_plans()
    choices = [r for r in peaks if r["control_run_accession"] == control]
    if not choices:
        raise ValueError("Unknown authoritative Input")
    plan = execution_plan(choices[0]["analysis_id"], location(settings["reference_root"]))
    root = location(settings["output_root"]); dest = root / "shared_controls" / control
    with production_slot(settings), PRODUCTION.control_lock(root, control, exclusive=True):
        if dest.exists():
            if (dest / "control_manifest.json").exists():
                PRODUCTION.load_control(plan, dest); return
            remove_output_scaffold(dest)
        with scratch_context(settings, raw) as work:
            source = process_run(plan, "input", raw, work)
            write_json(source / "software_versions.json", versions)
            PRODUCTION.persist_control(plan, source, dest)
            shutil.copyfile(source / "software_versions.json", dest / "software_versions.json")


def overlap_count(bam, peaks):
    if not Path(peaks).stat().st_size:
        return 0
    proc = subprocess.Popen(["bedtools", "intersect", "-abam", str(bam), "-b", str(peaks), "-u", "-bed"], stdout=subprocess.PIPE, text=True)
    try:
        count = sum(1 for _ in proc.stdout)
    finally:
        proc.stdout.close()
    if proc.wait():
        raise RuntimeError("FRiP overlap counting failed")
    return count


def process_analysis(settings, analysis, raw):
    validate_contract(); versions = assert_versions()
    plan = execution_plan(analysis, location(settings["reference_root"]))
    root = location(settings["output_root"]); output = root / "analyses" / analysis
    remove_output_scaffold(output, ("", "artifacts", "artifacts/peaks"))
    control = plan["analysis"]["control_run_accession"]
    started = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with production_slot(settings), PRODUCTION.control_lock(root, control, exclusive=False), scratch_context(settings, raw) as work:
        work = Path(work); ip = process_run(plan, "ip", raw, work / "ip")
        control_dir = root / "shared_controls" / control
        runtime = work / "runtime.json"; PRODUCTION.make_runtime(plan, ip, control_dir, runtime)
        art = work / "artifacts"; peaks = art / "peaks"; fragment = art / "fragment"; frip = art / "frip"
        for d in (peaks, fragment, frip):
            d.mkdir(parents=True)
        fragment_path = None
        if plan["analysis"]["fragment_size_policy"] == "phantompeakqualtools":
            run(["run_spp.R", "-c=" + str(ip / "filtered.bam"), "-p=4", "-savp=" + str(fragment / "cross_correlation.pdf"), "-out=" + str(fragment / "phantompeakqualtools.tsv")], fragment / "phantompeakqualtools.log")
            if not (fragment / "cross_correlation.pdf").stat().st_size:
                raise ValueError("Missing cross-correlation plot")
            fragment_path = fragment / "fragment_estimate.json"
            helper("chipseq_peak_calling_dynamic", "parse-phantom", "--table", fragment / "phantompeakqualtools.tsv", "--runtime", runtime, "--analysis", analysis, "--output", fragment_path)
        parameters = art / "parameters.json"
        write_json(parameters, DYNAMIC.analysis_parameters(runtime, analysis, fragment_path))
        PRODUCTION.fragment_record(plan, parameters, fragment / "fragment_size_provenance.json")
        DYNAMIC.run_macs3(runtime, parameters, analysis, peaks, art / "macs3.log", art / "macs3_version.txt", fragment_path)
        narrow = plan["analysis"]["peak_mode"] == "narrow"
        primary = peaks / (analysis + ("_peaks.narrowPeak" if narrow else "_peaks.broadPeak"))
        secondary = peaks / (analysis + ("_summits.bed" if narrow else "_peaks.gappedPeak"))
        for role, bam in (("ip", ip / "filtered.bam"), ("input", control_dir / "filtered.bam")):
            count = subprocess.check_output(["samtools", "view", "-c", str(bam)], text=True)
            (frip / (role + "_total_reads.txt")).write_text(count)
            name = "ip_reads_in_peaks.txt" if role == "ip" else "input_reads_overlapping_ip_peaks.txt"
            (frip / name).write_text(str(overlap_count(bam, primary)) + "\n")
        treat = peaks / (analysis + "_treat_pileup.bdg"); lam = peaks / (analysis + "_control_lambda.bdg")
        DYNAMIC.summarize_analysis(runtime, parameters, analysis, primary, secondary, peaks / (analysis + "_peaks.xls"), treat, lam,
                                   art / "macs3.log", art / "macs3_version.txt", frip / "ip_total_reads.txt", frip / "input_total_reads.txt",
                                   frip / "ip_reads_in_peaks.txt", frip / "input_reads_overlapping_ip_peaks.txt", art / "peak_qc.json", art / "provenance.json", fragment_path)
        compressed = []
        for f in (treat, lam):
            z = Path(str(f) + ".gz")
            run(["gzip", "-n", "-c", f], output_binary=z)
            compressed.append(z)
        PRODUCTION.compression_record([treat, lam], compressed, art / "signal_compression.json")
        treat.unlink(); lam.unlink()
        shutil.copytree(ip, art / "ip_processing_qc", ignore=shutil.ignore_patterns("filtered.bam", "filtered.bam.csi"))
        shutil.copyfile(runtime, art / "runtime_manifest.json")
        write_json(art / "execution_plan.json", plan); write_json(art / "software_versions.json", versions)
        output.parent.mkdir(parents=True, exist_ok=True)
        stage = output.with_name("." + analysis + ".part")
        if stage.exists():
            raise ValueError("Analysis staging collision")
        stage.mkdir()
        try:
            shutil.copytree(art, stage / "artifacts")
            PRODUCTION.finalize(plan, stage / "artifacts", stage, started, "local_reproduction")
            stage.replace(output)
        except BaseException:
            shutil.rmtree(stage, ignore_errors=True); raise


def cohort_summary(settings, output):
    peaks, _ = PRODUCTION.validated_plans(); root = location(settings["output_root"])
    records = []
    for row in peaks:
        PRODUCTION.verify_analysis(root, row["analysis_id"])
        records.append(read_json(root / "analyses" / row["analysis_id"] / "artifacts/peak_qc.json"))
    value = DYNAMIC.cohort_summary(records)
    value["scope"] = "portable_reproduction_not_historical_production"
    write_json(output, value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/chipseq.yaml")
    parser.add_argument("--output-root")
    parser.add_argument("--reference-root")
    parser.add_argument("--scratch-root")
    parser.add_argument("--max-parallel-analyses", type=int)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("validate-plan", "fastq-manifest", "reference", "summary"):
        s = sub.add_parser(name); s.add_argument("--output", required=True)
    s = sub.add_parser("compare-plans"); s.add_argument("--analysis", required=True); s.add_argument("--processing", required=True); s.add_argument("--output", required=True)
    s = sub.add_parser("metadata-record")
    for name in ("plan", "outdir", "accession"):
        s.add_argument("--" + name, required=True)
    for name, identity in (("control", "control"), ("analysis", "analysis")):
        s = sub.add_parser(name); s.add_argument("--" + identity, required=True); s.add_argument("--raw", required=True)
    args = parser.parse_args()
    overrides = {key: getattr(args, key) for key in
                 ("output_root", "reference_root", "scratch_root", "max_parallel_analyses")
                 if getattr(args, key) is not None}
    settings = configuration(args.config, overrides)
    if args.command == "validate-plan": write_json(args.output, validate_contract())
    elif args.command == "fastq-manifest": fastq_manifest(args.output)
    elif args.command == "compare-plans": write_json(args.output, compare_plans(args.analysis, args.processing))
    elif args.command == "metadata-record": metadata_record(args.plan, args.outdir, args.accession)
    elif args.command == "reference": reference_build(args.output)
    elif args.command == "control": process_control(settings, args.control, args.raw)
    elif args.command == "analysis": process_analysis(settings, args.analysis, args.raw)
    else: cohort_summary(settings, args.output)


if __name__ == "__main__":
    main()
