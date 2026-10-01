#!/usr/bin/env python3
"""Fail-closed contracts for dynamic ChIP-seq peak calling and descriptive QC."""

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

MTDNA = "NC_007936.1"
ASSEMBLY_ACCESSION = "GCF_003668045.3"
ASSEMBLY_NAME = "CriGri-PICRH-1.0"
GENOME_SIZE_POLICY = (
    "derive_nuclear_reference_span_from_fai_excluding_mitochondrial_contig"
)
SAFE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")
DIGEST = re.compile(r"[0-9a-f]{64}")
PHANTOM_FIELDS = (
    "filename", "num_reads", "estimated_fragment_length", "fragment_correlation",
    "phantom_peak", "phantom_correlation", "minimum_cross_correlation_shift",
    "minimum_cross_correlation", "nsc", "rsc", "quality_tag",
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".part", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(json.dumps(value, indent=2, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def canonical_sha256(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def safe_identifier(value, label):
    if not isinstance(value, str) or not SAFE_ID.fullmatch(value) or value in {".", ".."}:
        raise ValueError(f"Unsafe {label}: {value!r}")
    return value


def valid_digest(value, label):
    if not isinstance(value, str) or not DIGEST.fullmatch(value):
        raise ValueError(f"Invalid SHA-256 for {label}")
    return value


def verify_hash(path, expected, label=None):
    valid_digest(expected, label or str(path))
    observed = sha256(path)
    if observed != expected:
        raise ValueError(
            f"SHA-256 mismatch for {label or path}: expected {expected}, observed {observed}"
        )
    return observed


def parse_fai(path, mitochondrial_accession=MTDNA):
    """Validate an FAI and return contig lengths and nuclear-span provenance."""
    path = Path(path)
    contigs = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                raise ValueError(f"Blank FAI line {line_number}")
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 5 or not fields[0]:
                raise ValueError(
                    f"Malformed FAI line {line_number}: expected five tab-separated fields"
                )
            name = fields[0]
            if name in contigs:
                raise ValueError(f"Duplicate FAI contig: {name}")
            try:
                length, offset, line_bases, line_width = map(int, fields[1:])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid FAI numeric field at line {line_number}") from exc
            if length <= 0 or str(length) != fields[1]:
                raise ValueError(f"Invalid FAI length at line {line_number}")
            if offset < 0 or str(offset) != fields[2]:
                raise ValueError(f"Invalid FAI offset at line {line_number}")
            if line_bases <= 0 or str(line_bases) != fields[3]:
                raise ValueError(f"Invalid FAI line bases at line {line_number}")
            if line_width < line_bases or str(line_width) != fields[4]:
                raise ValueError(f"Invalid FAI line geometry at line {line_number}")
            contigs[name] = length

    if mitochondrial_accession not in contigs:
        raise ValueError(f"Mitochondrial accession {mitochondrial_accession} absent from FAI")
    nuclear_span = sum(
        length for name, length in contigs.items() if name != mitochondrial_accession
    )
    if nuclear_span <= 0:
        raise ValueError("FAI contains no positive nuclear reference span")
    return {
        "fai_sha256": sha256(path),
        "contigs": contigs,
        "sequence_count": len(contigs),
        "total_span": sum(contigs.values()),
        "nuclear_span": nuclear_span,
        "mitochondrial_accession": mitochondrial_accession,
        "genome_size_policy": GENOME_SIZE_POLICY,
        "genome_size_description": (
            "nuclear-reference-span approximation for MACS3 genome size; "
            "not a mappability-derived effective genome size"
        ),
    }


def validate_reference(fai_path, metadata, expected_fai_sha256, fasta_path=None):
    result = parse_fai(fai_path)
    if not isinstance(metadata, dict):
        raise ValueError("Reference metadata must be an object")
    required = {
        "nuclear_accession", "nuclear_sha256", "mitochondrial_accession",
        "mitochondrial_sha256", "mapping_sha256", "nuclear_sequences", "purpose",
    }
    if not required.issubset(metadata):
        raise ValueError("Reference provenance is missing mandatory upstream identity fields")
    if metadata["nuclear_accession"] != ASSEMBLY_ACCESSION:
        raise ValueError("Unexpected reference nuclear accession")
    if metadata["mitochondrial_accession"] != MTDNA:
        raise ValueError("Unexpected reference mitochondrial accession")
    for key in ("nuclear_sha256", "mitochondrial_sha256", "mapping_sha256"):
        valid_digest(metadata[key], f"reference provenance {key}")
    if type(metadata["nuclear_sequences"]) is not int or metadata["nuclear_sequences"] <= 0:
        raise ValueError("Invalid nuclear sequence count in reference provenance")
    if metadata["purpose"] != "ChIP-seq nuclear and mitochondrial mapping":
        raise ValueError("Unexpected reference provenance purpose")
    verify_hash(fai_path, expected_fai_sha256, "reference FAI")
    if fasta_path is not None:
        fasta_path = Path(fasta_path)
        if not fasta_path.is_file():
            raise ValueError(f"Missing reference FASTA: {fasta_path}")
        result["fasta_sha256"] = sha256(fasta_path)
        if "fasta_sha256" in metadata:
            verify_hash(fasta_path, metadata["fasta_sha256"], "reference FASTA")
    result["assembly_accession"] = ASSEMBLY_ACCESSION
    result["assembly_name"] = ASSEMBLY_NAME
    result["source_reference_metadata"] = metadata
    return result


def _finite_float(value, label, allow_na=False):
    if allow_na and value in {"NA", "NaN", "nan", ""}:
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Malformed PhantomPeakQualTools {label}: {value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"Non-finite PhantomPeakQualTools {label}")
    return parsed


def _comma_list(value, label):
    if not isinstance(value, str) or not value or any(not item for item in value.split(",")):
        raise ValueError(f"Empty or malformed PhantomPeakQualTools {label} list")
    return value.split(",")


def _quality_tag(value):
    if not re.fullmatch(r"-2|-1|0|1|2", value):
        raise ValueError("Malformed PhantomPeakQualTools quality tag")
    return int(value)


def parse_phantompeak_table(
    path, expected_run, expected_bam_name=None, analysis_id=None,
    expected_bam_sha256=None,
):
    """Parse one formal 11-column PhantomPeakQualTools result row."""
    safe_identifier(expected_run, "PhantomPeakQualTools run")
    rows = []
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for raw in csv.reader(handle, delimiter="\t"):
            if not raw or not any(cell.strip() for cell in raw) or raw[0].startswith("#"):
                continue
            if tuple(cell.strip().lower() for cell in raw) == PHANTOM_FIELDS:
                continue
            if len(raw) != len(PHANTOM_FIELDS):
                raise ValueError("Malformed PhantomPeakQualTools row: expected 11 columns")
            rows.append([cell.strip() for cell in raw])
    if len(rows) != 1:
        raise ValueError(
            "PhantomPeakQualTools result must contain exactly one unambiguous data row"
        )
    values = dict(zip(PHANTOM_FIELDS, rows[0]))
    observed_name = Path(values["filename"]).name
    expected_name = Path(expected_bam_name).name if expected_bam_name else None
    if expected_name:
        matches = observed_name == expected_name
    else:
        matches = re.search(
            rf"(^|[^A-Za-z0-9]){re.escape(expected_run)}([^A-Za-z0-9]|$)", observed_name
        ) is not None
    if not matches:
        raise ValueError("PhantomPeakQualTools result does not match the expected IP BAM")
    fragment_texts = _comma_list(
        values["estimated_fragment_length"], "fragment candidate"
    )
    correlation_texts = _comma_list(
        values["fragment_correlation"], "candidate correlation"
    )
    if len(fragment_texts) != len(correlation_texts):
        raise ValueError(
            "PhantomPeakQualTools fragment candidates and correlations have unequal cardinality"
        )
    if any(not re.fullmatch(r"[1-9][0-9]*", item) for item in fragment_texts):
        raise ValueError("PhantomPeakQualTools fragment candidates must be positive integers")
    fragments = [int(item) for item in fragment_texts]
    correlations = [
        _finite_float(item, "candidate correlation") for item in correlation_texts
    ]
    num_reads = values["num_reads"]
    if not re.fullmatch(r"[1-9][0-9]*", num_reads):
        raise ValueError("PhantomPeakQualTools numReads must be positive")
    if analysis_id is None:
        analysis_id = expected_run
    safe_identifier(analysis_id, "PhantomPeakQualTools analysis_id")
    bam_identity = str(Path(expected_bam_name).resolve()) if expected_bam_name else None
    if expected_bam_sha256 is not None:
        if expected_bam_name is None:
            raise ValueError("BAM path is required when binding a Phantom result hash")
        verify_hash(expected_bam_name, expected_bam_sha256, "PhantomPeakQualTools IP BAM")
    result = {
        "schema_version": 2,
        "scope": "phantompeakqualtools_fragment_estimate",
        "analysis_id": analysis_id,
        "run_accession": expected_run,
        "ip_bam_path": bam_identity,
        "ip_bam_sha256": expected_bam_sha256,
        "source_filename": values["filename"],
        "num_reads": int(num_reads),
        "num_reads_semantics": (
            "PhantomPeakQualTools-reported eligible reads; preserved separately from "
            "the filtering retained-read count because formal semantics do not justify equality"
        ),
        "selected_fragment_candidate_bp": fragments[0],
        "fragment_size_bp": fragments[0],
        "fragment_candidates_bp": fragments,
        "fragment_candidate_correlations": correlations,
        "fragment_correlation": correlations[0],
        "phantom_peak": _finite_float(values["phantom_peak"], "phantom peak", True),
        "phantom_correlation": _finite_float(values["phantom_correlation"], "phantom correlation", True),
        "minimum_cross_correlation_shift": _finite_float(values["minimum_cross_correlation_shift"], "minimum shift", True),
        "minimum_cross_correlation": _finite_float(values["minimum_cross_correlation"], "minimum correlation", True),
        "nsc": _finite_float(values["nsc"], "NSC", True),
        "rsc": _finite_float(values["rsc"], "RSC", True),
        "quality_tag": _quality_tag(values["quality_tag"]),
        "quality_tag_interpretation": "descriptive_tool_output_not_a_universal_pass_fail_threshold",
        "source_result_path": str(Path(path).resolve()),
        "source_result_sha256": sha256(path),
    }
    result["record_sha256"] = canonical_sha256(result)
    return result


def validate_phantom_result(
    plan_row, phantom_result, runtime_run=None, fragment_path=None, raw_table_path=None,
):
    if not isinstance(phantom_result, dict):
        raise ValueError("PRJEB9291 requires its own valid fragment estimate; no fallback")
    if (phantom_result.get("schema_version") != 2
            or phantom_result.get("scope") != "phantompeakqualtools_fragment_estimate"):
        raise ValueError("Unsupported fragment-estimation result schema")
    required = {
        "analysis_id", "run_accession", "ip_bam_path", "ip_bam_sha256",
        "source_filename", "num_reads", "num_reads_semantics",
        "selected_fragment_candidate_bp", "fragment_size_bp",
        "fragment_candidates_bp", "fragment_candidate_correlations",
        "fragment_correlation", "phantom_peak", "phantom_correlation",
        "minimum_cross_correlation_shift", "minimum_cross_correlation", "nsc", "rsc",
        "quality_tag", "quality_tag_interpretation", "source_result_path",
        "source_result_sha256", "record_sha256",
    }
    if not required.issubset(phantom_result):
        raise ValueError("Fragment-estimation result was tampered or is incomplete")
    recorded = phantom_result.get("record_sha256")
    payload = {key: value for key, value in phantom_result.items() if key != "record_sha256"}
    if not isinstance(recorded, str) or recorded != canonical_sha256(payload):
        raise ValueError("Fragment-estimation result was tampered or is incomplete")
    if phantom_result["analysis_id"] != plan_row.get("analysis_id"):
        raise ValueError("Fragment estimate analysis_id does not match analysis")
    if phantom_result["run_accession"] != plan_row.get("ip_run_accession"):
        raise ValueError("Fragment estimate run does not match analysis IP")
    if type(phantom_result["num_reads"]) is not int or phantom_result["num_reads"] <= 0:
        raise ValueError("Invalid fragment-estimation read count")
    candidates = phantom_result["fragment_candidates_bp"]
    correlations = phantom_result["fragment_candidate_correlations"]
    if (not isinstance(candidates, list) or not candidates
            or any(type(value) is not int or value <= 0 for value in candidates)):
        raise ValueError("Invalid external fragment candidate list; no fallback")
    if (not isinstance(correlations, list) or len(correlations) != len(candidates)
            or any(type(value) not in {int, float} or not math.isfinite(value)
                   for value in correlations)):
        raise ValueError("Invalid external fragment-correlation list; no fallback")
    selected = phantom_result["selected_fragment_candidate_bp"]
    if (type(selected) is not int or selected != candidates[0]
            or phantom_result["fragment_size_bp"] != selected
            or phantom_result["fragment_correlation"] != correlations[0]):
        raise ValueError("Selected fragment must be the documented first candidate")
    _quality_tag(str(phantom_result["quality_tag"]))
    for key in (
        "phantom_peak", "phantom_correlation", "minimum_cross_correlation_shift",
        "minimum_cross_correlation", "nsc", "rsc",
    ):
        value = phantom_result[key]
        if value is not None and (type(value) not in {int, float} or not math.isfinite(value)):
            raise ValueError(f"Invalid PhantomPeakQualTools field: {key}")
    valid_digest(phantom_result["source_result_sha256"], "PhantomPeakQualTools table")
    valid_digest(phantom_result["ip_bam_sha256"], "PhantomPeakQualTools IP BAM")
    if runtime_run is not None:
        if runtime_run.get("run_accession") != plan_row.get("ip_run_accession"):
            raise ValueError("Runtime IP run identity mismatch")
        if phantom_result["ip_bam_sha256"] != runtime_run.get("bam_sha256"):
            raise ValueError("Fragment estimate IP BAM SHA-256 mismatch")
        if Path(phantom_result["ip_bam_path"]) != Path(runtime_run.get("bam_path", "")):
            raise ValueError("Fragment estimate IP BAM identity mismatch")
        verify_hash(runtime_run["bam_path"], runtime_run["bam_sha256"], "fragment IP BAM")
    if raw_table_path is not None:
        expected_path = Path(raw_table_path).resolve()
        if Path(phantom_result["source_result_path"]) != expected_path:
            raise ValueError("Fragment result points to the wrong raw Phantom table")
        verify_hash(expected_path, phantom_result["source_result_sha256"], "raw Phantom table")
        reparsed = parse_phantompeak_table(
            expected_path,
            plan_row["ip_run_accession"],
            runtime_run["bam_path"] if runtime_run else phantom_result["ip_bam_path"],
            plan_row["analysis_id"],
            runtime_run["bam_sha256"] if runtime_run else phantom_result["ip_bam_sha256"],
        )
        if reparsed != phantom_result:
            raise ValueError("Fragment record does not match the bound raw Phantom table")
    return selected


def resolve_fragment_size(
    plan_row, phantom_result=None, runtime_run=None, fragment_path=None,
    raw_table_path=None,
):
    policy = plan_row.get("fragment_size_policy")
    study = plan_row.get("study_accession")
    if policy == "fixed":
        if study != "PRJNA865478" or str(plan_row.get("fragment_size_bp")) != "147":
            raise ValueError("Fixed 147 bp policy is restricted to PRJNA865478")
        if phantom_result is not None:
            raise ValueError("PRJNA865478 must not consume PhantomPeakQualTools results")
        return 147, "study_fixed_publication_supported"
    if policy == "phantompeakqualtools":
        if study != "PRJEB9291":
            raise ValueError("PhantomPeakQualTools policy is restricted to PRJEB9291")
        value = validate_phantom_result(
            plan_row, phantom_result, runtime_run, fragment_path, raw_table_path
        )
        return value, "phantompeakqualtools_external_estimate"
    raise ValueError(f"Unsupported fragment-size policy: {policy!r}")


def macs3_contract(
    plan_row, nuclear_span, phantom_result=None, runtime_run=None,
    fragment_path=None, raw_table_path=None,
):
    if type(nuclear_span) is not int or nuclear_span <= 0:
        raise ValueError("MACS3 requires a positive FAI-derived nuclear span")
    fragment_size, fragment_source = resolve_fragment_size(
        plan_row, phantom_result, runtime_run, fragment_path, raw_table_path
    )
    mode = plan_row.get("peak_mode")
    if mode not in {"narrow", "broad"}:
        raise ValueError(f"Unsupported peak mode: {mode!r}")
    if plan_row.get("keep_dup") != "all":
        raise ValueError("Dynamic MACS3 requires --keep-dup all")
    if plan_row.get("format") != "BAM":
        raise ValueError("Dynamic MACS3 supports BAM input only")
    if plan_row.get("scale_to") not in {"small", "large"}:
        raise ValueError("Unsupported MACS3 scale policy")
    if plan_row.get("ip_run_accession") == plan_row.get("control_run_accession"):
        raise ValueError("MACS3 IP and Input must be distinct")
    if plan_row.get("effective_genome_size_policy") != GENOME_SIZE_POLICY:
        raise ValueError("Unsupported MACS3 genome-size policy")
    booleans = {}
    for key in ("spmr", "store_bdg", "call_summits", "cutoff_analysis"):
        value = plan_row.get(key)
        if value not in {"true", "false"}:
            raise ValueError(f"Plan field {key} must be normalized true/false")
        booleans[key] = value == "true"
    if booleans["spmr"] and not booleans["store_bdg"]:
        raise ValueError("SPMR requires bedGraph output")
    try:
        qvalue = float(plan_row["qvalue"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Invalid MACS3 q-value") from exc
    if not 0 < qvalue < 1:
        raise ValueError("MACS3 q-value must be between zero and one")
    flags = ["--nomodel", "--extsize", str(fragment_size), "--keep-dup", "all"]
    if booleans["store_bdg"]:
        flags.append("-B")
    if booleans["spmr"]:
        flags.append("--SPMR")
    if booleans["call_summits"] and mode == "narrow":
        flags.append("--call-summits")
    if booleans["cutoff_analysis"]:
        flags.append("--cutoff-analysis")
    outputs = ["narrowPeak", "summits"] if mode == "narrow" else ["broadPeak", "gappedPeak"]
    if mode == "broad":
        cutoff = plan_row.get("broad_cutoff")
        try:
            cutoff_number = float(cutoff)
        except (TypeError, ValueError) as exc:
            raise ValueError("Broad mode requires a valid broad cutoff") from exc
        if not 0 < cutoff_number < 1:
            raise ValueError("Broad cutoff must be between zero and one")
        flags.extend(["--broad", "--broad-cutoff", str(cutoff)])
    return {
        "analysis_id": safe_identifier(plan_row.get("analysis_id"), "analysis_id"),
        "ip_run_accession": plan_row.get("ip_run_accession"),
        "control_run_accession": plan_row.get("control_run_accession"),
        "peak_mode": mode,
        "format": "BAM",
        "qvalue": qvalue,
        "scale_to": plan_row["scale_to"],
        "keep_dup": "all",
        "macs3_genome_size": nuclear_span,
        "macs3_genome_size_basis": GENOME_SIZE_POLICY,
        "fragment_size_bp": fragment_size,
        "fragment_size_source": fragment_source,
        "flags": flags,
        "required_peak_outputs": outputs,
    }


def _number(value, label, integer=False):
    try:
        result = int(value) if integer else float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Malformed {label}: {value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"Non-finite {label}")
    return result


def parse_peak_file(path, mode, contigs, mitochondrial_accession=MTDNA):
    """Validate narrowPeak/broadPeak/gappedPeak coordinates and summarize peaks."""
    minimum_columns = {"narrow": 10, "broad": 9, "gapped": 15}
    if mode not in minimum_columns:
        raise ValueError(f"Unsupported peak-file mode: {mode}")
    if mitochondrial_accession not in contigs:
        raise ValueError("Peak validation contigs omit the mitochondrial accession")
    intervals = defaultdict(list)
    widths, signals, q_scores = [], [], []
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#") or raw.startswith("track") or raw.startswith("browser"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < minimum_columns[mode]:
                raise ValueError(f"Malformed {mode} peak line {line_number}")
            chrom = fields[0]
            if chrom not in contigs:
                raise ValueError(f"Unknown peak contig at line {line_number}: {chrom}")
            if chrom == mitochondrial_accession:
                raise ValueError("Mitochondrial peak detected")
            start = _number(fields[1], "peak start", True)
            end = _number(fields[2], "peak end", True)
            if start < 0 or end <= start or end > contigs[chrom]:
                raise ValueError(f"Peak coordinates out of range at line {line_number}")
            signal_column, q_column = (12, 14) if mode == "gapped" else (6, 8)
            signal = _number(fields[signal_column], "peak signal")
            q_score = _number(fields[q_column], "peak q-value score")
            if signal < 0 or q_score < 0:
                raise ValueError(f"Negative peak metric at line {line_number}")
            if mode == "narrow":
                summit = _number(fields[9], "summit offset", True)
                if summit < 0 or summit >= end - start:
                    raise ValueError(f"Summit outside peak at line {line_number}")
            if mode == "gapped":
                block_count = _number(fields[9], "block count", True)
                sizes = [x for x in fields[10].rstrip(",").split(",") if x]
                starts = [x for x in fields[11].rstrip(",").split(",") if x]
                if block_count <= 0 or len(sizes) != block_count or len(starts) != block_count:
                    raise ValueError(f"Malformed gappedPeak blocks at line {line_number}")
                for size, offset in zip(sizes, starts):
                    size, offset = _number(size, "block size", True), _number(offset, "block offset", True)
                    if size <= 0 or offset < 0 or offset + size > end - start:
                        raise ValueError(f"GappedPeak block outside peak at line {line_number}")
            intervals[chrom].append((start, end))
            widths.append(end - start)
            signals.append(signal)
            q_scores.append(q_score)
    union_bp = 0
    for values in intervals.values():
        values.sort()
        if not values:
            continue
        left, right = values[0]
        for start, end in values[1:]:
            if start <= right:
                right = max(right, end)
            else:
                union_bp += right - left
                left, right = start, end
        union_bp += right - left
    def median(values):
        ordered = sorted(values)
        n = len(ordered)
        return None if not n else (ordered[n // 2] if n % 2 else (ordered[n // 2 - 1] + ordered[n // 2]) / 2)
    return {
        "peak_count": len(widths),
        "chromosomes_with_peaks": len(intervals),
        "sum_peak_bp": sum(widths),
        "union_peak_bp": union_bp,
        "mean_peak_width_bp": sum(widths) / len(widths) if widths else None,
        "median_peak_width_bp": median(widths),
        "median_fold_enrichment": median(signals),
        "median_minus_log10_qvalue": median(q_scores),
    }


def descriptive_frip(ip_total, input_total, ip_overlap, input_overlap):
    values = (ip_total, input_total, ip_overlap, input_overlap)
    if any(type(value) is not int or value < 0 for value in values):
        raise ValueError("FRiP counts must be non-negative integers")
    if ip_overlap > ip_total or input_overlap > input_total:
        raise ValueError("Peak-overlap count exceeds its read denominator")
    ip_fraction = ip_overlap / ip_total if ip_total else None
    input_fraction = input_overlap / input_total if input_total else None
    ratio = (
        ip_fraction / input_fraction
        if ip_fraction is not None and input_fraction not in {None, 0}
        else None
    )
    return {
        "ip_total_reads": ip_total,
        "ip_reads_in_peaks": ip_overlap,
        "ip_frip": ip_fraction,
        "input_total_reads": input_total,
        "input_reads_overlapping_ip_peaks": input_overlap,
        "input_peak_overlap_fraction": input_fraction,
        "ip_to_input_peak_overlap_ratio": ratio,
        "interpretation": "descriptive_only_no_universal_pass_fail_threshold",
    }


def validate_output_destination(path):
    path = Path(path)
    if path.exists():
        raise ValueError(f"Refusing unsafe overwrite of existing output: {path}")
    if path.name in {"", ".", ".."}:
        raise ValueError(f"Unsafe output destination: {path}")
    return path


def read_peak_plan(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not reader.fieldnames or not rows:
        raise ValueError("Peak plan is empty or has no header")
    required = {
        "analysis_id", "study_accession", "ip_run_accession",
        "control_run_accession", "declared_target", "peak_mode",
        "fragment_size_policy", "fragment_size_bp", "qvalue", "broad_cutoff",
        "keep_dup", "scale_to", "spmr", "store_bdg", "call_summits",
        "cutoff_analysis", "effective_genome_size_policy", "execution_status",
    }
    missing = required - set(reader.fieldnames)
    if missing:
        raise ValueError(f"Peak plan missing required fields: {sorted(missing)}")
    seen_ids, seen_ips = set(), set()
    for row in rows:
        analysis_id = safe_identifier(row["analysis_id"], "analysis_id")
        ip = safe_identifier(row["ip_run_accession"], "IP run")
        safe_identifier(row["control_run_accession"], "Input run")
        safe_identifier(row["study_accession"], "study")
        if analysis_id in seen_ids:
            raise ValueError(f"Duplicate peak-plan analysis_id: {analysis_id}")
        if ip in seen_ips:
            raise ValueError(f"Duplicate peak-plan IP analysis: {ip}")
        seen_ids.add(analysis_id)
        seen_ips.add(ip)
    return rows


def checksum_manifest(path):
    result = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split(None, 1)
            if len(fields) != 2:
                raise ValueError(f"Malformed checksum manifest line {line_number}")
            digest, name = fields[0], fields[1].strip()
            valid_digest(digest, f"manifest line {line_number}")
            path_value = Path(name)
            if path_value.is_absolute() or ".." in path_value.parts or name in result:
                raise ValueError(f"Unsafe or duplicate checksum path: {name!r}")
            result[name] = digest
    if not result:
        raise ValueError("Empty checksum manifest")
    return result


def build_runtime_manifest(
    peak_plan_path,
    peak_summary_path,
    policy_path,
    analysis_plan_path,
    filtering_job_dir,
    verify_large_files=True,
):
    """Join execution artifacts to exact upstream pairs without changing pairing."""
    summary = json.loads(Path(peak_summary_path).read_text(encoding="utf-8"))
    expected_hashes = summary.get("input_sha256", {})
    verify_hash(analysis_plan_path, expected_hashes.get("analysis_plan", ""), "analysis plan")
    verify_hash(policy_path, expected_hashes.get("policy", ""), "peak policy")
    verify_hash(peak_plan_path, summary.get("peak_calling_plan_sha256", ""), "peak plan")
    rows = read_peak_plan(peak_plan_path)
    with Path(analysis_plan_path).open(encoding="utf-8", newline="") as handle:
        upstream_ready = [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row.get("analysis_status") == "ready"
        ]
    upstream_pairs = {
        row["analysis_id"]: (row["ip_run_accession"], row["control_run_accession"])
        for row in upstream_ready
    }
    planned_pairs = {
        row["analysis_id"]: (row["ip_run_accession"], row["control_run_accession"])
        for row in rows
    }
    if len(upstream_pairs) != len(upstream_ready) or planned_pairs != upstream_pairs:
        raise ValueError("Peak plan does not exactly inherit authoritative ready IP/Input pairs")

    filtering_job_dir = Path(filtering_job_dir).resolve()
    status_path = filtering_job_dir / "job_status.tsv"
    status = {}
    if not status_path.is_file():
        raise ValueError("Filtering job has no job_status.tsv")
    for raw in status_path.read_text(encoding="utf-8").splitlines():
        fields = raw.split("\t")
        if len(fields) != 2 or fields[0] in status:
            raise ValueError("Malformed filtering job status")
        status[fields[0]] = fields[1]
    if status.get("stage") != "completed" or status.get("exit_status") != "0":
        raise ValueError("Dynamic filtering job is not successfully completed")
    manifest_path = filtering_job_dir / "output.sha256"
    if not manifest_path.is_file():
        raise ValueError("Filtering job has no output checksum manifest")
    hashes = checksum_manifest(manifest_path)

    reference_name = "outputs/reference/genome_plus_mt.fa.fai"
    provenance_name = "outputs/reference/reference_provenance.json"
    for name in (
        reference_name,
        provenance_name,
        "outputs/input_provenance.json",
        "outputs/filtering_parameters.json",
    ):
        if name not in hashes or not (filtering_job_dir / name).is_file():
            raise ValueError(f"Filtering output is missing required provenance: {name}")
        verify_hash(filtering_job_dir / name, hashes[name], name)
    reference_metadata = json.loads(
        (filtering_job_dir / provenance_name).read_text(encoding="utf-8")
    )
    filtering_input_provenance = json.loads(
        (filtering_job_dir / "outputs/input_provenance.json").read_text(encoding="utf-8")
    )
    if filtering_input_provenance.get("schema_version") != 1:
        raise ValueError("Unsupported filtering input-provenance schema")
    if filtering_input_provenance.get("reference_provenance") != reference_metadata:
        raise ValueError("Filtering and reference provenance identities disagree")
    upstream_small = filtering_input_provenance.get("small_files")
    if not isinstance(upstream_small, dict):
        raise ValueError("Filtering input provenance lacks upstream reference hashes")
    for upstream_name in (reference_name, provenance_name):
        valid_digest(upstream_small.get(upstream_name), upstream_name)
        if upstream_small[upstream_name] != hashes[upstream_name]:
            raise ValueError(f"Filtering reference hash disagrees with upstream: {upstream_name}")
    valid_digest(
        filtering_input_provenance.get("source_manifest_sha256"),
        "upstream alignment output manifest",
    )
    reference = validate_reference(
        filtering_job_dir / reference_name, reference_metadata, hashes[reference_name]
    )
    reference["source_path"] = str((filtering_job_dir / reference_name).resolve())
    reference["provenance_sha256"] = hashes[provenance_name]

    required_runs = {
        row[key]
        for row in rows
        for key in ("ip_run_accession", "control_run_accession")
    }
    runs = {}
    for run in sorted(required_runs):
        names = {
            "bam": f"outputs/{run}/filtered.bam",
            "csi": f"outputs/{run}/filtered.bam.csi",
            "validation": f"outputs/{run}/filtered_validation.json",
            "qc": f"outputs/{run}/filtering_qc.json",
        }
        for label, name in names.items():
            path = filtering_job_dir / name
            if name not in hashes or not path.is_file():
                raise ValueError(f"Missing filtered {label} for {run}")
            if label not in {"bam", "csi"} or verify_large_files:
                verify_hash(path, hashes[name], name)
        validation = json.loads((filtering_job_dir / names["validation"]).read_text())
        qc = json.loads((filtering_job_dir / names["qc"]).read_text())
        if validation.get("run_accession") != run or qc.get("run_accession") != run:
            raise ValueError(f"Filtering run identity mismatch: {run}")
        if validation.get("full_read_verified") is not True:
            raise ValueError(f"Filtering full-read validation is absent: {run}")
        if hashes[names["bam"]] != validation.get("bam_sha256"):
            raise ValueError(f"Filtered BAM hash provenance mismatch: {run}")
        if hashes[names["csi"]] != validation.get("index_sha256"):
            raise ValueError(f"Filtered CSI hash provenance mismatch: {run}")
        retained = validation.get("retained_reads")
        if type(retained) is not int or retained <= 0 or retained != qc.get("retained_reads"):
            raise ValueError(f"Filtered retained-read provenance mismatch: {run}")
        runs[run] = {
            "run_accession": run,
            "role": qc.get("role"),
            "bam_path": str((filtering_job_dir / names["bam"]).resolve()),
            "csi_path": str((filtering_job_dir / names["csi"]).resolve()),
            "bam_sha256": hashes[names["bam"]],
            "csi_sha256": hashes[names["csi"]],
            "validation_sha256": hashes[names["validation"]],
            "qc_sha256": hashes[names["qc"]],
            "retained_reads": retained,
        }
    for row in rows:
        if runs[row["ip_run_accession"]]["role"] != "ip":
            raise ValueError(f"Authoritative IP has non-IP filtering role: {row['analysis_id']}")
        if runs[row["control_run_accession"]]["role"] != "input":
            raise ValueError(f"Authoritative Input has non-Input filtering role: {row['analysis_id']}")

    return {
        "schema_version": 1,
        "scope": "dynamic_chipseq_peak_calling_runtime",
        "analysis_count": len(rows),
        "run_count": len(runs),
        "analyses": rows,
        "runs": runs,
        "reference": reference,
        "input_sha256": {
            "analysis_plan": sha256(analysis_plan_path),
            "peak_policy": sha256(policy_path),
            "peak_plan": sha256(peak_plan_path),
            "peak_plan_summary": sha256(peak_summary_path),
            "filtering_output_manifest": sha256(manifest_path),
            "filtering_job_status": sha256(status_path),
            "filtering_input_provenance": hashes["outputs/input_provenance.json"],
            "filtering_policy": hashes["outputs/filtering_parameters.json"],
            "reference_fai": reference["fai_sha256"],
            "reference_provenance": reference["provenance_sha256"],
            "upstream_alignment_manifest": filtering_input_provenance["source_manifest_sha256"],
            "reference_mapping": reference_metadata["mapping_sha256"],
        },
        "pairing_source": str(Path(analysis_plan_path)),
        "pairing_policy": "exact_inheritance_no_downstream_repairing",
        "full_bam_and_csi_hashes_verified": bool(verify_large_files),
    }


def cohort_summary(records):
    ids = [record.get("analysis_id") for record in records]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate analysis record in cohort summary")
    return {
        "schema_version": 1,
        "scope": "dynamic_chipseq_peak_calling_cohort_summary",
        "analysis_count": len(records),
        "study_counts": dict(sorted(Counter(x["study_accession"] for x in records).items())),
        "peak_mode_counts": dict(sorted(Counter(x["peak_mode"] for x in records).items())),
        "target_counts": dict(sorted(Counter(x["declared_target"] for x in records).items())),
        "completion_status_counts": dict(sorted(
            Counter(x["peak_calling_status"] for x in records).items()
        )),
        "total_peak_count": sum(int(x["peak_count"]) for x in records),
        "analysis_ids": sorted(ids),
        "blacklist_filtering": "not_performed_no_validated_CHO_blacklist_available",
        "replicate_idr": "not_performed_no_validated_replicate_groups_in_authoritative_plan",
        "qc_interpretation": "descriptive_no_universal_pass_fail_thresholds",
    }


def _runtime_analysis(runtime, analysis_id):
    if (runtime.get("schema_version") != 1
            or runtime.get("scope") != "dynamic_chipseq_peak_calling_runtime"):
        raise ValueError("Unsupported dynamic peak runtime manifest")
    if runtime.get("full_bam_and_csi_hashes_verified") is not True:
        raise ValueError("Runtime BAM/CSI hashes were not fully verified")
    matches = [row for row in runtime.get("analyses", []) if row.get("analysis_id") == analysis_id]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one runtime analysis: {analysis_id}")
    row = matches[0]
    ip_run = row.get("ip_run_accession")
    if not isinstance(runtime.get("runs"), dict) or ip_run not in runtime["runs"]:
        raise ValueError("Runtime manifest lacks the authoritative IP run")
    return row, runtime["runs"][ip_run]


def analysis_parameters(runtime_path, analysis_id, fragment_path=None):
    runtime_path = Path(runtime_path)
    runtime = json.loads(runtime_path.read_text(encoding="utf-8"))
    row, ip = _runtime_analysis(runtime, analysis_id)
    reference = runtime.get("reference", {})
    fai_path = reference.get("source_path")
    if not fai_path:
        raise ValueError("Runtime manifest lacks the verified reference FAI path")
    fai = parse_fai(fai_path)
    if (fai["fai_sha256"] != reference.get("fai_sha256")
            or fai["nuclear_span"] != reference.get("nuclear_span")):
        raise ValueError("Runtime reference FAI identity or nuclear span changed")
    fragment = None
    fragment_hash = None
    raw_table_path = None
    if fragment_path is not None:
        fragment_path = Path(fragment_path)
        fragment = json.loads(fragment_path.read_text(encoding="utf-8"))
        fragment_hash = sha256(fragment_path)
        raw_table_path = fragment_path.with_name("phantompeakqualtools.tsv")
    contract = macs3_contract(
        row, fai["nuclear_span"], fragment, ip, fragment_path, raw_table_path
    )
    contract["runtime_manifest_sha256"] = sha256(runtime_path)
    contract["input_sha256"] = dict(runtime["input_sha256"])
    contract["fragment_result_sha256"] = fragment_hash
    contract["fragment_result_path"] = (
        str(fragment_path.resolve()) if fragment_path is not None else None
    )
    contract["blacklist_filtering"] = "not_performed_no_validated_CHO_blacklist_available"
    contract["replicate_idr"] = "not_performed_no_validated_replicate_groups_in_authoritative_plan"
    return contract


def run_macs3(
    runtime_path, parameters_path, analysis_id, output_dir, log_path, version_path,
    fragment_path=None,
):
    runtime = json.loads(Path(runtime_path).read_text(encoding="utf-8"))
    parameters = json.loads(Path(parameters_path).read_text(encoding="utf-8"))
    row, _ = _runtime_analysis(runtime, analysis_id)
    reconstructed = analysis_parameters(runtime_path, analysis_id, fragment_path)
    if parameters != reconstructed:
        raise ValueError("MACS3 parameter record differs from reconstructed verified contract")
    parameters = reconstructed
    ip = runtime["runs"][row["ip_run_accession"]]
    control = runtime["runs"][row["control_run_accession"]]
    for label, item in (("IP BAM", ip), ("IP CSI", ip), ("Input BAM", control), ("Input CSI", control)):
        is_index = "CSI" in label
        verify_hash(
            item["csi_path"] if is_index else item["bam_path"],
            item["csi_sha256"] if is_index else item["bam_sha256"],
            label,
        )
    output_dir = Path(output_dir)
    expected = [
        output_dir / f"{analysis_id}_peaks.xls",
        output_dir / f"{analysis_id}_treat_pileup.bdg",
        output_dir / f"{analysis_id}_control_lambda.bdg",
    ]
    if parameters["peak_mode"] == "narrow":
        expected.extend([
            output_dir / f"{analysis_id}_peaks.narrowPeak",
            output_dir / f"{analysis_id}_summits.bed",
        ])
    else:
        expected.extend([
            output_dir / f"{analysis_id}_peaks.broadPeak",
            output_dir / f"{analysis_id}_peaks.gappedPeak",
        ])
    existing = [str(path) for path in expected if path.exists()]
    if existing:
        raise ValueError(f"Refusing unsafe MACS3 overwrite: {existing}")
    output_dir.mkdir(parents=True, exist_ok=True)
    command = [
        "macs3", "callpeak", "-t", ip["bam_path"], "-c", control["bam_path"],
        "-f", "BAM", "-g", str(parameters["macs3_genome_size"]),
        "-n", analysis_id, "--outdir", str(output_dir), "-q", str(parameters["qvalue"]),
        "--scale-to", parameters["scale_to"], *parameters["flags"],
    ]
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    version = subprocess.check_output(["macs3", "--version"], text=True).strip()
    Path(version_path).write_text(version + "\n", encoding="utf-8")
    missing = [str(path) for path in expected if not path.is_file()]
    if missing:
        raise ValueError(f"MACS3 did not produce required mode-specific outputs: {missing}")


def integer_file(path):
    text = Path(path).read_text(encoding="utf-8").strip()
    if not re.fullmatch(r"[0-9]+", text):
        raise ValueError(f"Expected non-negative integer count in {path}")
    return int(text)


def summarize_analysis(
    runtime_path, parameters_path, analysis_id, primary_peak, secondary_peak,
    xls, treat_bdg, control_bdg, log, version, ip_total_path, input_total_path,
    ip_overlap_path, input_overlap_path, qc_path, provenance_path, fragment_path=None,
):
    runtime = json.loads(Path(runtime_path).read_text(encoding="utf-8"))
    parameters = json.loads(Path(parameters_path).read_text(encoding="utf-8"))
    rows = [row for row in runtime["analyses"] if row["analysis_id"] == analysis_id]
    if len(rows) != 1:
        raise ValueError("Analysis identity mismatch during QC")
    if parameters != analysis_parameters(runtime_path, analysis_id, fragment_path):
        raise ValueError("QC parameter record differs from reconstructed verified contract")
    row = rows[0]
    mode = row["peak_mode"]
    primary_mode = "narrow" if mode == "narrow" else "broad"
    peak_stats = parse_peak_file(primary_peak, primary_mode, runtime["reference"]["contigs"])
    secondary_mode = None if mode == "narrow" else "gapped"
    if mode == "narrow":
        summit_count = 0
        with Path(secondary_peak).open(encoding="utf-8") as handle:
            for raw in handle:
                if raw.strip() and not raw.startswith(("#", "track", "browser")):
                    summit_count += 1
        if summit_count != peak_stats["peak_count"]:
            raise ValueError("Narrow peak/summit count disagreement")
        secondary_stats = {"summit_count": summit_count}
    else:
        secondary_stats = parse_peak_file(
            secondary_peak, secondary_mode, runtime["reference"]["contigs"]
        )
    ip = runtime["runs"][row["ip_run_accession"]]
    control = runtime["runs"][row["control_run_accession"]]
    ip_total, input_total = integer_file(ip_total_path), integer_file(input_total_path)
    if ip_total != ip["retained_reads"] or input_total != control["retained_reads"]:
        raise ValueError("FRiP denominator disagrees with filtering provenance")
    frip = descriptive_frip(
        ip_total, input_total, integer_file(ip_overlap_path), integer_file(input_overlap_path)
    )
    required_outputs = [primary_peak, secondary_peak, xls, treat_bdg, control_bdg, log, version]
    for path in required_outputs:
        if not Path(path).is_file():
            raise ValueError(f"Missing analysis output: {path}")
    qc = {
        "schema_version": 1,
        "analysis_id": analysis_id,
        "study_accession": row["study_accession"],
        "ip_run_accession": row["ip_run_accession"],
        "control_run_accession": row["control_run_accession"],
        "declared_target": row["declared_target"],
        "peak_mode": mode,
        "fragment_size_bp": parameters["fragment_size_bp"],
        "fragment_size_source": parameters["fragment_size_source"],
        **peak_stats,
        "secondary_output_statistics": secondary_stats,
        "nuclear_reference_span": runtime["reference"]["nuclear_span"],
        "peak_nuclear_span_coverage_fraction": (
            peak_stats["union_peak_bp"] / runtime["reference"]["nuclear_span"]
        ),
        **frip,
        "peak_calling_status": "completed_with_peaks" if peak_stats["peak_count"] else "completed_no_peaks",
        "blacklist_filtering": "not_performed_no_validated_CHO_blacklist_available",
        "replicate_idr": "not_performed_no_validated_replicate_groups_in_authoritative_plan",
    }
    provenance = {
        "schema_version": 1,
        "analysis_id": analysis_id,
        "pairing_source": runtime["pairing_source"],
        "pairing_policy": runtime["pairing_policy"],
        "input_sha256": dict(runtime["input_sha256"]),
        "ip_bam_sha256": ip["bam_sha256"],
        "ip_csi_sha256": ip["csi_sha256"],
        "input_bam_sha256": control["bam_sha256"],
        "input_csi_sha256": control["csi_sha256"],
        "fragment_result_sha256": parameters.get("fragment_result_sha256"),
        "parameters_sha256": sha256(parameters_path),
        "output_sha256": {Path(path).name: sha256(path) for path in required_outputs},
    }
    atomic_json(qc_path, qc)
    atomic_json(provenance_path, provenance)
    return qc


def publish_outputs(source, destination):
    source, destination = Path(source), Path(destination)
    if not source.is_dir():
        raise ValueError(f"Missing dynamic peak output directory: {source}")
    output = destination / "outputs"
    validate_output_destination(output)
    shutil.copytree(source, output)
    files = sorted(path for path in output.rglob("*") if path.is_file())
    if not files:
        raise ValueError("Dynamic peak output directory is empty")
    manifest = destination / "output.sha256"
    validate_output_destination(manifest)
    manifest.write_text("".join(
        f"{sha256(path)}  {path.relative_to(destination).as_posix()}\n" for path in files
    ), encoding="utf-8")
    for digest, relative in checksum_manifest(manifest).items():
        verify_hash(destination / relative, digest, relative)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    phantom = sub.add_parser("parse-phantom")
    phantom.add_argument("--table", required=True)
    phantom.add_argument("--runtime", required=True)
    phantom.add_argument("--analysis", required=True)
    phantom.add_argument("--output", required=True)
    fai = sub.add_parser("validate-fai")
    fai.add_argument("--fai", required=True)
    fai.add_argument("--output", required=True)
    parameters = sub.add_parser("analysis-parameters")
    parameters.add_argument("--runtime", required=True)
    parameters.add_argument("--analysis", required=True)
    parameters.add_argument("--fragment")
    parameters.add_argument("--output", required=True)
    call = sub.add_parser("call-macs3")
    call.add_argument("--runtime", required=True)
    call.add_argument("--parameters", required=True)
    call.add_argument("--analysis", required=True)
    call.add_argument("--fragment")
    call.add_argument("--output-dir", required=True)
    call.add_argument("--log", required=True)
    call.add_argument("--version", required=True)
    summary = sub.add_parser("summarize-analysis")
    for name in (
        "runtime", "parameters", "analysis", "primary-peak", "secondary-peak",
        "xls", "treat-bdg", "control-bdg", "log", "version", "ip-total",
        "input-total", "ip-overlap", "input-overlap", "qc", "provenance",
    ):
        summary.add_argument("--" + name, required=True)
    summary.add_argument("--fragment")
    cohort = sub.add_parser("summarize-cohort")
    cohort.add_argument("--output", required=True)
    cohort.add_argument("records", nargs="+")
    publish = sub.add_parser("publish")
    publish.add_argument("--source", required=True)
    publish.add_argument("--destination", required=True)
    args = parser.parse_args()
    if args.command == "parse-phantom":
        runtime = json.loads(Path(args.runtime).read_text(encoding="utf-8"))
        row, ip = _runtime_analysis(runtime, args.analysis)
        result = parse_phantompeak_table(
            args.table, row["ip_run_accession"], ip["bam_path"],
            args.analysis, ip["bam_sha256"],
        )
    elif args.command == "validate-fai":
        result = parse_fai(args.fai)
        result["contigs"] = dict(sorted(result["contigs"].items()))
    elif args.command == "analysis-parameters":
        result = analysis_parameters(args.runtime, args.analysis, args.fragment)
    elif args.command == "call-macs3":
        run_macs3(
            args.runtime, args.parameters, args.analysis, args.output_dir,
            args.log, args.version, args.fragment,
        )
        return
    elif args.command == "summarize-analysis":
        summarize_analysis(
            args.runtime, args.parameters, args.analysis, args.primary_peak,
            args.secondary_peak, args.xls, args.treat_bdg, args.control_bdg,
            args.log, args.version, args.ip_total, args.input_total,
            args.ip_overlap, args.input_overlap, args.qc, args.provenance, args.fragment,
        )
        return
    elif args.command == "summarize-cohort":
        records = [json.loads(Path(path).read_text(encoding="utf-8")) for path in args.records]
        result = cohort_summary(records)
    else:
        publish_outputs(args.source, args.destination)
        return
    atomic_json(args.output, result)


if __name__ == "__main__":
    main()
