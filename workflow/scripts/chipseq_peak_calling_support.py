#!/usr/bin/env python3
"""Support functions for pilot ChIP-seq peak calling and signal QC."""

import argparse
import csv
import hashlib
import json
import re
import shutil
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

CONFIG = Path("config/chipseq_peak_calling.json")
PLAN = Path("config/chipseq_peak_calling_inputs.json")
OUT = Path("results/chipseq/peak_calling")


def load(path):
    with Path(path).open() as handle:
        return json.load(handle)


def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check(path, expected):
    path = Path(path)
    observed = sha(path)
    if observed != expected:
        raise ValueError(
            f"SHA-256 mismatch for {path}: expected {expected}, observed {observed}"
        )
    return observed


def verify_manifest(directory):
    directory = Path(directory)
    manifest = directory / "SHA256SUMS.txt"

    if not manifest.is_file():
        raise ValueError(f"Missing validation manifest: {manifest}")

    checked = 0

    for raw in manifest.read_text().splitlines():
        if not raw.strip():
            continue

        fields = raw.split(None, 1)
        if len(fields) != 2:
            raise ValueError(f"Malformed checksum line: {raw}")

        expected, relative = fields
        relative = relative.strip()

        if relative.startswith("./"):
            relative = relative[2:]

        target = directory / relative

        if not target.is_file():
            raise ValueError(f"Missing snapshot file: {target}")

        check(target, expected)
        checked += 1

    if checked == 0:
        raise ValueError("Empty snapshot checksum manifest")

    return checked


def configuration(root=ROOT):
    cfg = load(Path(root) / CONFIG)

    required = [
        "source_job",
        "filtering_evidence",
        "treatment_run",
        "control_run",
        "analysis_name",
        "format",
        "peak_mode",
        "effective_genome_size",
        "mitochondrial_accession",
        "qvalue",
        "keep_dup",
        "scale_to",
        "store_bdg",
        "spmr",
        "model",
    ]

    missing = [key for key in required if key not in cfg]
    if missing:
        raise ValueError(f"Missing configuration keys: {missing}")

    if cfg["treatment_run"] == cfg["control_run"]:
        raise ValueError("Treatment and control must be different runs")

    if cfg["format"] != "BAM":
        raise ValueError("Pilot is SINGLE-END and must use MACS3 format BAM")

    if cfg["peak_mode"] != "narrow":
        raise ValueError("Pilot H3K4me3 implementation expects narrow peak calling")

    if cfg["keep_dup"] != "all":
        raise ValueError("Filtered BAMs require MACS3 --keep-dup all")

    if cfg["scale_to"] not in {"small", "large"}:
        raise ValueError("scale_to must be 'small' or 'large'")

    if not 0 < float(cfg["qvalue"]) < 1:
        raise ValueError("qvalue must be between 0 and 1")

    if int(cfg["effective_genome_size"]) <= 0:
        raise ValueError("effective_genome_size must be positive")

    if cfg["spmr"] and not cfg["store_bdg"]:
        raise ValueError("SPMR requires bedGraph output")

    if cfg["model"] not in {"auto", "fixed"}:
        raise ValueError("model must be 'auto' or 'fixed'")

    if cfg["model"] == "fixed":
        if int(cfg.get("extsize", 0)) <= 0:
            raise ValueError("Fixed model requires a positive extsize")
        if not cfg.get("fragment_size_basis"):
            raise ValueError("Fixed model requires fragment_size_basis")

    return cfg


def fai_nuclear_span(path, mitochondrial_accession):
    total = 0
    nuclear = 0
    mt_seen = False
    sequences = 0

    with Path(path).open() as handle:
        for line in handle:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Malformed FAI line: {line.rstrip()}")

            name = fields[0]
            length = int(fields[1])

            total += length
            sequences += 1

            if name == mitochondrial_accession:
                mt_seen = True
            else:
                nuclear += length

    if not mt_seen:
        raise ValueError(
            f"Mitochondrial accession {mitochondrial_accession} absent from FAI"
        )

    return {
        "total_span": total,
        "nuclear_span": nuclear,
        "sequences": sequences,
    }


def prepare(root=ROOT):
    root = Path(root)
    cfg = configuration(root)

    evidence = root / cfg["filtering_evidence"]
    source = root / cfg["source_job"]

    if not evidence.is_dir():
        raise ValueError(f"Missing filtering evidence: {evidence}")

    snapshot_files_verified = verify_manifest(evidence)

    job_status = source / "job_status.tsv"
    if not job_status.is_file():
        raise ValueError(f"Missing filtering job status: {job_status}")

    status_text = job_status.read_text()

    if "stage\tcompleted" not in status_text or "exit_status\t0" not in status_text:
        raise ValueError("Filtering source job is not recorded as successfully completed")

    pilot = load(root / "config/chipseq_pilot.json")

    if pilot.get("ip_run_accession") != cfg["treatment_run"]:
        raise ValueError("Treatment run disagrees with chipseq_pilot.json")

    if pilot.get("input_run_accession") != cfg["control_run"]:
        raise ValueError("Control run disagrees with chipseq_pilot.json")

    roles = {
        cfg["treatment_run"]: "ip",
        cfg["control_run"]: "input",
    }

    runs = []

    for run, expected_role in roles.items():
        qc_path = evidence / run / "filtering_qc.json"
        validation_path = evidence / run / "filtered_validation.json"

        qc = load(qc_path)
        validation = load(validation_path)

        if qc["run_accession"] != run:
            raise ValueError(f"QC accession mismatch for {run}")

        if qc["role"] != expected_role:
            raise ValueError(
                f"Role mismatch for {run}: expected {expected_role}, got {qc['role']}"
            )

        if validation["run_accession"] != run:
            raise ValueError(f"Validation accession mismatch for {run}")

        if validation.get("full_read_verified") is not True:
            raise ValueError(f"Filtered BAM full-read validation missing for {run}")

        if int(validation["retained_reads"]) != int(qc["retained_reads"]):
            raise ValueError(f"Retained-read disagreement for {run}")

        bam_relative = f"outputs/{run}/filtered.bam"
        index_relative = f"outputs/{run}/filtered.bam.csi"

        bam = source / bam_relative
        index = source / index_relative

        if not bam.is_file() or bam.stat().st_size == 0:
            raise ValueError(f"Missing filtered BAM: {bam}")

        if not index.is_file() or index.stat().st_size == 0:
            raise ValueError(f"Missing CSI index: {index}")

        runs.append({
            "run_accession": run,
            "role": expected_role,
            "source_relative": bam_relative,
            "index_relative": index_relative,
            "bam_sha256": validation["bam_sha256"],
            "index_sha256": validation["index_sha256"],
            "retained_reads": int(qc["retained_reads"]),
            "retained_pct_of_input": float(qc["retained_pct_of_input"]),
            "picard_percent_duplication": float(qc["picard_percent_duplication"]),
        })

    reference = source / "outputs/reference/genome_plus_mt.fa.fai"

    if not reference.is_file():
        raise ValueError(f"Missing source reference FAI: {reference}")

    span = fai_nuclear_span(reference, cfg["mitochondrial_accession"])

    if span["nuclear_span"] != int(cfg["effective_genome_size"]):
        raise ValueError(
            "Configured effective genome size disagrees with nuclear reference span: "
            f"{cfg['effective_genome_size']} vs {span['nuclear_span']}"
        )

    return {
        "schema_version": 1,
        "configuration_sha256": sha(root / CONFIG),
        "filtering_evidence": cfg["filtering_evidence"],
        "filtering_evidence_files_verified": snapshot_files_verified,
        "source_root": str(source.resolve()),
        "analysis_name": cfg["analysis_name"],
        "effective_genome_size": int(cfg["effective_genome_size"]),
        "effective_genome_size_basis": cfg["effective_genome_size_basis"],
        "runs": runs,
        "reference": {
            "source_path": str(reference.resolve()),
            "sha256": sha(reference),
            "total_span": span["total_span"],
            "nuclear_span": span["nuclear_span"],
            "sequence_count": span["sequences"],
            "mitochondrial_accession": cfg["mitochondrial_accession"],
        },
    }


def stage(run, plan_path, bam_output, index_output, report):
    plan = load(plan_path)

    check(CONFIG, plan["configuration_sha256"])

    matches = [item for item in plan["runs"] if item["run_accession"] == run]

    if len(matches) != 1:
        raise ValueError(f"Expected exactly one input plan entry for {run}")

    item = matches[0]
    source_root = Path(plan["source_root"])

    bam_source = source_root / item["source_relative"]
    index_source = source_root / item["index_relative"]

    bam_output = Path(bam_output)
    index_output = Path(index_output)

    bam_output.parent.mkdir(parents=True, exist_ok=True)
    index_output.parent.mkdir(parents=True, exist_ok=True)

    try:
        shutil.copyfile(bam_source, bam_output)
        check(bam_output, item["bam_sha256"])

        shutil.copyfile(index_source, index_output)
        check(index_output, item["index_sha256"])

    except Exception:
        bam_output.unlink(missing_ok=True)
        index_output.unlink(missing_ok=True)
        raise

    dump(report, {
        "run_accession": run,
        "role": item["role"],
        "source_bam": str(bam_source),
        "source_index": str(index_source),
        "bam_sha256": item["bam_sha256"],
        "index_sha256": item["index_sha256"],
        "retained_reads": item["retained_reads"],
        "full_input_hash_verified": True,
    })


def parse_narrowpeak(path, mitochondrial_accession):
    intervals = defaultdict(list)
    widths = []
    enrichments = []
    q_scores = []
    chromosomes = set()
    peak_count = 0

    path = Path(path)

    if not path.is_file():
        raise ValueError(f"Missing narrowPeak: {path}")

    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                raise ValueError(
                    f"narrowPeak line {line_number} has {len(fields)} columns, expected >=10"
                )

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            signal = float(fields[6])
            qvalue_score = float(fields[8])
            summit_offset = int(fields[9])

            if start < 0 or end <= start:
                raise ValueError(f"Invalid interval at narrowPeak line {line_number}")

            width = end - start

            if summit_offset < 0 or summit_offset >= width:
                raise ValueError(
                    f"Summit outside peak at narrowPeak line {line_number}"
                )

            if chrom == mitochondrial_accession:
                raise ValueError("Mitochondrial peak detected after nuclear filtering")

            if signal < 0 or qvalue_score < 0:
                raise ValueError(f"Negative peak metric at line {line_number}")

            intervals[chrom].append((start, end))
            chromosomes.add(chrom)
            widths.append(width)
            enrichments.append(signal)
            q_scores.append(qvalue_score)
            peak_count += 1

    union_bp = 0

    for chrom_intervals in intervals.values():
        chrom_intervals.sort()

        current_start = None
        current_end = None

        for start, end in chrom_intervals:
            if current_start is None:
                current_start, current_end = start, end
                continue

            if start <= current_end:
                current_end = max(current_end, end)
            else:
                union_bp += current_end - current_start
                current_start, current_end = start, end

        if current_start is not None:
            union_bp += current_end - current_start

    return {
        "peak_count": peak_count,
        "chromosomes_with_peaks": len(chromosomes),
        "sum_peak_bp": sum(widths),
        "union_peak_bp": union_bp,
        "mean_peak_width_bp": statistics.mean(widths) if widths else None,
        "median_peak_width_bp": statistics.median(widths) if widths else None,
        "median_fold_enrichment": statistics.median(enrichments)
        if enrichments else None,
        "median_minus_log10_qvalue": statistics.median(q_scores)
        if q_scores else None,
    }


def count_bed_records(path):
    path = Path(path)

    if not path.is_file():
        raise ValueError(f"Missing BED file: {path}")

    count = 0

    with path.open() as handle:
        for line in handle:
            if line.strip() and not line.startswith("#") and not line.startswith("track"):
                count += 1

    return count


def parse_fragment_size(*paths):
    patterns = [
        re.compile(r"predicted fragment length is\s+(\d+)", re.I),
        re.compile(r"#\s*d\s*=\s*(\d+)", re.I),
        re.compile(r"fragment size\s*=\s*(\d+)", re.I),
    ]

    for path in paths:
        text = Path(path).read_text(errors="replace")

        for pattern in patterns:
            match = pattern.search(text)
            if match:
                return int(match.group(1))

    return None


def integer_file(path):
    text = Path(path).read_text().strip()

    if not re.fullmatch(r"[0-9]+", text):
        raise ValueError(f"Expected integer count in {path}: {text!r}")

    return int(text)


def summarize(
    peaks,
    summits,
    xls,
    treat_bdg,
    control_bdg,
    log,
    version,
    ip_total_file,
    control_total_file,
    ip_overlap_file,
    control_overlap_file,
    qc_json,
    qc_tsv,
    parameters_json,
    provenance_json,
):
    cfg = configuration(Path.cwd())
    plan = load(PLAN)

    check(CONFIG, plan["configuration_sha256"])

    for required in [xls, treat_bdg, control_bdg, log, version]:
        if not Path(required).is_file():
            raise ValueError(f"Missing MACS3 output: {required}")

    peak_stats = parse_narrowpeak(peaks, cfg["mitochondrial_accession"])
    summit_count = count_bed_records(summits)

    if summit_count != peak_stats["peak_count"]:
        raise ValueError(
            f"Peak/summit count disagreement: {peak_stats['peak_count']} vs {summit_count}"
        )

    ip_total = integer_file(ip_total_file)
    control_total = integer_file(control_total_file)
    ip_overlap = integer_file(ip_overlap_file)
    control_overlap = integer_file(control_overlap_file)

    by_role = {item["role"]: item for item in plan["runs"]}

    if ip_total != int(by_role["ip"]["retained_reads"]):
        raise ValueError(
            f"IP denominator disagreement: {ip_total} vs {by_role['ip']['retained_reads']}"
        )

    if control_total != int(by_role["input"]["retained_reads"]):
        raise ValueError(
            "Input denominator disagreement: "
            f"{control_total} vs {by_role['input']['retained_reads']}"
        )

    if ip_overlap > ip_total:
        raise ValueError("IP reads in peaks exceed total IP reads")

    if control_overlap > control_total:
        raise ValueError("Input reads in peaks exceed total Input reads")

    ip_frip = ip_overlap / ip_total if ip_total else 0.0
    input_overlap_fraction = (
        control_overlap / control_total if control_total else 0.0
    )

    enrichment_ratio = (
        ip_frip / input_overlap_fraction
        if input_overlap_fraction > 0
        else None
    )

    if cfg["model"] == "fixed":
        fragment_size = int(cfg["extsize"])
        fragment_size_source = "configured_source_informed_extsize"
    else:
        fragment_size = parse_fragment_size(xls, log)
        fragment_size_source = "macs3_cross_correlation_model"

    coverage_fraction = (
        peak_stats["union_peak_bp"] / int(cfg["effective_genome_size"])
        if cfg["effective_genome_size"]
        else None
    )

    qc = {
        "analysis_name": cfg["analysis_name"],
        "histone_mark": cfg["histone_mark"],
        "condition": cfg["condition"],
        "treatment_run": cfg["treatment_run"],
        "control_run": cfg["control_run"],
        "peak_mode": cfg["peak_mode"],
        "macs3_format": cfg["format"],
        "qvalue_cutoff": float(cfg["qvalue"]),
        "effective_genome_size": int(cfg["effective_genome_size"]),
        "effective_genome_size_basis": cfg["effective_genome_size_basis"],
        "duplicate_policy": "duplicates_removed_upstream_keep_dup_all",
        "macs3_model": cfg["model"],
        "fragment_size_bp": fragment_size,
        "fragment_size_source": fragment_size_source,
        "peak_count": peak_stats["peak_count"],
        "chromosomes_with_peaks": peak_stats["chromosomes_with_peaks"],
        "sum_peak_bp": peak_stats["sum_peak_bp"],
        "union_peak_bp": peak_stats["union_peak_bp"],
        "peak_genome_coverage_fraction": coverage_fraction,
        "mean_peak_width_bp": peak_stats["mean_peak_width_bp"],
        "median_peak_width_bp": peak_stats["median_peak_width_bp"],
        "median_fold_enrichment": peak_stats["median_fold_enrichment"],
        "median_minus_log10_qvalue": peak_stats["median_minus_log10_qvalue"],
        "ip_total_reads": ip_total,
        "ip_reads_in_peaks": ip_overlap,
        "ip_frip": ip_frip,
        "input_total_reads": control_total,
        "input_reads_overlapping_ip_peaks": control_overlap,
        "input_peak_overlap_fraction": input_overlap_fraction,
        "ip_to_input_peak_overlap_ratio": enrichment_ratio,
        "signal_quality": "evaluated_descriptively",
        "peak_calling_status": (
            "completed_with_peaks"
            if peak_stats["peak_count"] > 0
            else "completed_no_peaks"
        ),
        "blacklist_filtering": "not_performed",
        "replicate_concordance": "not_evaluated_single_pilot_pair",
    }

    dump(qc_json, qc)

    columns = list(qc.keys())

    with Path(qc_tsv).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerow({
            key: "" if value is None else value
            for key, value in qc.items()
        })

    parameters = dict(cfg)
    parameters["macs3_version"] = Path(version).read_text().strip()
    parameters["fragment_size_bp_observed"] = fragment_size
    dump(parameters_json, parameters)

    dump(provenance_json, {
        "configuration_sha256": plan["configuration_sha256"],
        "source_root": plan["source_root"],
        "filtering_evidence": plan["filtering_evidence"],
        "runs": plan["runs"],
        "reference": plan["reference"],
    })


def publish(save):
    save = Path(save)

    source = OUT

    if not source.is_dir():
        raise ValueError(f"Missing peak-calling output directory: {source}")

    destination = save / "outputs"

    if destination.exists():
        raise ValueError(f"Refusing to overwrite existing outputs: {destination}")

    shutil.copytree(source, destination)

    files = sorted(path for path in destination.rglob("*") if path.is_file())

    if not files:
        raise ValueError("No outputs were copied")

    manifest = save / "output.sha256"

    manifest.write_text(
        "".join(
            f"{sha(path)}  {path.relative_to(save).as_posix()}\n"
            for path in files
        )
    )

    for line in manifest.read_text().splitlines():
        digest, relative = line.split(None, 1)
        check(save / relative.strip(), digest)

    print(
        f"[OK] Peak-calling outputs copied and SHA-256 verified: "
        f"{len(files)} files"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    stage_parser = sub.add_parser("stage")
    stage_parser.add_argument("--run", required=True)
    stage_parser.add_argument("--plan", default=str(PLAN))
    stage_parser.add_argument("--bam", required=True)
    stage_parser.add_argument("--index", required=True)
    stage_parser.add_argument("--report", required=True)

    summary = sub.add_parser("summarize")
    summary.add_argument("--peaks", required=True)
    summary.add_argument("--summits", required=True)
    summary.add_argument("--xls", required=True)
    summary.add_argument("--treat-bdg", required=True)
    summary.add_argument("--control-bdg", required=True)
    summary.add_argument("--log", required=True)
    summary.add_argument("--version", required=True)
    summary.add_argument("--ip-total", required=True)
    summary.add_argument("--control-total", required=True)
    summary.add_argument("--ip-overlap", required=True)
    summary.add_argument("--control-overlap", required=True)
    summary.add_argument("--qc-json", required=True)
    summary.add_argument("--qc-tsv", required=True)
    summary.add_argument("--parameters-json", required=True)
    summary.add_argument("--provenance-json", required=True)

    publish_parser = sub.add_parser("publish")
    publish_parser.add_argument("--save", required=True)

    args = parser.parse_args()

    if args.command == "stage":
        stage(
            args.run,
            args.plan,
            args.bam,
            args.index,
            args.report,
        )

    elif args.command == "summarize":
        summarize(
            args.peaks,
            args.summits,
            args.xls,
            args.treat_bdg,
            args.control_bdg,
            args.log,
            args.version,
            args.ip_total,
            args.control_total,
            args.ip_overlap,
            args.control_overlap,
            args.qc_json,
            args.qc_tsv,
            args.parameters_json,
            args.provenance_json,
        )

    elif args.command == "publish":
        publish(args.save)


if __name__ == "__main__":
    main()
