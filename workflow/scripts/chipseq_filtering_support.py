#!/usr/bin/env python3
"""Verify upstream BAMs, audit duplicate marking and filter SINGLE ChIP-seq runs."""
import argparse
import csv
import hashlib
import json
import re
import shutil
from collections import Counter
from pathlib import Path

CONFIG = "config/chipseq_filtering.json"
PLAN = "config/chipseq_filtering_inputs.json"
OUT = Path("results/chipseq/filtering")
REASONS = ("nonprimary", "unmapped", "qc_fail", "mitochondrial", "mapq_unknown", "low_mapq", "duplicate")


def sha(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def dump(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def load(path):
    return json.loads(Path(path).read_text())


def relative(value):
    path = Path(value)
    if path.is_absolute() or not path.parts or ".." in path.parts:
        raise ValueError(f"Expected a relative project path: {value}")
    return path


def checksums(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        digest, name = line.split("  ", 1)
        name = relative(name).as_posix()
        if not re.fullmatch(r"[0-9a-f]{64}", digest) or name in result:
            raise ValueError(f"Invalid checksum entry: {name}")
        result[name] = digest
    return result


def check(path, digest):
    if sha(path) != digest:
        raise ValueError(f"SHA-256 mismatch: {path}")


def configuration(root=Path(".")):
    cfg = load(root / CONFIG)
    if (cfg["schema_version"] != 1 or cfg["min_mapq"] != 30
            or cfg["exclude_flags"] != 3844 or cfg["exclude_mapq_255"] is not True
            or cfg["optical_duplicate_detection"] is not False
            or cfg["duplicate_scoring_strategy"] != "SUM_OF_BASE_QUALITIES"
            or cfg["threads"] != 4 or cfg["picard_heap_mb"] != 12000
            or cfg["mitochondrial_accession"] != "NC_007936.1"):
        raise ValueError("Unsupported pilot policy; review code and tests before changing it")
    return cfg


def prepare(root):
    root = Path(root).resolve()
    cfg = configuration(root)
    source = root / relative(cfg["source_job"])
    evidence = root / relative(cfg["alignment_evidence"])
    # Anchor the live results to the previously archived validation evidence.
    archived = checksums(evidence / "SHA256SUMS.txt")
    for name in ("job/output.sha256", "job/job_status.tsv"):
        check(evidence / name, archived[name])
    check(source / "output.sha256", sha(evidence / "job/output.sha256"))
    check(source / "job_status.tsv", sha(evidence / "job/job_status.tsv"))
    status = dict(line.split() for line in (source / "job_status.tsv").read_text().splitlines())
    if status.get("stage") != "completed" or status.get("exit_status") != "0":
        raise ValueError("Upstream alignment did not complete successfully")
    hashes = checksums(source / "output.sha256")
    pilot = load(root / "config/chipseq_pilot.json")
    roles = {pilot["input_run_accession"]: "input", pilot["ip_run_accession"]: "ip"}
    if len(roles) != 2 or any(not re.fullmatch(r"[DES]RR[0-9]+", run) for run in roles):
        raise ValueError("Expected one input and one distinct IP")
    small = ["outputs/reference/reference_provenance.json", "outputs/reference/genome_plus_mt.fa.fai",
             "outputs/input_provenance.json", "outputs/alignment_parameters.json", "outputs/alignment_qc.tsv"]
    small += [f"outputs/{run}/alignment_qc.json" for run in roles]
    for name in small:
        check(source / name, hashes[name])
    provenance = load(source / small[0])
    if provenance["mitochondrial_accession"] != cfg["mitochondrial_accession"]:
        raise ValueError("Unexpected mitochondrial accession")
    upstream = load(source / "outputs/input_provenance.json")
    if len(upstream["runs"]) != 2 or {r["run_accession"]: r["role"] for r in upstream["runs"]} != roles:
        raise ValueError("Upstream BAM selection differs from the current pilot")
    runs = []
    for run, role in sorted(roles.items()):
        qc = load(source / f"outputs/{run}/alignment_qc.json")
        row = next(r for r in upstream["runs"] if r["run_accession"] == run)
        if (qc["run_accession"] != run or qc["role"] != role
                or qc["duplicate_marking"] != "not_performed" or qc["bam_filtering"] != "not_performed"
                or qc["nonprimary_records"] != 0 or qc["primary_reads"] != row["reads"]
                or qc["alignment_records"] != row["reads"] or row["reads"] <= 0
                or qc["diagnostic_mapq"] != cfg["min_mapq"]):
            raise ValueError(f"Unexpected upstream SINGLE pilot QC: {run}")
        name = f"outputs/{run}/raw.sorted.bam"
        bam = source / name
        if not bam.is_file() or bam.stat().st_size <= 0:
            raise ValueError(f"Missing upstream BAM: {bam}")
        runs.append(dict(run_accession=run, role=role, source_relative=name,
                         sha256=hashes[name], bytes=bam.stat().st_size,
                         input_reads=row["reads"], mapped_reads=qc["mapped_reads"],
                         nuclear_mapq_ge_threshold=qc["nuclear_mapq_ge_threshold"]))
    return dict(schema_version=1, source_root=str(source), runs=runs,
                configuration_sha256=sha(root / CONFIG), source_manifest_sha256=sha(source / "output.sha256"),
                alignment_evidence_sha256=sha(evidence / "SHA256SUMS.txt"),
                small_files={name: hashes[name] for name in small},
                reference_provenance=provenance)


def item_for(run):
    plan = load(PLAN)
    if plan["configuration_sha256"] != sha(CONFIG):
        raise ValueError("Filtering configuration differs from the submission plan")
    matches = [r for r in plan["runs"] if r["run_accession"] == run]
    if len(matches) != 1:
        raise ValueError("Unknown or duplicate filtering run")
    return plan, matches[0]


def validate_header(bam, run, fai, mt):
    import pysam
    sequences = {}
    for line in Path(fai).read_text().splitlines():
        fields = line.split("\t")
        if fields[0] in sequences:
            raise ValueError("Duplicate contig in reference index")
        sequences[fields[0]] = int(fields[1])
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        header = handle.header.to_dict()
        if header.get("HD", {}).get("SO") != "coordinate":
            raise ValueError("BAM must be coordinate sorted")
        if dict(zip(handle.references, handle.lengths)) != sequences or mt not in sequences:
            raise ValueError("BAM sequence dictionary differs from the mapping reference")
        groups = header.get("RG", [])
        if len(groups) != 1 or any(groups[0].get(key) != run for key in ("ID", "SM", "LB")):
            raise ValueError("Unexpected BAM sample/library read group")


def stage(run, output, report):
    import pysam
    cfg = configuration()
    plan, item = item_for(run)
    source = Path(plan["source_root"]) / relative(item["source_relative"])
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(".part")
    shutil.copyfile(source, partial)
    check(partial, item["sha256"])
    if partial.stat().st_size != item["bytes"]:
        raise ValueError("BAM size changed")
    pysam.quickcheck(str(partial))
    validate_header(partial, run, "inputs/reference/genome_plus_mt.fa.fai", cfg["mitochondrial_accession"])
    partial.replace(target)
    dump(report, {"run_accession": run, "bam_sha256": item["sha256"], "bytes": item["bytes"],
                  "reference_dictionary_verified": True})
    print(f"[OK] BAM copied, SHA-256 and reference dictionary verified: {run}")


def picard_metrics(path, run):
    lines = Path(path).read_text().splitlines()
    for index, line in enumerate(lines):
        if line.startswith("LIBRARY\t"):
            rows = []
            for entry in lines[index + 1:]:
                if not entry.strip() or entry.startswith("#"):
                    break
                rows.append(entry)
            parsed = list(csv.DictReader([line] + rows, delimiter="\t"))
            if len(parsed) != 1 or parsed[0]["LIBRARY"] != run:
                raise ValueError("Expected Picard metrics for exactly one SINGLE library")
            row = parsed[0]
            if int(row["READ_PAIRS_EXAMINED"]) != 0 or int(row["READ_PAIR_DUPLICATES"]) != 0:
                raise ValueError("Unexpected paired-end Picard metrics")
            return row
    raise ValueError("Missing Picard duplication metrics header")


def failures(flag, reference, mapq, mt, threshold):
    if flag & 1:
        raise ValueError("Paired reads require a separate reviewed filtering policy")
    return dict(nonprimary=bool(flag & 0x900), unmapped=bool(flag & 4),
                qc_fail=bool(flag & 512), mitochondrial=not bool(flag & 4) and reference == mt,
                mapq_unknown=not bool(flag & 4) and mapq == 255,
                low_mapq=not bool(flag & 4) and mapq < threshold, duplicate=bool(flag & 1024))


def filter_bam(run, bam, metrics, output, report, flow):
    import pysam
    cfg = configuration()
    _, item = item_for(run)
    picard = picard_metrics(metrics, run)
    counts, independent = Counter(), Counter()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_suffix(".part")
    with pysam.AlignmentFile(bam, "rb") as source:
        with pysam.AlignmentFile(str(partial), "wb", template=source) as target:
            for record in source.fetch(until_eof=True):
                if not record.has_tag("RG") or record.get_tag("RG") != run:
                    raise ValueError("Missing or unexpected record read group")
                checks = failures(record.flag, record.reference_name, record.mapping_quality,
                                  cfg["mitochondrial_accession"], cfg["min_mapq"])
                counts["input_records"] += 1
                if not record.is_secondary and not record.is_supplementary:
                    counts["primary_reads"] += 1
                    if not record.is_unmapped:
                        counts["mapped_reads"] += 1
                        counts["mapped_primary_duplicates"] += int(record.is_duplicate)
                        if record.reference_name != cfg["mitochondrial_accession"] and record.mapping_quality >= cfg["min_mapq"]:
                            counts["nuclear_mapq_ge_threshold_before_filtering"] += 1
                independent.update(name for name, failed in checks.items() if failed)
                reason = next((name for name in REASONS if checks[name]), None)
                if reason:
                    counts["removed_" + reason] += 1
                else:
                    target.write(record)
                    counts["retained_reads"] += 1
    if (counts["input_records"] != item["input_reads"] or counts["primary_reads"] != item["input_reads"]
            or counts["mapped_reads"] != item["mapped_reads"]
            or counts["nuclear_mapq_ge_threshold_before_filtering"] != item["nuclear_mapq_ge_threshold"]):
        raise ValueError("Duplicate marking changed the expected upstream alignment counts")
    if counts["mapped_primary_duplicates"] != int(picard["UNPAIRED_READ_DUPLICATES"]):
        raise ValueError("Marked BAM and Picard duplicate counts disagree")
    if counts["retained_reads"] + sum(counts["removed_" + name] for name in REASONS) != counts["input_records"]:
        raise ValueError("Filtering accounting does not balance")
    if not counts["retained_reads"]:
        raise ValueError("No reads survive filtering; review the run before downstream analysis")
    pysam.quickcheck(str(partial))
    partial.replace(output)
    data = dict(run_accession=run, role=item["role"], min_mapq=cfg["min_mapq"], exclude_flags=cfg["exclude_flags"],
                input_reads=counts["input_records"], retained_reads=counts["retained_reads"],
                retained_pct_of_input=100 * counts["retained_reads"] / counts["input_records"],
                mapped_primary_duplicates=counts["mapped_primary_duplicates"],
                picard_unpaired_reads_examined=int(picard["UNPAIRED_READS_EXAMINED"]),
                picard_percent_duplication=float(picard["PERCENT_DUPLICATION"]),
                optical_duplicate_detection="not_performed", signal_quality="not_evaluated",
                blacklist_filtering="not_performed", independent_counts_may_overlap=True)
    for reason in REASONS:
        data["removed_" + reason] = counts["removed_" + reason]
        data["independent_" + reason] = independent[reason]
    dump(report, data)
    with Path(flow).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["run_accession", "step", "criterion", "reads_entering", "reads_removed", "reads_remaining"])
        remaining = counts["input_records"]
        for step, reason in enumerate(REASONS, 1):
            removed = counts["removed_" + reason]
            writer.writerow([run, step, reason, remaining, removed, remaining - removed])
            remaining -= removed
    print(f"[OK] Filtering counts reconcile with upstream BAM and Picard: {run}")


def verify(run, bam, report, output):
    import pysam
    cfg = configuration()
    count = 0
    # Independent full output scan: never rely only on BAM quickcheck.
    with pysam.AlignmentFile(bam, "rb") as handle:
        if not handle.has_index():
            raise ValueError("Filtered BAM index is missing")
        for record in handle.fetch(until_eof=True):
            if (record.flag & (1 | 3844) or record.mapping_quality < cfg["min_mapq"]
                    or record.mapping_quality == 255 or record.reference_name == cfg["mitochondrial_accession"]):
                raise ValueError("Filtered BAM contains an excluded record")
            count += 1
    qc = load(report)
    if qc["run_accession"] != run or count != qc["retained_reads"]:
        raise ValueError("Filtered BAM count differs from the filtering report")
    dump(output, {"run_accession": run, "full_read_verified": True, "retained_reads": count,
                  "bam_sha256": sha(bam), "index_sha256": sha(bam + ".csi")})


def summarize(output, reports):
    rows = [load(path) for path in reports]
    with Path(output).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def publish(save):
    plan = load(PLAN)
    for item in plan["runs"]:
        run = item["run_accession"]
        verified = load(OUT / run / "filtered_validation.json")
        check(OUT / run / "filtered.bam", verified["bam_sha256"])
        check(OUT / run / "filtered.bam.csi", verified["index_sha256"])
    for src, dst in ((CONFIG, "filtering_parameters.json"), (PLAN, "input_provenance.json")):
        shutil.copyfile(src, OUT / dst)
    shutil.copytree("inputs/reference", OUT / "reference", dirs_exist_ok=True)
    manifest = []
    for path in sorted(OUT.rglob("*")):
        if not path.is_file():
            continue
        name = Path("outputs") / path.relative_to(OUT)
        target = Path(save) / name
        target.parent.mkdir(parents=True, exist_ok=True)
        digest = sha(path)
        shutil.copyfile(path, target)
        check(target, digest)
        manifest.append(f"{digest}  {name.as_posix()}\n")
    (Path(save) / "output.sha256").write_text("".join(manifest))
    print("[OK] Filtered BAMs and reports copied and SHA-256 verified")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("stage")
    for key in ("run", "output", "report"):
        p.add_argument("--" + key, required=True)
    p = sub.add_parser("filter")
    for key in ("run", "bam", "metrics", "output", "report", "flow"):
        p.add_argument("--" + key, required=True)
    p = sub.add_parser("verify")
    for key in ("run", "bam", "report", "output"):
        p.add_argument("--" + key, required=True)
    p = sub.add_parser("summarize")
    p.add_argument("--output", required=True)
    p.add_argument("reports", nargs="+")
    p = sub.add_parser("publish")
    p.add_argument("--save", required=True)
    args = parser.parse_args()
    if args.action == "stage": stage(args.run, args.output, args.report)
    elif args.action == "filter": filter_bam(args.run, args.bam, args.metrics, args.output, args.report, args.flow)
    elif args.action == "verify": verify(args.run, args.bam, args.report, args.output)
    elif args.action == "summarize": summarize(args.output, args.reports)
    elif args.action == "publish": publish(args.save)


if __name__ == "__main__":
    main()
