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


def _job_status(path):
    path = Path(path)

    if not path.is_file():
        raise ValueError(
            f"Missing upstream alignment job status: {path}"
        )

    result = {}

    for line in path.read_text(
        encoding="utf-8"
    ).splitlines():

        if not line.strip():
            continue

        fields = line.split()

        if len(fields) != 2:
            raise ValueError(
                f"Malformed upstream alignment job status: {line!r}"
            )

        key, value = fields

        if key in result:
            raise ValueError(
                f"Duplicate upstream alignment job-status key: {key}"
            )

        result[key] = value

    return result


def configuration(root=Path(".")):
    cfg = load(
        Path(root)
        / CONFIG
    )

    expected_scope = (
        "dynamic_single_end_duplicate_marking_and_"
        "nuclear_filtering"
    )

    if (
        cfg["schema_version"] != 1
        or cfg["min_mapq"] != 30
        or cfg["exclude_flags"] != 3844
        or cfg["exclude_mapq_255"] is not True
        or cfg["optical_duplicate_detection"] is not False
        or cfg["duplicate_scoring_strategy"] != "SUM_OF_BASE_QUALITIES"
        or cfg["threads"] != 4
        or cfg["picard_heap_mb"] != 12000
        or cfg["mitochondrial_accession"] != "NC_007936.1"
        or cfg.get("scope") != expected_scope
    ):
        raise ValueError(
            "Unsupported filtering policy; review code and tests "
            "before changing it"
        )

    pointer = cfg.get(
        "alignment_latest_submission"
    )

    if not isinstance(
        pointer,
        str,
    ):
        raise ValueError(
            "Filtering policy has no alignment_latest_submission"
        )

    relative(
        pointer
    )

    return cfg


def prepare(root):
    root = Path(
        root
    ).resolve()

    cfg = configuration(
        root
    )

    pointer = (
        root
        / relative(
            cfg[
                "alignment_latest_submission"
            ]
        )
    )

    if not pointer.is_file():
        raise ValueError(
            "No submitted dynamic ChIP alignment is available: "
            f"{pointer}"
        )

    pointer_value = pointer.read_text(
        encoding="utf-8"
    ).strip()

    if not pointer_value:
        raise ValueError(
            "Alignment latest-submission pointer is empty"
        )

    submission = Path(
        pointer_value
    )

    if not submission.is_absolute():
        raise ValueError(
            "Alignment latest-submission pointer must be absolute"
        )

    submission = submission.resolve()

    expected_parent = pointer.parent.resolve()

    if (
        submission.parent != expected_parent
        or not re.fullmatch(
            r"submission_[A-Za-z0-9_]+",
            submission.name,
        )
    ):
        raise ValueError(
            "Alignment latest-submission pointer does not identify "
            "a valid submission directory"
        )

    if not submission.is_dir():
        raise ValueError(
            f"Alignment submission directory is missing: {submission}"
        )

    job_id_file = (
        submission
        / "job_id.txt"
    )

    if not job_id_file.is_file():
        raise ValueError(
            "Latest alignment submission has no job_id.txt and "
            "is therefore not a completed productive submission"
        )

    job_id = job_id_file.read_text(
        encoding="utf-8"
    ).strip()

    if not re.fullmatch(
        r"[0-9]+",
        job_id,
    ):
        raise ValueError(
            f"Invalid alignment job id: {job_id!r}"
        )

    source = (
        submission
        / f"job_{job_id}"
    )

    if not source.is_dir():
        raise ValueError(
            f"Alignment job directory is missing: {source}"
        )

    status_path = (
        source
        / "job_status.tsv"
    )

    status = _job_status(
        status_path
    )

    if (
        status.get("stage") != "completed"
        or status.get("exit_status") != "0"
    ):
        raise ValueError(
            "Upstream alignment did not complete successfully"
        )

    output_validation = (
        source
        / "output_copy_validation.log"
    )

    if not output_validation.is_file():
        raise ValueError(
            "Successful alignment publication log is missing"
        )

    if (
        "[OK] Outputs copied and SHA-256 verified"
        not in output_validation.read_text(
            encoding="utf-8"
        )
    ):
        raise ValueError(
            "Alignment output publication was not verified"
        )

    manifest = (
        source
        / "output.sha256"
    )

    if not manifest.is_file():
        raise ValueError(
            "Successful alignment output manifest is missing"
        )

    hashes = checksums(
        manifest
    )

    base_small = [
        "outputs/reference/reference_provenance.json",
        "outputs/reference/genome_plus_mt.fa.fai",
        "outputs/input_provenance.json",
        "outputs/alignment_parameters.json",
        "outputs/alignment_qc.tsv",
        "outputs/verified_fastq_inputs.json",
    ]

    for name in base_small:

        if name not in hashes:
            raise ValueError(
                f"Alignment manifest is missing required output: {name}"
            )

        path = (
            source
            / name
        )

        if not path.is_file():
            raise ValueError(
                f"Alignment output is missing: {path}"
            )

        check(
            path,
            hashes[name],
        )

    provenance = load(
        source
        / "outputs/reference/reference_provenance.json"
    )

    if (
        provenance.get(
            "mitochondrial_accession"
        )
        != cfg[
            "mitochondrial_accession"
        ]
    ):
        raise ValueError(
            "Unexpected mitochondrial accession in alignment provenance"
        )

    alignment_parameters = load(
        source
        / "outputs/alignment_parameters.json"
    )

    if (
        alignment_parameters.get(
            "diagnostic_mapq"
        )
        != cfg[
            "min_mapq"
        ]
        or alignment_parameters.get(
            "mitochondrial_accession"
        )
        != cfg[
            "mitochondrial_accession"
        ]
    ):
        raise ValueError(
            "Alignment/filtering MAPQ or mitochondrial policy mismatch"
        )

    upstream = load(
        source
        / "outputs/input_provenance.json"
    )

    if upstream.get(
        "schema_version"
    ) != 1:
        raise ValueError(
            "Unsupported alignment input-provenance schema"
        )

    if upstream.get(
        "library_layout"
    ) != "SINGLE":
        raise ValueError(
            "Dynamic filtering currently requires SINGLE alignment input"
        )

    if upstream.get(
        "instrument_platform"
    ) != "ILLUMINA":
        raise ValueError(
            "Dynamic filtering currently requires ILLUMINA alignment input"
        )

    upstream_runs = upstream.get(
        "runs"
    )

    if (
        not isinstance(
            upstream_runs,
            list,
        )
        or not upstream_runs
    ):
        raise ValueError(
            "Alignment input provenance contains no runs"
        )

    declared_run_count = upstream.get(
        "run_count"
    )

    if (
        declared_run_count is not None
        and int(
            declared_run_count
        )
        != len(
            upstream_runs
        )
    ):
        raise ValueError(
            "Alignment input-provenance run_count does not match runs"
        )

    observed_roles = Counter()

    seen = set()

    runs = []

    small = list(
        base_small
    )

    for row in upstream_runs:

        if not isinstance(
            row,
            dict,
        ):
            raise ValueError(
                "Alignment run provenance entry must be an object"
            )

        run = row.get(
            "run_accession"
        )

        role = row.get(
            "role"
        )

        if (
            not isinstance(
                run,
                str,
            )
            or not re.fullmatch(
                r"[DES]RR[0-9]+",
                run,
            )
        ):
            raise ValueError(
                f"Invalid alignment run accession: {run!r}"
            )

        if run in seen:
            raise ValueError(
                f"Duplicate alignment run: {run}"
            )

        seen.add(
            run
        )

        if role not in {
            "ip",
            "input",
        }:
            raise ValueError(
                f"Invalid alignment role for {run}: {role!r}"
            )

        observed_roles[
            role
        ] += 1

        try:
            input_reads = int(
                row[
                    "reads"
                ]
            )

        except (
            KeyError,
            TypeError,
            ValueError,
        ) as exc:
            raise ValueError(
                f"Invalid alignment input read count for {run}"
            ) from exc

        if input_reads <= 0:
            raise ValueError(
                f"Non-positive alignment input read count for {run}"
            )

        qc_name = (
            f"outputs/{run}/alignment_qc.json"
        )

        if qc_name not in hashes:
            raise ValueError(
                f"Alignment manifest is missing QC for {run}"
            )

        qc_path = (
            source
            / qc_name
        )

        if not qc_path.is_file():
            raise ValueError(
                f"Alignment QC output is missing: {qc_path}"
            )

        check(
            qc_path,
            hashes[
                qc_name
            ],
        )

        small.append(
            qc_name
        )

        qc = load(
            qc_path
        )

        try:
            mapped_reads = int(
                qc[
                    "mapped_reads"
                ]
            )

            nuclear_mapq = int(
                qc[
                    "nuclear_mapq_ge_threshold"
                ]
            )

        except (
            KeyError,
            TypeError,
            ValueError,
        ) as exc:
            raise ValueError(
                f"Invalid alignment QC counts for {run}"
            ) from exc

        if (
            qc.get(
                "run_accession"
            )
            != run
            or qc.get(
                "role"
            )
            != role
            or qc.get(
                "duplicate_marking"
            )
            != "not_performed"
            or qc.get(
                "bam_filtering"
            )
            != "not_performed"
            or int(
                qc.get(
                    "nonprimary_records",
                    -1,
                )
            )
            != 0
            or int(
                qc.get(
                    "primary_reads",
                    -1,
                )
            )
            != input_reads
            or int(
                qc.get(
                    "alignment_records",
                    -1,
                )
            )
            != input_reads
            or qc.get(
                "diagnostic_mapq"
            )
            != cfg[
                "min_mapq"
            ]
        ):
            raise ValueError(
                f"Unexpected upstream SINGLE alignment QC: {run}"
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

        bam_name = (
            f"outputs/{run}/raw.sorted.bam"
        )

        index_name = (
            f"outputs/{run}/raw.sorted.bam.csi"
        )

        for required in (
            bam_name,
            index_name,
        ):

            if required not in hashes:
                raise ValueError(
                    f"Alignment manifest is missing output: {required}"
                )

            candidate = (
                source
                / required
            )

            if (
                not candidate.is_file()
                or candidate.stat().st_size <= 0
            ):
                raise ValueError(
                    f"Missing upstream alignment output: {candidate}"
                )

        bam = (
            source
            / bam_name
        )

        runs.append({
            "run_accession": run,
            "role": role,
            "source_relative": bam_name,
            "sha256": hashes[
                bam_name
            ],
            "bytes": bam.stat().st_size,
            "input_reads": input_reads,
            "mapped_reads": mapped_reads,
            "nuclear_mapq_ge_threshold": nuclear_mapq,
        })

    if observed_roles["ip"] < 1:
        raise ValueError(
            "Alignment cohort contains no IP runs"
        )

    if observed_roles["input"] < 1:
        raise ValueError(
            "Alignment cohort contains no Input runs"
        )

    declared_roles = upstream.get(
        "role_counts"
    )

    if (
        declared_roles is not None
        and {
            "ip": int(
                declared_roles.get(
                    "ip",
                    -1,
                )
            ),
            "input": int(
                declared_roles.get(
                    "input",
                    -1,
                )
            ),
        }
        != {
            "ip": observed_roles["ip"],
            "input": observed_roles["input"],
        }
    ):
        raise ValueError(
            "Alignment input-provenance role_counts do not match runs"
        )

    return {
        "schema_version": 1,
        "source_root": str(
            source
        ),
        "alignment_submission": str(
            submission
        ),
        "alignment_job_id": job_id,
        "library_layout": "SINGLE",
        "instrument_platform": "ILLUMINA",
        "run_count": len(
            runs
        ),
        "role_counts": {
            "ip": observed_roles[
                "ip"
            ],
            "input": observed_roles[
                "input"
            ],
        },
        "runs": sorted(
            runs,
            key=lambda item: item[
                "run_accession"
            ],
        ),
        "configuration_sha256": sha(
            root
            / CONFIG
        ),
        "source_manifest_sha256": sha(
            manifest
        ),
        "job_status_sha256": sha(
            status_path
        ),
        "small_files": {
            name: hashes[
                name
            ]
            for name in small
        },
        "reference_provenance": provenance,
        "full_bam_sha256_verification_deferred_to_compute": True,
        "upstream_raw_bam_cleanup_authorized": False,
    }


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
