#!/usr/bin/env python3
"""Stage verified pilot FASTQs, construct mapping reference and audit alignments."""
import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dump(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def table(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def relative_path(value):
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ValueError(f"Expected a relative project path: {value}")
    return path


def checksums(path):
    records = {}
    for line in Path(path).read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        name = str(relative_path(name.lstrip("*")))
        if not re.fullmatch(r"[a-f0-9]{64}", digest) or name in records:
            raise ValueError(f"Invalid checksum entry: {name}")
        records[name] = digest
    return records


def check(path, expected):
    if sha(path) != expected:
        raise ValueError(f"SHA-256 mismatch: {path}")


def prepare(root):
    root = Path(root).resolve()
    config = json.loads((root / "config/chipseq_alignment.json").read_text())
    if config["schema_version"] != 1:
        raise ValueError("Unsupported schema")
    source = root / relative_path(config["source_project"])
    hashfile = root / relative_path(config["source_checksums"])
    statusfile = root / relative_path(config["source_job_status"])
    status = dict(line.split() for line in statusfile.read_text().splitlines() if line.strip())
    if status.get("exit_status") != "0" or status.get("workflow_exit_status") != "0":
        raise ValueError("Upstream preprocessing job did not complete successfully")
    hashes = checksums(hashfile)
    report = "results/preprocessing/fastp/reports/preprocessing_qc_by_fastq.tsv"
    job_report = "results/preprocessing/fastp/reports/preprocessing_qc_by_job.tsv"
    for name in (report, job_report):
        check(source / name, hashes[name])
    rows, jobs = table(source / report), table(source / job_report)
    pilot = json.loads((root / "config/chipseq_pilot.json").read_text())
    roles = {pilot["ip_run_accession"]: "ip", pilot["input_run_accession"]: "input"}
    for collection in (rows, jobs):
        if len(collection) != 2 or {r["run_accession"] for r in collection} != set(roles):
            raise ValueError("Upstream QC tables do not match the current pilot")
    plan = []
    for row in sorted(rows, key=lambda r: r["run_accession"]):
        run = row["run_accession"]
        if not re.fullmatch(r"[DES]RR[0-9]+", run) or row["fastq_role"] != "SINGLE":
            raise ValueError("This pilot workflow currently supports SINGLE runs only")
        job = next(r for r in jobs if r["run_accession"] == run)
        if row["study_accession"] != pilot["study_accession"] or row["omics"] != "ChIP-seq":
            raise ValueError("Upstream study or modality mismatch")
        processed = str(relative_path(row["processed_fastq"]))
        fastq = source / processed
        if not fastq.is_file() or fastq.stat().st_size <= 0:
            raise ValueError(f"Missing processed FASTQ: {fastq}")
        reads = int(job["reads_after"])
        if reads <= 0 or reads != int(row["post_total_sequences"]):
            raise ValueError("Inconsistent processed read count")
        plan.append(dict(run_accession=run, role=roles[run], source_relative=processed,
                         sha256=hashes[processed], bytes=fastq.stat().st_size, reads=reads,
                         destination=f"inputs/{run}.fastq.gz"))
    return dict(schema_version=1, source_root=str(source.resolve()), runs=plan,
                upstream_checksums_sha256=sha(hashfile), upstream_status_sha256=sha(statusfile),
                reports_sha256={name: hashes[name] for name in (report, job_report)})


def stage():
    plan = json.loads(Path("config/chipseq_alignment_inputs.json").read_text())
    verified = []
    for row in plan["runs"]:
        source = Path(plan["source_root"]) / relative_path(row["source_relative"])
        target = relative_path(row["destination"])
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_suffix(".part")
        shutil.copyfile(source, temporary)
        check(temporary, row["sha256"])
        if temporary.stat().st_size != row["bytes"]:
            raise ValueError(f"Unexpected FASTQ size: {source}")
        temporary.replace(target)
        verified.append({"run_accession": row["run_accession"], "sha256": row["sha256"]})
        print(f"[OK] Processed FASTQ copied and verified: {row['run_accession']}")
    dump("inputs/verified.json", {"runs": verified})


def mapping_reference():
    import yaml
    reference = yaml.safe_load(Path("config/config.yaml").read_text())["reference"]
    config = json.loads(Path("config/chipseq_alignment.json").read_text())
    nuclear = Path(reference["fasta"])
    for name, digest in checksums(reference["sha256"]).items():
        check(name, digest)
    accession = config["mitochondrial_accession"]
    if not re.fullmatch(r"[A-Z]+_[0-9]+\.[0-9]+", accession):
        raise ValueError("A versioned mitochondrial accession is required")
    params = urllib.parse.urlencode(dict(db="nuccore", id=accession, rettype="fasta", retmode="text"))
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + params
    with urllib.request.urlopen(url, timeout=90) as response:
        mt_text = response.read().decode("utf-8")
    lines = mt_text.splitlines()
    if not lines or not lines[0].startswith(">") or lines[0][1:].split()[0] != accession:
        raise ValueError("NCBI mitochondrial accession mismatch")
    if any(line.startswith(">") for line in lines[1:]):
        raise ValueError("Expected one mitochondrial record")
    sequence = "".join(lines[1:]).upper()
    if not sequence or not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", sequence):
        raise ValueError("Invalid mitochondrial sequence")
    identifiers = set()
    with nuclear.open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in identifiers or name == accession:
                    raise ValueError(f"Duplicate reference identifier: {name}")
                identifiers.add(name)
    if not identifiers:
        raise ValueError("Empty nuclear FASTA")
    directory = Path(reference["dir"]) / "chipseq"
    directory.mkdir(parents=True, exist_ok=True)
    mt = directory / "mitochondrial.fa"
    mt.write_text(mt_text.rstrip() + "\n")
    combined = directory / "genome_plus_mt.fa"
    temporary = combined.with_suffix(".part")
    with temporary.open("wb") as target, nuclear.open("rb") as source:
        shutil.copyfileobj(source, target, length=1024 * 1024)
        source.seek(-1, 2)
        if source.read(1) != b"\n":
            target.write(b"\n")
        target.write(mt.read_bytes())
    temporary.replace(combined)
    dump(directory / "reference_provenance.json", dict(
        nuclear_accession=reference["refseq_accession"], nuclear_sha256=sha(nuclear),
        mitochondrial_accession=accession, mitochondrial_sha256=sha(mt),
        mapping_sha256=sha(combined), nuclear_sequences=len(identifiers),
        mitochondrial_length=len(sequence), mitochondrial_url=url,
        configured_annotation_release=reference["annotation_release"],
        annotation_release_independently_verified=False,
        purpose="ChIP-seq nuclear and mitochondrial mapping"))
    print("[OK] Mapping reference built; checksums and unique identifiers verified")


def scan_sam(lines, mt, threshold):
    counts = Counter()
    for line in lines:
        if line.startswith("@"):
            continue
        fields = line.rstrip().split("\t")
        if len(fields) < 11:
            raise ValueError("Truncated SAM record")
        flag = int(fields[1])
        counts["alignment_records"] += 1
        if flag & 0x900:
            counts["nonprimary_records"] += 1
            continue
        counts["primary_reads"] += 1
        if flag & 1:
            raise ValueError("Unexpected paired record in SINGLE pilot")
        if flag & 4:
            counts["unmapped_reads"] += 1
            continue
        counts["mapped_reads"] += 1
        counts["mitochondrial_reads" if fields[2] == mt else "nuclear_reads"] += 1
        if int(fields[4]) >= threshold:
            counts["mapped_mapq_ge_threshold"] += 1
            if fields[2] != mt:
                counts["nuclear_mapq_ge_threshold"] += 1
        ops = re.findall(r"(\d+)([MIDNSHP=X])", fields[5])
        if not ops or "".join(n + op for n, op in ops) != fields[5]:
            raise ValueError("Invalid CIGAR")
        if any(op == "H" for _, op in ops):
            raise ValueError("Unexpected hard clipping in primary Bowtie2 output")
        left = int(ops[0][0]) if ops[0][1] == "S" else 0
        right = int(ops[-1][0]) if ops[-1][1] == "S" else 0
        five, three = (right, left) if flag & 16 else (left, right)
        counts["soft_clipped_reads"] += int(five + three > 0)
        counts["five_prime_clipped_reads"] += int(five > 0)
        counts["three_prime_clipped_reads"] += int(three > 0)
        counts["five_prime_clipped_bases"] += five
        counts["three_prime_clipped_bases"] += three
    return counts


def qc(run, bam, output):
    config = json.loads(Path("config/chipseq_alignment.json").read_text())
    plan = json.loads(Path("config/chipseq_alignment_inputs.json").read_text())
    item = next(r for r in plan["runs"] if r["run_accession"] == run)
    process = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE, text=True)
    try:
        counts = scan_sam(process.stdout, config["mitochondrial_accession"], config["diagnostic_mapq"])
    except BaseException:
        process.kill()
        process.wait()
        raise
    finally:
        process.stdout.close()
    if process.wait() != 0:
        raise RuntimeError("Full BAM read failed")
    if counts["primary_reads"] != item["reads"]:
        raise ValueError(f"BAM/FASTQ read count mismatch: {run}")
    fields = ("alignment_records", "nonprimary_records", "primary_reads", "mapped_reads",
              "unmapped_reads", "mitochondrial_reads", "nuclear_reads", "mapped_mapq_ge_threshold",
              "nuclear_mapq_ge_threshold", "soft_clipped_reads", "five_prime_clipped_reads",
              "three_prime_clipped_reads", "five_prime_clipped_bases", "three_prime_clipped_bases")
    data = dict(run_accession=run, role=item["role"], diagnostic_mapq=config["diagnostic_mapq"],
                **{key: counts[key] for key in fields})
    for key in ("mapped_reads", "mitochondrial_reads", "nuclear_mapq_ge_threshold"):
        data[key + "_pct_of_input"] = 100 * counts[key] / counts["primary_reads"]
    data["soft_clipped_pct_of_mapped"] = (100 * counts["soft_clipped_reads"] / counts["mapped_reads"]
                                           if counts["mapped_reads"] else None)
    data["bam_filtering"] = "not_performed"
    data["duplicate_marking"] = "not_performed"
    data["signal_quality"] = "not_evaluated"
    dump(output, data)
    print(f"[OK] Full BAM read and input read count verified: {run}")


def publish(save):
    import yaml
    save = Path(save)
    ref = yaml.safe_load(Path("config/config.yaml").read_text())["reference"]
    rdir = Path(ref["dir"]) / "chipseq"
    out = Path("results/chipseq/alignment")
    prov = out / "reference"
    prov.mkdir(exist_ok=True)
    for path in (rdir / "reference_provenance.json", rdir / "genome_plus_mt.fa.fai",
                 rdir / "mitochondrial.fa", Path(ref["metadata"]), Path(ref["sequence_report"])):
        shutil.copyfile(path, prov / path.name)
    shutil.copyfile(ref["sha256"], prov / "original_nuclear_resources.sha256")
    with (rdir / "genome_plus_mt.fa").open("rb") as source, gzip.open(prov / "genome_plus_mt.fa.gz", "wb", compresslevel=1) as target:
        shutil.copyfileobj(source, target, length=1024 * 1024)
    shutil.copyfile("inputs/verified.json", out / "verified_fastq_inputs.json")
    shutil.copyfile("config/chipseq_alignment_inputs.json", out / "input_provenance.json")
    shutil.copyfile("config/chipseq_alignment.json", out / "alignment_parameters.json")
    manifest = []
    for path in sorted(out.rglob("*")):
        if path.is_file():
            relative = path.relative_to(out)
            target = save / "outputs" / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            digest = sha(path)
            shutil.copyfile(path, target)
            check(target, digest)
            manifest.append(f"{digest}  outputs/{relative.as_posix()}\n")
    (save / "output.sha256").write_text("".join(manifest))
    print(f"[OK] Outputs copied and SHA-256 verified: {save / 'outputs'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("prepare"); p.add_argument("--root", required=True); p.add_argument("--output", required=True)
    sub.add_parser("stage"); sub.add_parser("reference")
    p = sub.add_parser("qc"); p.add_argument("--run", required=True); p.add_argument("--bam", required=True); p.add_argument("--output", required=True)
    p = sub.add_parser("summarize"); p.add_argument("--output", required=True); p.add_argument("reports", nargs="+")
    p = sub.add_parser("publish"); p.add_argument("--save", required=True)
    args = parser.parse_args()
    if args.action == "prepare": dump(args.output, prepare(args.root))
    elif args.action == "stage": stage()
    elif args.action == "reference": mapping_reference()
    elif args.action == "qc": qc(args.run, args.bam, args.output)
    elif args.action == "publish": publish(args.save)
    elif args.action == "summarize":
        rows = [json.loads(Path(path).read_text()) for path in args.reports]
        with Path(args.output).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)


if __name__ == "__main__":
    main()
