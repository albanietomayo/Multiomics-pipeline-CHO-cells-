#!/usr/bin/env python3

import argparse
import csv
import gzip
import hashlib
import json
import math
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def assert_nonempty(path, label):
    path = Path(path)
    if not path.is_file() or path.stat().st_size <= 0:
        raise RuntimeError(f"{label} missing or empty: {path}")


def load_metadata(path):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    required = {
        "experiment_accession",
        "study_accession",
        "sample_accession",
    }

    if not rows:
        raise RuntimeError("Metadata is empty")

    missing = required.difference(rows[0])
    if missing:
        raise RuntimeError(
            f"Metadata missing required columns: {sorted(missing)}"
        )

    experiments = [r["experiment_accession"] for r in rows]

    if any(not x for x in experiments):
        raise RuntimeError("Empty experiment_accession in metadata")

    if len(experiments) != len(set(experiments)):
        raise RuntimeError("Duplicate experiment_accession in metadata")

    if any(not r["study_accession"] for r in rows):
        raise RuntimeError("Empty study_accession in metadata")

    return rows


def build_study_structure(metadata):
    exp_to_study = {}
    study_to_exps = defaultdict(list)
    study_to_samples = defaultdict(set)

    for row in metadata:
        exp = row["experiment_accession"]
        study = row["study_accession"]
        sample = row["sample_accession"]

        exp_to_study[exp] = study
        study_to_exps[study].append(exp)

        if sample:
            study_to_samples[study].add(sample)

    study_order = sorted(study_to_exps)

    return (
        exp_to_study,
        study_to_exps,
        study_to_samples,
        study_order,
    )


def group_indices(experiment_header, exp_to_study, study_order):
    if len(experiment_header) != len(set(experiment_header)):
        raise RuntimeError("Duplicate experiment columns in matrix")

    if set(experiment_header) != set(exp_to_study):
        missing = sorted(set(exp_to_study) - set(experiment_header))
        extra = sorted(set(experiment_header) - set(exp_to_study))

        raise RuntimeError(
            "Matrix/metadata experiment mismatch. "
            f"Missing={missing[:10]} Extra={extra[:10]}"
        )

    grouped = {study: [] for study in study_order}

    for i, exp in enumerate(experiment_header):
        grouped[exp_to_study[exp]].append(i)

    indices = []

    for study in study_order:
        idx = grouped[study]
        if not idx:
            raise RuntimeError(f"Study has no matrix columns: {study}")
        indices.append(idx)

    return indices


def median_values(values, grouped_indices):
    result = []

    for idx in grouped_indices:
        x = [values[i] for i in idx]
        result.append(statistics.median(x))

    return result


def aggregate_matrix(
    input_path,
    output_path,
    exp_to_study,
    study_order,
    expected_experiment_header=None,
    expected_genes=None,
    require_nonnegative=False,
):
    genes = []
    column_sums = [0.0] * len(study_order)

    with gzip.open(
        input_path,
        "rt",
        newline="",
        encoding="utf-8",
    ) as src, gzip.open(
        output_path,
        "wt",
        newline="",
        encoding="utf-8",
    ) as dst:

        reader = csv.reader(src, delimiter="\t")
        writer = csv.writer(
            dst,
            delimiter="\t",
            lineterminator="\n",
        )

        try:
            header = next(reader)
        except StopIteration:
            raise RuntimeError(f"Empty matrix: {input_path}")

        if not header or header[0] != "Geneid":
            raise RuntimeError(
                f"First matrix column must be Geneid: {input_path}"
            )

        experiment_header = header[1:]

        if expected_experiment_header is not None:
            if experiment_header != expected_experiment_header:
                raise RuntimeError(
                    "Experiment column order differs between matrices"
                )

        grouped_indices = group_indices(
            experiment_header,
            exp_to_study,
            study_order,
        )

        writer.writerow(["Geneid", *study_order])

        seen_genes = set()

        for row_index, row in enumerate(reader):
            if len(row) != len(header):
                raise RuntimeError(
                    f"Invalid column count at row {row_index + 2}: "
                    f"{len(row)} != {len(header)}"
                )

            gene = row[0]

            if not gene:
                raise RuntimeError(
                    f"Empty Geneid at row {row_index + 2}"
                )

            if gene in seen_genes:
                raise RuntimeError(f"Duplicate Geneid: {gene}")

            seen_genes.add(gene)

            if expected_genes is not None:
                if row_index >= len(expected_genes):
                    raise RuntimeError(
                        "Matrix contains more genes than expected"
                    )

                if gene != expected_genes[row_index]:
                    raise RuntimeError(
                        "Gene order mismatch at row "
                        f"{row_index + 2}: "
                        f"{gene} != {expected_genes[row_index]}"
                    )

            values = []

            for text in row[1:]:
                value = float(text)

                if not math.isfinite(value):
                    raise RuntimeError(
                        f"Non-finite value for gene {gene}"
                    )

                if require_nonnegative and value < 0:
                    raise RuntimeError(
                        f"Negative TPM value for gene {gene}"
                    )

                values.append(value)

            aggregated = median_values(
                values,
                grouped_indices,
            )

            if any(not math.isfinite(x) for x in aggregated):
                raise RuntimeError(
                    f"Non-finite study median for gene {gene}"
                )

            if require_nonnegative and any(
                x < 0 for x in aggregated
            ):
                raise RuntimeError(
                    f"Negative study-median TPM for gene {gene}"
                )

            writer.writerow(
                [gene] + [format(x, ".17g") for x in aggregated]
            )

            if require_nonnegative:
                for j, x in enumerate(aggregated):
                    column_sums[j] += x

            genes.append(gene)

    if expected_genes is not None:
        if genes != expected_genes:
            raise RuntimeError(
                "Gene set/order differs between matrices"
            )

    return experiment_header, genes, column_sums


def write_study_metadata(
    path,
    study_order,
    study_to_exps,
    study_to_samples,
):
    with open(
        path,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "study_accession",
                "n_experiments",
                "n_samples",
                "experiment_accessions",
            ]
        )

        for study in study_order:
            experiments = sorted(study_to_exps[study])

            writer.writerow(
                [
                    study,
                    len(experiments),
                    len(study_to_samples[study]),
                    ";".join(experiments),
                ]
            )


def validate_gzip(path):
    result = subprocess.run(
        ["gzip", "-t", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"gzip integrity failure: {path}\n{result.stderr}"
        )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--logcpm", required=True)
    parser.add_argument("--tpm", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--outdir", required=True)

    args = parser.parse_args()

    logcpm = Path(args.logcpm).resolve()
    tpm = Path(args.tpm).resolve()
    metadata_path = Path(args.metadata).resolve()
    outdir = Path(args.outdir).resolve()

    assert_nonempty(logcpm, "logCPM matrix")
    assert_nonempty(tpm, "TPM matrix")
    assert_nonempty(metadata_path, "experiment metadata")

    if outdir.exists():
        raise RuntimeError(
            f"Output directory already exists: {outdir}"
        )

    outdir.parent.mkdir(parents=True, exist_ok=True)

    tmp = outdir.with_name(
        outdir.name + f".tmp-{os.getpid()}"
    )

    if tmp.exists():
        raise RuntimeError(
            f"Temporary output directory already exists: {tmp}"
        )

    tmp.mkdir()

    try:
        metadata = load_metadata(metadata_path)

        (
            exp_to_study,
            study_to_exps,
            study_to_samples,
            study_order,
        ) = build_study_structure(metadata)

        n_experiments = len(exp_to_study)
        n_studies = len(study_order)

        log_name = (
            f"rnaseq_tmm_logcpm_by_study_median_"
            f"{n_studies}.tsv.gz"
        )
        tpm_name = (
            f"rnaseq_tpm_by_study_median_"
            f"{n_studies}.tsv.gz"
        )

        log_out = tmp / log_name
        tpm_out = tmp / tpm_name

        experiment_header, genes, _ = aggregate_matrix(
            logcpm,
            log_out,
            exp_to_study,
            study_order,
            expected_experiment_header=None,
            expected_genes=None,
            require_nonnegative=False,
        )

        _, tpm_genes, tpm_column_sums = aggregate_matrix(
            tpm,
            tpm_out,
            exp_to_study,
            study_order,
            expected_experiment_header=experiment_header,
            expected_genes=genes,
            require_nonnegative=True,
        )

        if tpm_genes != genes:
            raise RuntimeError(
                "TPM/logCPM gene order mismatch"
            )

        study_metadata = tmp / "study_metadata.tsv"

        write_study_metadata(
            study_metadata,
            study_order,
            study_to_exps,
            study_to_samples,
        )

        study_sizes = [
            len(study_to_exps[s])
            for s in study_order
        ]

        summary = {
            "input_dimensions": {
                "genes": len(genes),
                "experiments": n_experiments,
                "studies": n_studies,
            },
            "aggregation_level": "study",
            "aggregation_statistic": "median",
            "aggregation_definition": (
                "Gene-wise median across normalized "
                "experiment-level values within each study."
            ),
            "study_ordering": (
                "lexicographic study_accession"
            ),
            "study_size_experiments": {
                "min": min(study_sizes),
                "median": statistics.median(study_sizes),
                "max": max(study_sizes),
            },
            "post_aggregation_renormalization": False,
            "tpm_interpretation": (
                "Study-median TPM summaries are not "
                "renormalized after gene-wise median "
                "aggregation and are therefore not expected "
                "to sum exactly to 1e6."
            ),
            "study_median_tpm_column_sum": {
                "min": min(tpm_column_sums),
                "median": statistics.median(
                    tpm_column_sums
                ),
                "max": max(tpm_column_sums),
            },
        }

        summary_path = (
            tmp / "study_aggregation_qc_summary.json"
        )

        with open(
            summary_path,
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump(
                summary,
                handle,
                indent=2,
                sort_keys=True,
            )
            handle.write("\n")

        repo_root = (
            Path(__file__).resolve().parents[2]
        )

        git_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            text=True,
        ).strip()

        provenance = {
            "git_commit": git_commit,
            "aggregation_script": str(
                Path(__file__).resolve().relative_to(
                    repo_root
                )
            ),
            "aggregation_script_sha256": sha256_file(
                Path(__file__).resolve()
            ),
            "python_version": sys.version,
            "inputs": {
                "logcpm": {
                    "path": str(logcpm),
                    "sha256": sha256_file(logcpm),
                },
                "tpm": {
                    "path": str(tpm),
                    "sha256": sha256_file(tpm),
                },
                "experiment_metadata": {
                    "path": str(metadata_path),
                    "sha256": sha256_file(
                        metadata_path
                    ),
                },
            },
            "aggregation": {
                "level": "study",
                "statistic": "median",
                "post_aggregation_renormalization": False,
                "study_weighting": (
                    "One output column per study; "
                    "equal study weighting is applied "
                    "when summarizing across study columns."
                ),
            },
        }

        provenance_path = tmp / "provenance.json"

        with open(
            provenance_path,
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump(
                provenance,
                handle,
                indent=2,
                sort_keys=True,
            )
            handle.write("\n")

        output_files = [
            log_out,
            tpm_out,
            study_metadata,
            summary_path,
            provenance_path,
        ]

        for path in output_files:
            assert_nonempty(path, path.name)

        validate_gzip(log_out)
        validate_gzip(tpm_out)

        checksum_path = tmp / "SHA256SUMS.txt"

        with open(
            checksum_path,
            "w",
            encoding="utf-8",
        ) as handle:
            for path in output_files:
                handle.write(
                    f"{sha256_file(path)}  {path.name}\n"
                )

        assert_nonempty(
            checksum_path,
            "SHA256SUMS.txt",
        )

        expected_names = {
            log_name,
            tpm_name,
            "study_metadata.tsv",
            "study_aggregation_qc_summary.json",
            "provenance.json",
            "SHA256SUMS.txt",
        }

        observed_names = {
            p.name for p in tmp.iterdir()
            if p.is_file()
        }

        if observed_names != expected_names:
            raise RuntimeError(
                "Unexpected output set. "
                f"Expected={sorted(expected_names)} "
                f"Observed={sorted(observed_names)}"
            )

        tmp.rename(outdir)

        print(
            f"Validated and wrote "
            f"{len(genes)} genes x "
            f"{n_studies} studies to {outdir}"
        )

    except Exception:
        shutil.rmtree(tmp, ignore_errors=True)
        raise


if __name__ == "__main__":
    import os
    main()
