"""Small deterministic end-to-end RNA normalization test; writes only to /tmp."""

import argparse
import csv
import gzip
import json
import math
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT_DIR = ROOT / "benchmarks/rnaseq/2026-09-22/rna_experiment_matrix_v1"
MATRIX = EXPERIMENT_DIR / "rnaseq_gene_counts_by_experiment_1473.tsv.gz"
METADATA = EXPERIMENT_DIR / "experiment_metadata.tsv"
SCRIPT = ROOT / "workflow/scripts/normalize_rnaseq_experiments.R"


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_matrix(path):
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = list(reader)
    return header, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="rnaseq-normalization-smoke-") as temp:
        temp = Path(temp)
        selected_nonzero = []
        selected_zero = []
        with gzip.open(MATRIX, "rt", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            experiments = header[1:5]
            for row in reader:
                if all(value == "0" for value in row[1:]):
                    target, limit = selected_zero, 2
                elif any(value != "0" for value in row[1:5]):
                    target, limit = selected_nonzero, 8
                else:
                    continue
                if len(target) < limit:
                    target.append([row[0], *row[1:5]])
        assert len(selected_nonzero) == 8 and len(selected_zero) == 2
        selected = selected_nonzero + selected_zero
        genes = [row[0] for row in selected]
        counts = temp / "counts.tsv.gz"
        with gzip.open(counts, "wt", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["Geneid", *experiments])
            writer.writerows(selected)

        metadata_rows = read_tsv(METADATA)
        selected_metadata = {row["experiment_accession"]: row for row in metadata_rows
                             if row["experiment_accession"] in experiments}
        assert set(selected_metadata) == set(experiments)
        metadata = temp / "metadata.tsv"
        with open(metadata, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["experiment_accession", "study_accession"],
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows({"experiment_accession": experiment,
                             "study_accession": selected_metadata[experiment]["study_accession"]}
                            for experiment in experiments)

        annotation_rows = {}
        with gzip.open(args.annotation, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if row["Geneid"] in genes:
                    annotation_rows[row["Geneid"]] = row
        assert set(annotation_rows) == set(genes)
        annotation = temp / "annotation.tsv.gz"
        with gzip.open(annotation, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["Geneid", "Chr", "Start", "End", "Strand", "Length"],
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(annotation_rows[gene] for gene in genes)

        libraries = [sum(int(row[i + 1]) for row in selected) for i in range(4)]
        assert all(size > 0 for size in libraries)
        qc = temp / "raw_qc.tsv"
        with open(qc, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["experiment_accession", "total_assigned_gene_counts"])
            writer.writerows(zip(experiments, libraries))

        command = [args.rscript, str(SCRIPT), "--counts", str(counts), "--metadata", str(metadata),
                   "--qc", str(qc), "--annotation", str(annotation), "--outdir", str(temp / "normalized"),
                   "--smoke-test"]
        subprocess.run(command, check=True, cwd=ROOT)
        out = temp / "normalized"
        factor_rows = read_tsv(out / "tmm_factors_4.tsv")
        assert [row["experiment_accession"] for row in factor_rows] == experiments
        assert all(math.isfinite(float(row["tmm_norm_factor"])) and
                   float(row["tmm_norm_factor"]) > 0 for row in factor_rows)
        assert all(abs(float(row["effective_library_size"]) -
                       float(row["raw_library_size"]) * float(row["tmm_norm_factor"])) < 1e-6
                   for row in factor_rows)

        tpm_header, tpm_rows = read_matrix(out / "rnaseq_tpm_by_experiment_4.tsv.gz")
        log_header, log_rows = read_matrix(out / "rnaseq_tmm_logcpm_by_experiment_4.tsv.gz")
        for result_header, result_rows in ((tpm_header, tpm_rows), (log_header, log_rows)):
            assert result_header == ["Geneid", *experiments]
            assert [row[0] for row in result_rows] == genes
            assert len(result_rows) == 10 and all(len(row) == 5 for row in result_rows)
            assert all(math.isfinite(float(value)) for row in result_rows for value in row[1:])
        tpm_sums = [sum(float(row[i + 1]) for row in tpm_rows) for i in range(4)]
        assert all(abs(total - 1e6) <= 1e-4 for total in tpm_sums)
        assert all(float(value) == 0 for row in tpm_rows[-2:] for value in row[1:])
        gene_qc = read_tsv(out / "gene_normalization_qc.tsv")
        assert [row["Geneid"] for row in gene_qc] == genes
        assert [row["all_zero_across_experiments"].upper() for row in gene_qc[-2:]] == ["TRUE", "TRUE"]
        summary = json.loads((out / "normalization_qc_summary.json").read_text())
        assert summary["globally_all_zero_genes"] == 2
        assert summary["input_dimensions"]["genes"] == 10
        assert summary["input_dimensions"]["experiments"] == 4
        provenance = json.loads((out / "provenance.json").read_text())
        assert provenance["tmm_method"] == "TMM" and provenance["logcpm_prior_count"] == 2
        assert set(provenance["inputs"]) == {"matrix", "metadata", "raw_count_qc", "annotation"}
        assert set(provenance["package_versions"]) == {"edgeR", "data.table", "jsonlite"}
        assert len(provenance["normalization_script_sha256"]) == 64
        subprocess.run(["sha256sum", "--check", "SHA256SUMS.txt"], cwd=out, check=True,
                       stdout=subprocess.DEVNULL)

        # A mismatched raw-count QC total must fail before any result is published.
        bad_qc = temp / "bad_raw_qc.tsv"
        bad_qc.write_text("experiment_accession\ttotal_assigned_gene_counts\n" +
                          "\n".join(f"{experiment}\t{size + (i == 0)}"
                                    for i, (experiment, size) in enumerate(zip(experiments, libraries))) + "\n")
        bad_command = command.copy()
        bad_command[bad_command.index("--qc") + 1] = str(bad_qc)
        bad_command[bad_command.index("--outdir") + 1] = str(temp / "must_not_exist")
        rejected = subprocess.run(bad_command, cwd=ROOT, capture_output=True, text=True)
        assert rejected.returncode != 0 and not (temp / "must_not_exist").exists()

        print("SMOKE PASS: 10 genes x 4 experiments; 2 all-zero genes retained; TPM sums within 1e-4")
        print("SMOKE PASS: TMM factors/effective sizes finite and positive; order and SHA256SUMS verified")
        print("SMOKE PASS: inconsistent raw-count QC rejected with non-zero exit")


if __name__ == "__main__":
    main()
