#!/usr/bin/env python3

import argparse
import csv
import gzip
import hashlib
import json
import os
import statistics
import subprocess
from collections import Counter


def sha256_file(path):
    h = hashlib.sha256()

    with open(path, "rb") as fh:
        for chunk in iter(
            lambda: fh.read(1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(
            path,
            "rt",
            encoding="utf-8",
            newline="",
        )

    return open(
        path,
        "r",
        encoding="utf-8",
        newline="",
    )


def load_by_key(path, key):
    rows = {}

    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        if key not in (reader.fieldnames or []):
            raise RuntimeError(
                f"{path}: missing key column {key}"
            )

        for row in reader:
            value = row[key]

            if value in rows:
                raise RuntimeError(
                    f"{path}: duplicate {key}={value}"
                )

            rows[value] = row

    return rows


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--locus-audit", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--raw-support", required=True)
    parser.add_argument("--sensitivity", required=True)

    parser.add_argument("--output", required=True)
    parser.add_argument("--qc", required=True)
    parser.add_argument("--provenance", required=True)

    args = parser.parse_args()

    # --------------------------------------------------------
    # 1. Load all validated inputs
    # --------------------------------------------------------

    benchmark = load_by_key(
        args.benchmark,
        "canonical_locus_id",
    )

    locus = load_by_key(
        args.locus_audit,
        "canonical_locus_id",
    )

    consensus = load_by_key(
        args.consensus,
        "Geneid",
    )

    raw = load_by_key(
        args.raw_support,
        "canonical_locus_id",
    )

    sensitivity = load_by_key(
        args.sensitivity,
        "canonical_locus_id",
    )

    # --------------------------------------------------------
    # 2. Validate dimensions and locus identity
    # --------------------------------------------------------

    if len(benchmark) != 37:
        raise RuntimeError(
            f"Expected 37 benchmark loci; found {len(benchmark)}"
        )

    if len(locus) != 37:
        raise RuntimeError(
            f"Expected 37 locus-audit rows; found {len(locus)}"
        )

    if len(raw) != 37:
        raise RuntimeError(
            f"Expected 37 raw-support rows; found {len(raw)}"
        )

    if len(sensitivity) != 37:
        raise RuntimeError(
            f"Expected 37 sensitivity rows; found {len(sensitivity)}"
        )

    if len(consensus) != 27775:
        raise RuntimeError(
            f"Expected 27,775 consensus genes; "
            f"found {len(consensus)}"
        )

    ids = set(benchmark)

    for label, data in [
        ("locus", locus),
        ("raw", raw),
        ("sensitivity", sensitivity),
    ]:
        if set(data) != ids:
            raise RuntimeError(
                f"Benchmark locus set differs from {label}"
            )

    # --------------------------------------------------------
    # 3. Final output schema
    # --------------------------------------------------------

    columns = [
        # Benchmark identity
        "canonical_locus_id",
        "locus_name",
        "benchmark_role",
        "best_evidence_tier",
        "n_observations",
        "target_assembly",
        "target_seqname",
        "target_start_1based",
        "target_end_1based",
        "locus_width_bp",
        "locus_semantic_class",
        "mapping_confidence",
        "supporting_studies",

        # RNA spatial context
        "rna_nearest_gene_id",
        "rna_nearest_gene_tss_1based",
        "rna_nearest_gene_strand",
        "rna_locus_relative_to_tss",
        "rna_min_tss_distance_bp",
        "rna_tss_inside_locus",
        "rna_tss_inside_count",

        # 92-study consensus
        "rna_n_studies",
        "rna_tmm_logcpm_median",
        "rna_tmm_logcpm_q1",
        "rna_tmm_logcpm_q3",
        "rna_tmm_logcpm_iqr",
        "rna_tpm_median",
        "rna_tpm_q1",
        "rna_tpm_q3",
        "rna_tpm_iqr",
        "rna_study_median_tpm_positive_n",
        "rna_study_median_tpm_positive_fraction",

        # Raw experiment support
        "rna_n_experiments",
        "rna_raw_count_total",
        "rna_experiments_raw_count_positive_n",
        "rna_experiments_raw_count_positive_fraction",
        "rna_raw_all_zero",

        # PRJEB35183 sensitivity
        "rna_tmm_logcpm_median_without_PRJEB35183",
        "rna_tmm_logcpm_delta_without_PRJEB35183",
        "rna_tmm_logcpm_abs_delta_without_PRJEB35183",
        "rna_tpm_median_without_PRJEB35183",
        "rna_tpm_delta_without_PRJEB35183",
        "rna_tpm_abs_delta_without_PRJEB35183",
        "rna_tpm_relative_change_without_PRJEB35183",
    ]

    final_rows = []

    # --------------------------------------------------------
    # 4. Preserve canonical benchmark order
    # --------------------------------------------------------

    with open_text(args.benchmark) as fh:

        reader = csv.DictReader(fh, delimiter="\t")

        for b in reader:

            locus_id = b["canonical_locus_id"]

            l = locus[locus_id]
            r = raw[locus_id]
            s = sensitivity[locus_id]

            # ------------------------------------------------
            # Coordinate consistency
            # ------------------------------------------------

            if l["seqname"] != b["target_seqname"]:
                raise RuntimeError(
                    f"{locus_id}: seqname mismatch"
                )

            if int(l["start_1based"]) != int(
                b["target_start"]
            ):
                raise RuntimeError(
                    f"{locus_id}: start mismatch"
                )

            if int(l["end_1based"]) != int(
                b["target_end"]
            ):
                raise RuntimeError(
                    f"{locus_id}: end mismatch"
                )

            gene = l["nearest_gene_ids"]

            if int(l["nearest_gene_count"]) != 1:
                raise RuntimeError(
                    f"{locus_id}: nearest gene not unique"
                )

            if ";" in gene:
                raise RuntimeError(
                    f"{locus_id}: multiple Geneids"
                )

            # ------------------------------------------------
            # Cross-input gene consistency
            # ------------------------------------------------

            if r["nearest_gene_id"] != gene:
                raise RuntimeError(
                    f"{locus_id}: raw-support gene mismatch"
                )

            if s["nearest_gene_id"] != gene:
                raise RuntimeError(
                    f"{locus_id}: sensitivity gene mismatch"
                )

            if gene not in consensus:
                raise RuntimeError(
                    f"{locus_id}: consensus missing {gene}"
                )

            c = consensus[gene]

            if int(c["n_studies"]) != 92:
                raise RuntimeError(
                    f"{gene}: consensus n_studies != 92"
                )

            if int(r["n_experiments"]) != 1473:
                raise RuntimeError(
                    f"{gene}: n_experiments != 1473"
                )

            # ------------------------------------------------
            # TSS uniqueness already validated upstream
            # ------------------------------------------------

            tss = l["nearest_tss_positions"]

            if ";" in tss:
                raise RuntimeError(
                    f"{locus_id}: multiple nearest TSS values"
                )

            tss_inside_count = int(
                l["tss_inside_count"]
            )

            # ------------------------------------------------
            # Assemble final row
            # ------------------------------------------------

            row = {
                "canonical_locus_id":
                    locus_id,
                "locus_name":
                    b["locus_name"],
                "benchmark_role":
                    b["benchmark_role"],
                "best_evidence_tier":
                    b["best_evidence_tier"],
                "n_observations":
                    b["n_observations"],
                "target_assembly":
                    b["target_assembly"],
                "target_seqname":
                    b["target_seqname"],
                "target_start_1based":
                    b["target_start"],
                "target_end_1based":
                    b["target_end"],
                "locus_width_bp":
                    l["width_bp"],
                "locus_semantic_class":
                    l["semantic_class"],
                "mapping_confidence":
                    b["mapping_confidence"],
                "supporting_studies":
                    b["supporting_studies"],

                "rna_nearest_gene_id":
                    gene,
                "rna_nearest_gene_tss_1based":
                    tss,
                "rna_nearest_gene_strand":
                    l["nearest_gene_strands"],
                "rna_locus_relative_to_tss":
                    l["nearest_locus_relative_to_tss"],
                "rna_min_tss_distance_bp":
                    l["min_tss_distance_bp"],
                "rna_tss_inside_locus":
                    "TRUE"
                    if tss_inside_count > 0
                    else "FALSE",
                "rna_tss_inside_count":
                    l["tss_inside_count"],

                "rna_n_studies":
                    c["n_studies"],
                "rna_tmm_logcpm_median":
                    c["tmm_logcpm_median"],
                "rna_tmm_logcpm_q1":
                    c["tmm_logcpm_q1"],
                "rna_tmm_logcpm_q3":
                    c["tmm_logcpm_q3"],
                "rna_tmm_logcpm_iqr":
                    c["tmm_logcpm_iqr"],
                "rna_tpm_median":
                    c["tpm_median"],
                "rna_tpm_q1":
                    c["tpm_q1"],
                "rna_tpm_q3":
                    c["tpm_q3"],
                "rna_tpm_iqr":
                    c["tpm_iqr"],
                "rna_study_median_tpm_positive_n":
                    c["study_median_tpm_positive_n"],
                "rna_study_median_tpm_positive_fraction":
                    c[
                        "study_median_tpm_positive_fraction"
                    ],

                "rna_n_experiments":
                    r["n_experiments"],
                "rna_raw_count_total":
                    r["raw_count_total"],
                "rna_experiments_raw_count_positive_n":
                    r[
                        "experiments_raw_count_positive_n"
                    ],
                "rna_experiments_raw_count_positive_fraction":
                    r[
                        "experiments_raw_count_positive_fraction"
                    ],
                "rna_raw_all_zero":
                    r["raw_all_zero"],

                "rna_tmm_logcpm_median_without_PRJEB35183":
                    s[
                        "tmm_logcpm_median_91_without_PRJEB35183"
                    ],
                "rna_tmm_logcpm_delta_without_PRJEB35183":
                    s["tmm_logcpm_delta"],
                "rna_tmm_logcpm_abs_delta_without_PRJEB35183":
                    s["tmm_logcpm_abs_delta"],

                "rna_tpm_median_without_PRJEB35183":
                    s[
                        "tpm_median_91_without_PRJEB35183"
                    ],
                "rna_tpm_delta_without_PRJEB35183":
                    s["tpm_delta"],
                "rna_tpm_abs_delta_without_PRJEB35183":
                    s["tpm_abs_delta"],
                "rna_tpm_relative_change_without_PRJEB35183":
                    s["tpm_relative_change"],
            }

            final_rows.append(row)

    # --------------------------------------------------------
    # 5. Final structural validation
    # --------------------------------------------------------

    if len(final_rows) != 37:
        raise RuntimeError(
            f"Expected 37 final rows; "
            f"found {len(final_rows)}"
        )

    if len({
        r["canonical_locus_id"]
        for r in final_rows
    }) != 37:
        raise RuntimeError(
            "Duplicate canonical locus IDs"
        )

    if len({
        r["rna_nearest_gene_id"]
        for r in final_rows
    }) != 37:
        raise RuntimeError(
            "Expected 37 unique nearest genes"
        )

    # --------------------------------------------------------
    # 6. Write final RNA context table
    # --------------------------------------------------------

    os.makedirs(
        os.path.dirname(args.output),
        exist_ok=True,
    )

    with open(
        args.output,
        "w",
        encoding="utf-8",
        newline="",
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(final_rows)

    # --------------------------------------------------------
    # 7. QC statistics
    # --------------------------------------------------------

    roles = Counter(
        r["benchmark_role"]
        for r in final_rows
    )

    semantic = Counter(
        r["locus_semantic_class"]
        for r in final_rows
    )

    all_zero = sum(
        r["rna_raw_all_zero"] == "TRUE"
        for r in final_rows
    )

    tpm_positive = sum(
        float(r["rna_tpm_median"]) > 0
        for r in final_rows
    )

    tpm_zero = 37 - tpm_positive

    tss_inside = sum(
        r["rna_tss_inside_locus"] == "TRUE"
        for r in final_rows
    )

    distances = [
        int(r["rna_min_tss_distance_bp"])
        for r in final_rows
    ]

    positive_experiments = [
        int(
            r[
                "rna_experiments_raw_count_positive_n"
            ]
        )
        for r in final_rows
    ]

    tmm_sensitivity = [
        float(
            r[
                "rna_tmm_logcpm_abs_delta_without_PRJEB35183"
            ]
        )
        for r in final_rows
    ]

    with open(
        args.qc,
        "w",
        encoding="utf-8",
    ) as fh:

        fh.write("metric\tvalue\n")

        fh.write("benchmark_loci\t37\n")
        fh.write("unique_nearest_genes\t37\n")

        for role in sorted(roles):
            fh.write(
                f"benchmark_role_{role}\t"
                f"{roles[role]}\n"
            )

        for name in sorted(semantic):
            fh.write(
                f"semantic_{name}\t"
                f"{semantic[name]}\n"
            )

        fh.write(
            f"loci_with_tss_inside\t"
            f"{tss_inside}\n"
        )

        fh.write(
            f"nearest_genes_raw_all_zero\t"
            f"{all_zero}\n"
        )

        fh.write(
            f"nearest_genes_tpm_median_positive\t"
            f"{tpm_positive}\n"
        )

        fh.write(
            f"nearest_genes_tpm_median_zero\t"
            f"{tpm_zero}\n"
        )

        fh.write(
            f"minimum_tss_distance_bp\t"
            f"{min(distances)}\n"
        )

        fh.write(
            f"maximum_tss_distance_bp\t"
            f"{max(distances)}\n"
        )

        fh.write(
            f"minimum_positive_experiments\t"
            f"{min(positive_experiments)}\n"
        )

        fh.write(
            f"maximum_positive_experiments\t"
            f"{max(positive_experiments)}\n"
        )

        fh.write(
            "PRJEB35183_tmm_abs_delta_median\t"
            f"{statistics.median(tmm_sensitivity):.17g}\n"
        )

        fh.write(
            "PRJEB35183_tmm_abs_delta_max\t"
            f"{max(tmm_sensitivity):.17g}\n"
        )

    # --------------------------------------------------------
    # 8. Provenance
    # --------------------------------------------------------

    try:
        pipeline_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
        ).strip()
    except Exception:
        pipeline_commit = "unavailable"

    provenance = {
        "product":
            "benchmark_37_rna_context_v1",
        "description":
            (
                "Gene-level RNA-seq context for the 37 "
                "canonical CHO benchmark loci. RNA signal "
                "is associated by nearest gene-level TSS "
                "and must not be interpreted as localized "
                "transcriptional signal at the integration site."
            ),
        "pipeline_commit":
            pipeline_commit,
        "method": {
            "rna_resolution":
                "gene_id",
            "tss_definition":
                (
                    "5-prime boundary of GTF gene feature; "
                    "validated as most-upstream annotated "
                    "transcript TSS for all quantified genes"
                ),
            "locus_distance":
                (
                    "minimum distance from gene-level TSS "
                    "to complete canonical 1-based inclusive "
                    "locus interval"
                ),
            "study_consensus":
                (
                    "median across 92 study-level medians "
                    "with equal study weighting"
                ),
            "sensitivity":
                "repeat consensus excluding PRJEB35183",
        },
        "sources": {
            "benchmark": {
                "path": args.benchmark,
                "sha256": sha256_file(
                    args.benchmark
                ),
            },
            "locus_audit": {
                "path": args.locus_audit,
                "sha256": sha256_file(
                    args.locus_audit
                ),
            },
            "consensus": {
                "path": args.consensus,
                "sha256": sha256_file(
                    args.consensus
                ),
            },
            "raw_support": {
                "path": args.raw_support,
                "sha256": sha256_file(
                    args.raw_support
                ),
            },
            "sensitivity": {
                "path": args.sensitivity,
                "sha256": sha256_file(
                    args.sensitivity
                ),
            },
            "script": {
                "path": os.path.abspath(__file__),
                "sha256": sha256_file(__file__),
            },
        },
    }

    with open(
        args.provenance,
        "w",
        encoding="utf-8",
    ) as fh:

        json.dump(
            provenance,
            fh,
            indent=2,
            sort_keys=True,
        )

        fh.write("\n")

    # --------------------------------------------------------
    # 9. Console report
    # --------------------------------------------------------

    print("BENCHMARK_LOCI=37")
    print("UNIQUE_NEAREST_GENES=37")

    print()
    print("BENCHMARK_ROLES")

    for role in sorted(roles):
        print(f"{role}={roles[role]}")

    print()
    print(
        f"LOCI_WITH_TSS_INSIDE={tss_inside}"
    )

    print(
        f"NEAREST_GENES_RAW_ALL_ZERO={all_zero}"
    )

    print(
        "NEAREST_GENES_TPM_MEDIAN_POSITIVE="
        f"{tpm_positive}"
    )

    print(
        "NEAREST_GENES_TPM_MEDIAN_ZERO="
        f"{tpm_zero}"
    )

    print(
        "RAW_SUPPORT_RANGE="
        f"{min(positive_experiments)}-"
        f"{max(positive_experiments)}/1473"
    )

    print()
    print(
        "PASS: final 37-locus RNA context "
        "table built and validated"
    )


if __name__ == "__main__":
    main()
