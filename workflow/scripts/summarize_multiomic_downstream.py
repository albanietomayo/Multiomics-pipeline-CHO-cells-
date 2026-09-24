#!/usr/bin/env python3
import argparse
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list

RNA_FEATURES = [
    "rna_min_tss_distance_bp",
    "rna_tmm_logcpm_median",
    "rna_tpm_median",
    "rna_study_median_tpm_positive_fraction",
]
ATAC_FEATURES = [
    "atac_idr10_exact_n_unique_peaks",
    "atac_idr10_exact_overlap_fraction_of_locus",
    "atac_idr10_accessible_bp_pm1000bp",
    "atac_spmr_consensus_pm1000bp_mean",
    "atac_idr10_accessible_bp_pm10000bp",
    "atac_spmr_consensus_pm10000bp_mean",
    "atac_idr10_accessible_bp_pm50000bp",
    "atac_spmr_consensus_pm50000bp_mean",
]
MARKS = ["H3K27ac", "H3K4me3", "H3K36me3", "H3K4me1", "H3K27me3", "H3K9me3"]
CHIP_FEATURES = []
for mark in MARKS:
    CHIP_FEATURES.extend([
        f"{mark}__chip_any_peak_pm50000bp_fraction_positive_conditions_equal_study_weight",
        f"{mark}__chip_enriched_fraction_pm10000bp_consensus_median_across_studies",
    ])
HEATMAP_FEATURES = RNA_FEATURES + ATAC_FEATURES + CHIP_FEATURES
HEATMAP_LABELS = [
    "RNA: TSS distance", "RNA: logCPM", "RNA: TPM", "RNA: detected",
    "ATAC: exact peaks", "ATAC: exact overlap",
    "ATAC: bp +/-1 kb", "ATAC: SPMR +/-1 kb",
    "ATAC: bp +/-10 kb", "ATAC: SPMR +/-10 kb",
    "ATAC: bp +/-50 kb", "ATAC: SPMR +/-50 kb",
]
for mark in MARKS:
    HEATMAP_LABELS.extend([f"{mark}: peak +/-50 kb", f"{mark}: enrich +/-10 kb"])


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def median(df, col):
    return float(df[col].median())


def count_gt0(df, col):
    return int((df[col] > 0).sum())


def build_dhiman_summary(core: pd.DataFrame, out: Path):
    dh = core[core["canonical_locus_id"].str.startswith("CHO_BM_DHIMAN_")].copy()
    pos = dh[dh["benchmark_role"] == "positive"]
    neg = dh[dh["benchmark_role"] == "negative"]
    sup = dh[dh["benchmark_role"] == "support_only"]
    if (len(pos), len(neg), len(sup)) != (6, 4, 3):
        raise RuntimeError(f"Unexpected Dhiman composition: positive={len(pos)}, negative={len(neg)}, support_only={len(sup)}")

    rows = []
    def add(metric, pval, nval, unit, statistic):
        rows.append({
            "metric": metric,
            "statistic": statistic,
            "positive_n": len(pos),
            "positive_value": pval,
            "negative_n": len(neg),
            "negative_value": nval,
            "unit_or_scale": unit,
        })

    add("RNA TPM", median(pos, "rna_tpm_median"), median(neg, "rna_tpm_median"), "TPM", "median")
    add("RNA detectable-study fraction", median(pos, "rna_study_median_tpm_positive_fraction"), median(neg, "rna_study_median_tpm_positive_fraction"), "fraction", "median")
    add("RNA loci with TPM > 0", count_gt0(pos, "rna_tpm_median"), count_gt0(neg, "rna_tpm_median"), "loci", "count")
    add("ATAC accessible bp +/-50 kb", median(pos, "atac_idr10_accessible_bp_pm50000bp"), median(neg, "atac_idr10_accessible_bp_pm50000bp"), "bp", "median")
    add("ATAC SPMR +/-50 kb", median(pos, "atac_spmr_consensus_pm50000bp_mean"), median(neg, "atac_spmr_consensus_pm50000bp_mean"), "SPMR", "median")
    add("ATAC loci with accessible bp +/-50 kb > 0", count_gt0(pos, "atac_idr10_accessible_bp_pm50000bp"), count_gt0(neg, "atac_idr10_accessible_bp_pm50000bp"), "loci", "count")
    for mark in ["H3K27ac", "H3K4me3", "H3K36me3"]:
        col = f"{mark}__chip_nearest_peak_distance_to_locus_bp_consensus_median_across_studies"
        add(f"{mark} nearest peak distance", median(pos, col), median(neg, col), "bp", "median")
    for mark in ["H3K9me3", "H3K27me3"]:
        col_frac = f"{mark}__chip_enriched_fraction_pm10000bp_consensus_median_across_studies"
        col_any = f"{mark}__chip_any_peak_pm10000bp_fraction_positive_conditions_equal_study_weight"
        add(f"{mark} enriched fraction +/-10 kb", median(pos, col_frac), median(neg, col_frac), "fraction", "median")
        add(f"{mark} loci with any peak +/-10 kb", count_gt0(pos, col_any), count_gt0(neg, col_any), "loci", "count")

    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


def build_heatmap(core: pd.DataFrame, outdir: Path):
    missing = [c for c in HEATMAP_FEATURES if c not in core.columns]
    if missing:
        raise RuntimeError(f"Missing Core features: {missing}")
    X = core[HEATMAP_FEATURES].astype(float)
    std = X.std(axis=0, ddof=0)
    if (std == 0).any():
        raise RuntimeError(f"Zero-variance heatmap features: {list(std[std == 0].index)}")
    Z = (X - X.mean(axis=0)) / std

    order = leaves_list(linkage(Z.values, method="ward", metric="euclidean"))
    Z_ord = Z.iloc[order].copy()
    meta = core.iloc[order][["canonical_locus_id", "benchmark_role"]].copy()

    matrix = Z_ord.copy()
    matrix.insert(0, "benchmark_role", meta["benchmark_role"].values)
    matrix.insert(0, "canonical_locus_id", meta["canonical_locus_id"].values)
    matrix.to_csv(outdir / "Figure_5_heatmap_matrix_estandarizada.tsv", sep="\t", index=False)

    codes = {"positive": "P", "negative": "N", "support_only": "S"}
    row_labels = [
        f"{codes[r]}  {l.replace('CHO_BM_', '')}"
        for l, r in zip(meta["canonical_locus_id"], meta["benchmark_role"])
    ]

    fig, ax = plt.subplots(figsize=(17, 12))
    im = ax.imshow(Z_ord.values, aspect="auto")
    ax.set_xticks(np.arange(len(HEATMAP_LABELS)))
    ax.set_xticklabels(HEATMAP_LABELS, rotation=65, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_title("Contexto regulador multi-omico de los 37 loci del benchmark", fontsize=13, pad=10)
    ax.set_xlabel("Caracteristicas de la vista Core estandarizadas por columna")
    ax.set_ylabel("Loci (P = positive; N = negative; S = support_only)")
    for boundary in [3.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21.5]:
        ax.axvline(boundary, linewidth=0.8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("z-score por caracteristica")
    fig.text(
        0.5, 0.008,
        "Clustering jerarquico no supervisado de los loci (Ward, distancia euclidea). "
        "benchmark_role se muestra solo como anotacion y no interviene en el clustering.",
        ha="center", fontsize=8,
    )
    fig.tight_layout(rect=[0, 0.035, 1, 0.98])
    fig.savefig(outdir / "Figure_5_heatmap_multiomic_37_loci.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "Figure_5_heatmap_multiomic_37_loci.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--core", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    core = pd.read_csv(args.core, sep="\t")
    if core.shape != (37, 69):
        raise RuntimeError(f"Unexpected Core shape: {core.shape}; expected (37, 69)")
    expected_roles = {"positive": 27, "negative": 4, "support_only": 6}
    observed = core["benchmark_role"].value_counts().to_dict()
    if observed != expected_roles:
        raise RuntimeError(f"Unexpected benchmark roles: {observed}")

    build_dhiman_summary(core, args.outdir / "dhiman_within_study_summary_v1.tsv")
    build_heatmap(core, args.outdir)

    outputs = [
        args.outdir / "dhiman_within_study_summary_v1.tsv",
        args.outdir / "Figure_5_heatmap_matrix_estandarizada.tsv",
        args.outdir / "Figure_5_heatmap_multiomic_37_loci.png",
        args.outdir / "Figure_5_heatmap_multiomic_37_loci.pdf",
    ]
    with (args.outdir / "SHA256SUMS.txt").open("w") as fh:
        for p in outputs:
            fh.write(f"{sha256(p)}  {p.name}\n")
    print("DOWNSTREAM_MULTIOMIC_EXPLORATION=PASS")
    print(f"core_shape={core.shape[0]}x{core.shape[1]}")
    print(f"heatmap_features={len(HEATMAP_FEATURES)}")
    print("label_used_for_clustering=NO")
    print("dhiman_positive=6")
    print("dhiman_negative=4")
    print("dhiman_support_only=3")

if __name__ == "__main__":
    main()
