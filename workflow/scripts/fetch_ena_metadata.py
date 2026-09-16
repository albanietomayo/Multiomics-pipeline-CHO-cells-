#!/usr/bin/env python3

"""
Retrieve sequencing metadata from the ENA Portal API
and generate a structured metadata table.

TFM multi-omics pipeline.
Biological target and classification rules are defined externally.
"""

from pathlib import Path
import io
import sys

import pandas as pd
import requests
import yaml

from selection_policy import reclassify_dataframe, save_workflow_audit


# ===================================================================
# Paths - Localiza la raiz del proyecto y el archivo de configuracion
# ===================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

CONFIG_FILE = PROJECT_ROOT / "config" / "config.yaml"


# ============================================================
# Configuration- Lee los parametros de config.yaml
# ============================================================

def load_config(config_file):
    """Load pipeline configuration from YAML."""

    with open(config_file, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


# ============================================================
# ENA query - Construye la consulta a ENA
# ============================================================

def fetch_ena_metadata(taxid, result_type):
    """
    Query ENA Portal API using an NCBI taxonomy identifier.
    """

    url = "https://www.ebi.ac.uk/ena/portal/api/search"

    fields = [
        "study_accession",
        "sample_accession",
        "experiment_accession",
        "run_accession",
        "scientific_name",
        "tax_id",
        "study_title",
        "sample_title",
        "experiment_title",
        "library_name",
        "cell_line",
        "cell_type",
        "library_strategy",
        "library_source",
        "library_selection",
        "library_layout",
        "instrument_platform",
        "instrument_model",
        "fastq_ftp",
        "fastq_md5",
        "fastq_bytes",
    ]

    params = {
        "result": result_type,
        "query": f"tax_eq({taxid})",
        "fields": ",".join(fields),
        "format": "tsv",
        "limit": "0",
    }

    print(f"[INFO] Querying ENA for taxid {taxid}...")

    response = requests.get(
        url,
        params=params,
        timeout=120,
    )

    response.raise_for_status()

    if not response.text.strip():
        raise RuntimeError(
            f"ENA returned no records for taxid {taxid}."
        )

    dataframe = pd.read_csv(
        io.StringIO(response.text),
        sep="\t",
        dtype=str,
    )

    print(f"[INFO] Records retrieved from ENA: {len(dataframe)}")

    return dataframe


# ============================================================
# Preliminary filtering
# ============================================================

def filter_fastq_records(dataframe):
    """
    Retain records for which raw FASTQ files are available.
    """

    filtered = dataframe.copy()

    filtered = filtered[
        filtered["fastq_ftp"].notna()
        & (filtered["fastq_ftp"].str.strip() != "")
    ]

    print(
        "[INFO] Records with available FASTQ files: "
        f"{len(filtered)}"
    )

    return filtered


# ============================================================
# Omics classification
# ============================================================

def classify_omics(row):
    """
    Classify sequencing experiments using structured ENA metadata.

    Single-cell experiments are excluded because they require
    dedicated analytical workflows not included in this pipeline.
    """

    strategy = str(row.get("library_strategy", "")).strip().lower()
    source = str(row.get("library_source", "")).strip().lower()

    if "single cell" in source or "single-cell" in source:
        return "excluded_single_cell"

    if strategy == "rna-seq":
        return "RNA-seq"

    if strategy == "atac-seq":
        return "ATAC-seq"

    if strategy == "chip-seq":
        return "ChIP-seq"

    return "unclassified"


def add_omics_classification(dataframe):
    """Add normalized omics classification."""

    dataframe = dataframe.copy()

    dataframe["omics"] = dataframe.apply(
        classify_omics,
        axis=1,
    )

    return dataframe


def load_classification_rules(rules_file):
    """
    Load external biological classification rules.
    """

    rules = pd.read_csv(
        rules_file,
        sep="\t",
        dtype=str,
    ).fillna("")

    required_columns = {
        "rule_id",
        "stage",
        "field",
        "pattern",
        "match_type",
        "classification",
        "priority",
        "evidence",
        "reference",
    }

    missing = required_columns - set(rules.columns)

    if missing:
        raise ValueError(
            f"Missing required columns in classification rules: {missing}"
        )

    rules["priority"] = pd.to_numeric(
        rules["priority"],
        errors="coerce",
    ).fillna(0).astype(int)

    return rules

def build_rule_mask(metadata, rule):
    """
    Return a boolean mask identifying records matching one rule.
    """

    field = rule["field"]
    pattern = str(rule["pattern"]).strip()
    match_type = str(rule["match_type"]).strip().lower()

    if field not in metadata.columns:
        raise ValueError(
            f"Field '{field}' from rule {rule['rule_id']} "
            "is not present in the metadata."
        )

    values = (
        metadata[field]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    if match_type == "exact":
        mask = values.str.casefold().eq(
            pattern.casefold()
        )

    elif match_type == "prefix":
        mask = (
            values
            .str.casefold()
            .str.startswith(pattern.casefold())
        )

    elif match_type == "contains":
        mask = values.str.contains(
            pattern,
            case=False,
            regex=False,
            na=False,
        )

    elif match_type == "regex":
        mask = values.str.contains(
            pattern,
            case=False,
            regex=True,
            na=False,
        )

    else:
        raise ValueError(
            f"Unsupported match_type '{match_type}' "
            f"in rule {rule['rule_id']}."
        )

    return mask

def apply_explicit_evidence_rules(metadata, rules):
    """
    Apply rules belonging to the explicit_evidence stage.
    """

    metadata = metadata.copy()

    metadata["target_evidence"] = False
    metadata["target_status"] = "unresolved"

    stage_rules = (
        rules[
            rules["stage"] == "explicit_evidence"
        ]
        .sort_values(
            "priority",
            ascending=False,
            kind="stable",
        )
    )

    for _, rule in stage_rules.iterrows():

        mask = build_rule_mask(
            metadata,
            rule,
        )

        metadata.loc[
            mask,
            "target_evidence"
        ] = True

        metadata.loc[
            mask,
            "target_status"
        ] = rule["classification"]

    return metadata


def apply_classification_stage(
    metadata,
    rules,
    stage,
):
    """
    Apply one stage of external biological classification rules
    only to records that remain unresolved.
    """

    metadata = metadata.copy()

    stage_rules = (
        rules[
            rules["stage"] == stage
        ]
        .sort_values(
            "priority",
            ascending=False,
            kind="stable",
        )
    )

    for _, rule in stage_rules.iterrows():

        unresolved_mask = (
            metadata["target_status"] == "unresolved"
        )

        match_mask = build_rule_mask(
            metadata,
            rule,
        )

        metadata.loc[
            unresolved_mask & match_mask,
            "target_status"
        ] = rule["classification"]

    return metadata

def load_curated_studies(curation_file):
    """Load study-specific curation decisions from a TSV file."""

    curation = pd.read_csv(
        curation_file,
        sep="\t",
        dtype=str,
    )

    required_columns = {
        "curation_id",
        "study_accession",
        "field",
        "pattern",
        "match_type",
        "classification",
        "priority",
        "evidence",
        "reference",
    }

    missing_columns = (
        required_columns
        - set(curation.columns)
    )

    if missing_columns:
        raise ValueError(
            "Missing required columns in curated studies file: "
            + ", ".join(sorted(missing_columns))
        )

    curation["priority"] = (
        pd.to_numeric(
            curation["priority"],
            errors="coerce",
        )
        .fillna(0)
        .astype(int)
    )

    curation = curation.sort_values(
        "priority",
        ascending=False,
    )

    return curation

def apply_study_curation(metadata, curation_rules):
    """
    Apply study-specific curation decisions only to records
    that remain unresolved.
    """

    metadata = metadata.copy()

    for _, rule in curation_rules.iterrows():

        unresolved_mask = (
            metadata["target_status"]
            == "unresolved"
        )

        study_mask = (
            metadata["study_accession"]
            .fillna("")
            .astype(str)
            .eq(
                str(
                    rule["study_accession"]
                )
            )
        )

        rule_mask = build_rule_mask(
            metadata,
            rule,
        )

        final_mask = (
            unresolved_mask
            & study_mask
            & rule_mask
        )

        metadata.loc[
            final_mask,
            "target_status",
        ] = rule["classification"]

    return metadata


# ============================================================
# Output
# ============================================================

def save_metadata_table(dataframe, output_file, label):
    """
    Save an intermediate metadata table as TSV.
    """

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    dataframe.to_csv(
        output_file,
        sep="\t",
        index=False,
    )

    print(
        f"[INFO] {label} written to: {output_file} "
        f"({len(dataframe)} records)"
    )


def build_samples_table(dataframe, target_omics):
    """
    Build the final sample table for downstream processing.

    Only records belonging to the selected omics and finally
    classified as target are retained.
    """

    samples = dataframe[
        dataframe["omics"].isin(target_omics)
        & (dataframe["target_status"] == "target")
        & (dataframe["selection_status"] == "retained_by_rules")
    ].copy()

    return samples


def save_samples_table(dataframe, output_file):
    """Write final samples table as TSV."""

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    dataframe.to_csv(
        output_file,
        sep="\t",
        index=False,
    )

    print(
        f"[INFO] samples.tsv written to: {output_file} "
        f"({len(dataframe)} records)"
    )


# ============================================================
# Main
# ============================================================

def main():

    config = load_config(CONFIG_FILE)

    target_name = config["target"]["name"]

    taxid = config["target"]["organism"]["taxid"]

    result_type = config["ena"]["result_type"]

    target_omics = config["omics"]

    raw_output_file = (
         PROJECT_ROOT
         / config["metadata"]["raw_output"]
    )

    fastq_output_file = (
        PROJECT_ROOT
        / config["metadata"]["fastq_output"]
    )

    omics_output_file = (
        PROJECT_ROOT
       / config["metadata"]["omics_output"]
    )

    explicit_evidence_output_file = (
        PROJECT_ROOT
       / config["metadata"]["explicit_evidence_output"]
    )

    cell_line_output_file = (
       PROJECT_ROOT
       / config["metadata"]["cell_line_output"]
    )

    cell_type_output_file = (
       PROJECT_ROOT
       / config["metadata"]["cell_type_output"]
    )

    external_rules_output_file = (
       PROJECT_ROOT
       / config["metadata"]["external_rules_output"]
    )

    curated_studies_file = (
       PROJECT_ROOT
       / config["metadata"]["curated_studies"]
    )

    curation_output_file = (
       PROJECT_ROOT
       / config["metadata"]["curation_output"]
    )

    output_file = (
       PROJECT_ROOT
       / config["metadata"]["output"]
    )

    metadata = fetch_ena_metadata(taxid, result_type)

    save_metadata_table(
         metadata,
         raw_output_file,
         "Raw ENA metadata",
    )

    metadata = filter_fastq_records(metadata)

    save_metadata_table(
         metadata,
         fastq_output_file,
         "FASTQ-filtered metadata",
    )

    metadata = add_omics_classification(metadata)

    save_metadata_table(
        metadata,
        omics_output_file,
        "Omics-classified metadata",
    )


    rules_file = (
        PROJECT_ROOT
        / "config"
        / "classification_rules.tsv"
    )

    rules = load_classification_rules(
        rules_file
    )


# ------------------------------------------------------------
# Explicit textual evidence
# ------------------------------------------------------------

    metadata = apply_explicit_evidence_rules(
        metadata,
        rules,
    )

    save_metadata_table(
        metadata,
        explicit_evidence_output_file,
        "Explicit-evidence metadata",
    )


# ------------------------------------------------------------
# cell_line classification
# ------------------------------------------------------------

    metadata = apply_classification_stage(
        metadata,
        rules,
        stage="cell_line",
    )

    save_metadata_table(
        metadata,
        cell_line_output_file,
        "Metadata after cell_line classification",
    )


# ------------------------------------------------------------
# cell_type classification
# ------------------------------------------------------------

    metadata = apply_classification_stage(
        metadata,
        rules,
        stage="cell_type",
    )

    save_metadata_table(
        metadata,
        cell_type_output_file,
        "Metadata after cell_type classification",
    )


# ------------------------------------------------------------
# Reuse of previously verified external rules
# ------------------------------------------------------------

    metadata = apply_classification_stage(
        metadata,
        rules,
        stage="external",
    )

    save_metadata_table(
        metadata,
        external_rules_output_file,
        "Metadata after external classification rules",
    )

# ------------------------------------------------------------
# Study-specific curation
# ------------------------------------------------------------

    curation_rules = load_curated_studies(
        curated_studies_file
    )

    metadata = apply_study_curation(
        metadata,
        curation_rules,
    )

    # Previous stage tables remain legacy diagnostic evidence. The final table
    # independently resolves sample identity, conflicts and protocol eligibility.
    metadata = reclassify_dataframe(
        metadata, rules, curation_rules,
        PROJECT_ROOT / "config" / "selection_decisions.tsv", target_omics,
    )

    save_workflow_audit(
        metadata,
        PROJECT_ROOT / config["metadata"].get("selection_audit_dir", "results/metadata/selection_audit"),
        {
            "fastq_metadata": fastq_output_file,
            "classification_rules": rules_file,
            "study_curation": curated_studies_file,
            "run_decisions": PROJECT_ROOT / "config/selection_decisions.tsv",
            "policy_code": PROJECT_ROOT / "workflow/scripts/selection_policy.py",
            "config": PROJECT_ROOT / "config/config.yaml",
        },
        target_omics,
    )

    save_metadata_table(
        metadata,
        curation_output_file,
        "Metadata after study-specific curation",
    )


# ------------------------------------------------------------
# Final samples table for downstream processing
# ------------------------------------------------------------

    samples = build_samples_table(
        metadata,
        target_omics,
    )

    save_samples_table(
        samples,
        output_file,
    )
    target = metadata[
        metadata["omics"].isin(target_omics)
    ]

    print(
       f"\n[INFO] Final classification for target "
       f"'{target_name}' after study-specific curation:"
    )

    print(
       target["target_status"]
       .value_counts(dropna=False)
       .to_string()
    )

    print("\n[INFO] Omics classification:")

    print(
    metadata["omics"]
    .value_counts(dropna=False)
    .to_string()
    )


if __name__ == "__main__":

    try:
        main()

    except Exception as error:

        print(
            f"[ERROR] {error}",
            file=sys.stderr,
        )

        sys.exit(1)
