from pathlib import Path

import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "config.yaml"

METADATA_DIR = PROJECT_ROOT / "results" / "metadata"
REPORT_DIR = METADATA_DIR / "reports"

REPORT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

def load_config(config_file):
    """
    Load pipeline configuration from YAML.
    """

    with open(
        config_file,
        "r",
        encoding="utf-8",
    ) as handle:
        return yaml.safe_load(handle)

def generate_filtering_summary():
    """Summarize ENA retrieval and FASTQ filtering."""

    raw = pd.read_csv(
        METADATA_DIR / "ena_metadata_raw.tsv",
        sep="\t",
    )

    fastq = pd.read_csv(
        METADATA_DIR / "ena_metadata_fastq.tsv",
        sep="\t",
    )

    summary = pd.DataFrame(
        [
            {
                "stage": "ena_retrieval",
                "category": "retrieved",
                "n_records": len(raw),
            },
            {
                "stage": "fastq_filter",
                "category": "retained",
                "n_records": len(fastq),
            },
            {
                "stage": "fastq_filter",
                "category": "excluded",
                "n_records": len(raw) - len(fastq),
            },
        ]
    )

    output = REPORT_DIR / "filtering_summary.tsv"

    summary.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return summary


def generate_omics_summary():
    """Summarize classification by sequencing strategy."""

    metadata = pd.read_csv(
        METADATA_DIR / "ena_metadata_omics.tsv",
        sep="\t",
    )

    summary = (
        metadata["omics"]
        .value_counts(dropna=False)
        .rename_axis("category")
        .reset_index(name="n_records")
    )

    summary.insert(
        0,
        "stage",
        "omics_classification",
    )

    output = REPORT_DIR / "omics_summary.tsv"

    summary.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return summary


def generate_explicit_evidence_summary(target_omics):
    """
    Summarize explicit target evidence in target omics.
    """

    metadata = pd.read_csv(
        METADATA_DIR / "ena_metadata_explicit_evidence.tsv",
        sep="\t",
    )

    target = metadata[
        metadata["omics"].isin(target_omics)
    ].copy()

    rows = [
        {
            "stage": "target_omics",
            "category": "total",
            "n_records": len(target),
        },
        {
            "stage": "explicit_evidence",
            "category": "explicit_target",
            "n_records": int(
                target["target_evidence"].sum()
            ),
        },
        {
            "stage": "explicit_evidence",
            "category": "no_explicit_target",
            "n_records": int(
                (~target["target_evidence"]).sum()
            ),
        },
    ]

    summary = pd.DataFrame(rows)

    output = (
        REPORT_DIR
        / "explicit_evidence_summary.tsv"
    )

    summary.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return summary


def generate_explicit_negative_by_study(target_omics):
    """
    Group records without explicit target evidence by study.
    """

    metadata = pd.read_csv(
        METADATA_DIR / "ena_metadata_explicit_evidence.tsv",
        sep="\t",
    )

    target = metadata[
        metadata["omics"].isin(target_omics)
    ].copy()

    unresolved = target[
        target["target_evidence"] == False
    ].copy()

    summary = (
        unresolved
        .groupby(
            ["study_accession", "study_title"],
            dropna=False,
        )
        .size()
        .reset_index(name="n_records")
        .sort_values(
            "n_records",
            ascending=False,
        )
    )

    output = (
        REPORT_DIR
        / "explicit_negative_by_study.tsv"
    )

    summary.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return summary

def generate_cell_line_values_summary(target_omics):
    """
    Summarize informative cell_line values among records
    without explicit target evidence.
    """

    metadata = pd.read_csv(
        METADATA_DIR / "ena_metadata_explicit_evidence.tsv",
        sep="\t",
    )

    target = metadata[
        metadata["omics"].isin(target_omics)
    ].copy()

    unresolved = target[
       target["target_evidence"] == False
    ].copy()

    unresolved["cell_line"] = (
        unresolved["cell_line"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    informed = unresolved[
        unresolved["cell_line"] != ""
    ].copy()

    summary = (
        informed["cell_line"]
        .value_counts()
        .rename_axis("cell_line")
        .reset_index(name="n_records")
    )

    output = REPORT_DIR / "cell_line_values_summary.tsv"

    summary.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return summary, len(informed), len(unresolved) - len(informed)


def generate_cell_line_classification_summary(target_omics):
    """
    Summarize classification after the cell_line layer
    among records without explicit target evidence.
    """

    metadata = pd.read_csv(
        METADATA_DIR / "ena_metadata_after_cell_line.tsv",
        sep="\t",
    )

    target = metadata[
        metadata["omics"].isin(target_omics)
    ].copy()

    initial_unresolved = target[
       target["target_evidence"] == False
    ].copy()

    summary = (
        initial_unresolved["target_status"]
        .value_counts(dropna=False)
        .rename_axis("classification")
        .reset_index(name="n_records")
    )

    summary["percentage"] = (
        summary["n_records"]
        / len(initial_unresolved)
        * 100
    ).round(2)

    output = (
        REPORT_DIR
        / "cell_line_classification_summary.tsv"
    )

    summary.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return summary

def generate_cell_type_values_summary(target_omics):
    """
    Summarize informative cell_type values among records that
    remained unresolved after the cell_line classification layer.
    """

    metadata = pd.read_csv(
        METADATA_DIR / "ena_metadata_after_cell_line.tsv",
        sep="\t",
    )

    target = metadata[
        metadata["omics"].isin(target_omics)
    ].copy()

    unresolved_after_cell_line = target[
        (target["target_evidence"] == False)
        & (target["target_status"] == "unresolved")
    ].copy()

    unresolved_after_cell_line["cell_type"] = (
        unresolved_after_cell_line["cell_type"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    informed = unresolved_after_cell_line[
        unresolved_after_cell_line["cell_type"] != ""
    ].copy()

    summary = (
        informed["cell_type"]
        .value_counts()
        .rename_axis("cell_type")
        .reset_index(name="n_records")
    )

    summary["percentage"] = (
        summary["n_records"]
        / len(unresolved_after_cell_line)
        * 100
    ).round(2)

    output = REPORT_DIR / "cell_type_values_summary.tsv"

    summary.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return (
        summary,
        len(unresolved_after_cell_line),
        len(informed),
        len(unresolved_after_cell_line) - len(informed),
    )

def generate_cell_type_layer_summary(target_omics):
    """
    Summarize the effect of the cell_type classification layer
    on the 276 records that were unresolved after cell_line.
    """

    before = pd.read_csv(
        METADATA_DIR / "ena_metadata_after_cell_line.tsv",
        sep="\t",
    )

    after = pd.read_csv(
        METADATA_DIR / "ena_metadata_after_cell_type.tsv",
        sep="\t",
    )

    target_before = before[
        before["omics"].isin(target_omics)
    ].copy()

    unresolved_ids = set(
        target_before.loc[
            (target_before["target_evidence"] == False)
            & (target_before["target_status"] == "unresolved"),
            "run_accession",
        ]
    )

    evaluated = after[
        after["run_accession"].isin(unresolved_ids)
    ].copy()

    summary = (
        evaluated["target_status"]
        .value_counts(dropna=False)
        .rename_axis("classification")
        .reset_index(name="n_records")
    )

    summary["percentage"] = (
        summary["n_records"]
        / len(evaluated)
        * 100
    ).round(2)

    output = REPORT_DIR / "cell_type_layer_summary.tsv"

    summary.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return summary

def generate_cell_type_classification_summary(target_omics):
    """
    Summarize accumulated target classification after the
    cell_type layer across all target omics records.
    """

    metadata = pd.read_csv(
        METADATA_DIR / "ena_metadata_after_cell_type.tsv",
        sep="\t",
    )

    target = metadata[
        metadata["omics"].isin(target_omics)
    ].copy()

    summary = (
        target["target_status"]
        .value_counts(dropna=False)
        .rename_axis("classification")
        .reset_index(name="n_records")
    )

    summary["percentage"] = (
        summary["n_records"]
        / len(target)
        * 100
    ).round(2)

    output = (
        REPORT_DIR
        / "cell_type_classification_summary.tsv"
    )

    summary.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return summary

def generate_external_rules_summary(target_omics):
    """
    Summarize accumulated target classification after
    applying external verified classification rules.
    """

    metadata = pd.read_csv(
        METADATA_DIR / "ena_metadata_after_external_rules.tsv",
        sep="\t",
    )

    target = metadata[
        metadata["omics"].isin(target_omics)
    ].copy()

    summary = (
        target["target_status"]
        .value_counts(dropna=False)
        .rename_axis("classification")
        .reset_index(name="n_records")
    )

    summary["percentage"] = (
        summary["n_records"]
        / len(target)
        * 100
    ).round(2)

    output = (
        REPORT_DIR
        / "external_rules_classification_summary.tsv"
    )

    summary.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return summary

def generate_external_rules_comparison(target_omics):
    """
    Compare accumulated target classification before and after
    applying external verified classification rules.
    """

    before = pd.read_csv(
        METADATA_DIR / "ena_metadata_after_cell_type.tsv",
        sep="\t",
    )

    after = pd.read_csv(
        METADATA_DIR / "ena_metadata_after_external_rules.tsv",
        sep="\t",
    )

    target_before = before[
        before["omics"].isin(target_omics)
    ].copy()

    target_after = after[
        after["omics"].isin(target_omics)
    ].copy()

    before_counts = (
        target_before["target_status"]
        .value_counts()
        .reindex(
            ["target", "non_target", "unresolved"],
            fill_value=0,
        )
    )

    after_counts = (
        target_after["target_status"]
        .value_counts()
        .reindex(
            ["target", "non_target", "unresolved"],
            fill_value=0,
        )
    )

    summary = pd.DataFrame(
        {
            "classification": before_counts.index,
            "before_external": before_counts.values,
            "after_external": after_counts.values,
        }
    )

    summary["change"] = (
        summary["after_external"]
        - summary["before_external"]
    )

    output = (
        REPORT_DIR
        / "external_rules_comparison.tsv"
    )

    summary.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return summary

def generate_unresolved_after_external_by_study(target_omics):
    """
    Summarize unresolved records after external rules,
    grouped by study.
    """

    metadata = pd.read_csv(
        METADATA_DIR / "ena_metadata_after_external_rules.tsv",
        sep="\t",
    )

    target = metadata[
        metadata["omics"].isin(target_omics)
    ].copy()

    unresolved = target[
        target["target_status"] == "unresolved"
    ].copy()

    summary = (
        unresolved
        .groupby(
            ["study_accession", "study_title"],
            dropna=False,
        )
        .size()
        .reset_index(name="n_records")
        .sort_values(
            "n_records",
            ascending=False,
        )
    )

    output = (
        REPORT_DIR
        / "unresolved_after_external_by_study.tsv"
    )

    summary.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return summary


def generate_unresolved_after_external_patterns(target_omics):
    """
    Explore recurrent metadata values among records that remain
    unresolved after external classification rules.
    """

    metadata = pd.read_csv(
        METADATA_DIR / "ena_metadata_after_external_rules.tsv",
        sep="\t",
    )

    target = metadata[
        metadata["omics"].isin(target_omics)
    ].copy()

    unresolved = target[
        target["target_status"] == "unresolved"
    ].copy()

    fields = [
        "cell_line",
        "cell_type",
        "sample_title",
        "experiment_title",
        "study_title",
        "library_name",
    ]

    rows = []

    for field in fields:

        values = unresolved[
            [field, "study_accession"]
        ].copy()

        values[field] = (
            values[field]
            .fillna("")
            .astype(str)
            .str.strip()
        )

        values = values[
            values[field] != ""
        ]

        if values.empty:
            continue

        summary = (
            values
            .groupby(field)
            .agg(
                n_records=(field, "size"),
                n_studies=(
                    "study_accession",
                    "nunique",
                ),
            )
            .reset_index()
            .rename(
                columns={
                    field: "value"
                }
            )
        )

        summary.insert(
            0,
            "field",
            field,
        )

        rows.append(summary)

    if rows:
        result = pd.concat(
            rows,
            ignore_index=True,
        )
    else:
        result = pd.DataFrame(
            columns=[
                "field",
                "value",
                "n_records",
                "n_studies",
            ]
        )

    result = result.sort_values(
        [
            "field",
            "n_records",
            "n_studies",
        ],
        ascending=[
            True,
            False,
            False,
        ],
    )

    output = (
        REPORT_DIR
        / "unresolved_after_external_patterns.tsv"
    )

    result.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return result

def generate_curation_comparison(target_omics):
    """
    Compare target classification before and after
    study-specific curation.
    """

    before = pd.read_csv(
        METADATA_DIR / "ena_metadata_after_external_rules.tsv",
        sep="\t",
    )

    after = pd.read_csv(
        METADATA_DIR / "ena_metadata_after_curation.tsv",
        sep="\t",
    )

    target_before = before[
        before["omics"].isin(target_omics)
    ].copy()

    target_after = after[
        after["omics"].isin(target_omics)
    ].copy()

    categories = [
        "target",
        "non_target",
        "unresolved",
    ]

    before_counts = (
        target_before["target_status"]
        .value_counts()
        .reindex(
            categories,
            fill_value=0,
        )
    )

    after_counts = (
        target_after["target_status"]
        .value_counts()
        .reindex(
            categories,
            fill_value=0,
        )
    )

    summary = pd.DataFrame(
        {
            "classification": categories,
            "before_curation": [
                before_counts[category]
                for category in categories
            ],
            "after_curation": [
                after_counts[category]
                for category in categories
            ],
        }
    )

    summary["change"] = (
        summary["after_curation"]
        - summary["before_curation"]
    )

    output = (
        REPORT_DIR
        / "curation_comparison.tsv"
    )

    summary.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return summary

def generate_unresolved_after_curation(target_omics):
    """
    Export records that remain unresolved after
    study-specific curation.
    """

    metadata = pd.read_csv(
        METADATA_DIR / "ena_metadata_after_curation.tsv",
        sep="\t",
    )

    target = metadata[
        metadata["omics"].isin(target_omics)
    ].copy()

    unresolved = target[
        target["target_status"] == "unresolved"
    ].copy()

    columns = [
        "study_accession",
        "run_accession",
        "sample_accession",
        "sample_title",
        "experiment_title",
        "study_title",
        "cell_line",
        "cell_type",
        "omics",
        "target_status",
    ]

    available_columns = [
        column
        for column in columns
        if column in unresolved.columns
    ]

    unresolved = unresolved[
        available_columns
    ].sort_values(
        [
            "study_accession",
            "run_accession",
        ]
    )

    output = (
        REPORT_DIR
        / "unresolved_after_curation.tsv"
    )

    unresolved.to_csv(
        output,
        sep="\t",
        index=False,
    )

    return unresolved


def main():

    config = load_config(CONFIG_FILE)

    target_name = config["target"]["name"]
    target_omics = config["omics"]

    filtering = generate_filtering_summary()
    omics = generate_omics_summary()
    explicit = generate_explicit_evidence_summary(
       target_omics
    )
    studies = generate_explicit_negative_by_study(
       target_omics
    )
    cell_line_values, cell_line_informed, cell_line_missing = (
       generate_cell_line_values_summary(
         target_omics
       )
    )

    cell_line_classification = (
       generate_cell_line_classification_summary(
         target_omics
       )
    )

    cell_type_values, cell_type_total, cell_type_informed, cell_type_missing = (
       generate_cell_type_values_summary(
         target_omics
       )
    )

    cell_type_layer = (
       generate_cell_type_layer_summary(
         target_omics
       )
    )

    cell_type_classification = (
       generate_cell_type_classification_summary(
         target_omics
       )
    )

    external_rules_classification = (
       generate_external_rules_summary(
         target_omics
       )
    )

    external_rules_comparison = (
       generate_external_rules_comparison(
         target_omics
       )
    )

    curation_comparison = (
       generate_curation_comparison(
         target_omics
       )
    )

    unresolved_after_curation = (
       generate_unresolved_after_curation(
         target_omics
       )
    )

    unresolved_by_study = (
       generate_unresolved_after_external_by_study(
         target_omics
       )
    )

    unresolved_patterns = (
       generate_unresolved_after_external_patterns(
         target_omics
       )
    )

    print("\n[REPORT] Filtering summary")
    print(filtering.to_string(index=False))

    print("\n[REPORT] Omics summary")
    print(omics.to_string(index=False))

    print(
       f"\n[REPORT] Explicit evidence for target "
       f"'{target_name}'"
    )
    print(explicit.to_string(index=False))
    print(
       f"\n[REPORT] Records without explicit evidence "
       f"for target '{target_name}'"
    )
    print(f"Records: {studies['n_records'].sum()}")
    print(f"Studies: {studies['study_accession'].nunique()}")

    print("\n[REPORT] cell_line exploration")
    print(f"Distinct cell_line values: {len(cell_line_values)}")
    print(f"Records with cell_line: {cell_line_informed}")
    print(f"Records without cell_line: {cell_line_missing}")

    print("\n[REPORT] cell_line values")
    print(cell_line_values.to_string(index=False))

    print("\n[REPORT] Classification after cell_line")
    print(cell_line_classification.to_string(index=False))

    print("\n[REPORT] cell_type exploration")
    print(f"Records evaluated after cell_line: {cell_type_total}")
    print(f"Distinct cell_type values: {len(cell_type_values)}")
    print(f"Records with cell_type: {cell_type_informed}")
    print(f"Records without cell_type: {cell_type_missing}")

    print("\n[REPORT] cell_type values")
    print(cell_type_values.to_string(index=False))

    print("\n[REPORT] Classification within cell_type layer")
    print(cell_type_layer.to_string(index=False))

    print("\n[REPORT] Accumulated classification after cell_type")
    print(cell_type_classification.to_string(index=False))

    print("\n[REPORT] Accumulated classification after external rules")
    print(external_rules_classification.to_string(index=False))

    print("\n[REPORT] Change after external rules")
    print(external_rules_comparison.to_string(index=False))

    print("\n[REPORT] Change after study-specific curation")
    print(curation_comparison.to_string(index=False))

    print("\n[REPORT] Unresolved after study-specific curation")
    print(f"Records: {len(unresolved_after_curation)}")
    print(
        f"Studies: "
        f"{unresolved_after_curation['study_accession'].nunique()}"
    )

    print(
        unresolved_after_curation.to_string(
          index=False
        )
    )

    print("\n[REPORT] Unresolved after external rules")
    print(
       f"Records: "
       f"{unresolved_by_study['n_records'].sum()}"
    )
    print(
       f"Studies: "
       f"{unresolved_by_study['study_accession'].nunique()}"
    )

    print("\n[REPORT] Unresolved records by study")
    print(
       unresolved_by_study.to_string(index=False)
    )

    print("\n[REPORT] Top unresolved metadata patterns")

    for field in [
       "cell_line",
       "cell_type",
       "sample_title",
       "experiment_title",
       "study_title",
       "library_name",
    ]:
       subset = unresolved_patterns[
         unresolved_patterns["field"] == field
       ].head(10)

       print(f"\n{field}:")

       if subset.empty:
         print("No informative values")
       else:
         print(
            subset.to_string(index=False)
         )

    print(
        f"\n[INFO] Reports written to: {REPORT_DIR}"
    )


if __name__ == "__main__":
    main()
