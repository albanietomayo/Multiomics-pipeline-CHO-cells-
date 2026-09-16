# ============================================================
# Metadata retrieval, classification and curation
# ============================================================

rule metadata:
    input:
        classification_rules="config/classification_rules.tsv",
        curated_studies=config["metadata"]["curated_studies"],
        selection_policy="workflow/scripts/selection_policy.py",
        selection_decisions="config/selection_decisions.tsv",
        script="workflow/scripts/fetch_ena_metadata.py"

    output:
        raw=config["metadata"]["raw_output"],
        fastq=config["metadata"]["fastq_output"],
        omics=config["metadata"]["omics_output"],
        explicit=config["metadata"]["explicit_evidence_output"],
        cell_line=config["metadata"]["cell_line_output"],
        cell_type=config["metadata"]["cell_type_output"],
        external_rules=config["metadata"]["external_rules_output"],
        curated=config["metadata"]["curation_output"],
        samples=config["metadata"]["output"],
        selection_summary=config["metadata"].get("selection_audit_dir", "results/metadata/selection_audit") + "/summary.json",
        selection_excluded=config["metadata"].get("selection_audit_dir", "results/metadata/selection_audit") + "/excluded.tsv",
        selection_review=config["metadata"].get("selection_audit_dir", "results/metadata/selection_audit") + "/review_required.tsv"

    params:
        target_name=config["target"]["name"],
        taxid=str(config["target"]["organism"]["taxid"]),
        result_type=config["ena"]["result_type"],
        omics=",".join(config["omics"])

    shell:
        """
        python {input.script}
        """


# ============================================================
# Metadata audit reports
# ============================================================

rule metadata_reports:
    input:
        raw="results/metadata/ena_metadata_raw.tsv",
        fastq="results/metadata/ena_metadata_fastq.tsv",
        omics="results/metadata/ena_metadata_omics.tsv",
        explicit="results/metadata/ena_metadata_explicit_evidence.tsv",
        cell_line="results/metadata/ena_metadata_after_cell_line.tsv",
        cell_type="results/metadata/ena_metadata_after_cell_type.tsv",
        external_rules="results/metadata/ena_metadata_after_external_rules.tsv",
        curated="results/metadata/ena_metadata_after_curation.tsv",
        script="workflow/scripts/generate_metadata_report.py"

    output:
        filtering="results/metadata/reports/filtering_summary.tsv",
        omics_summary="results/metadata/reports/omics_summary.tsv",
        explicit_summary="results/metadata/reports/explicit_evidence_summary.tsv",
        explicit_by_study="results/metadata/reports/explicit_negative_by_study.tsv",
        cell_line_values="results/metadata/reports/cell_line_values_summary.tsv",
        cell_line_classification="results/metadata/reports/cell_line_classification_summary.tsv",
        cell_type_values="results/metadata/reports/cell_type_values_summary.tsv",
        cell_type_layer="results/metadata/reports/cell_type_layer_summary.tsv",
        cell_type_classification="results/metadata/reports/cell_type_classification_summary.tsv",
        external_classification="results/metadata/reports/external_rules_classification_summary.tsv",
        external_comparison="results/metadata/reports/external_rules_comparison.tsv",
        unresolved_by_study="results/metadata/reports/unresolved_after_external_by_study.tsv",
        unresolved_patterns="results/metadata/reports/unresolved_after_external_patterns.tsv",
        curation_comparison="results/metadata/reports/curation_comparison.tsv",
        unresolved_after_curation="results/metadata/reports/unresolved_after_curation.tsv"

    shell:
        """
        python {input.script}
        """
