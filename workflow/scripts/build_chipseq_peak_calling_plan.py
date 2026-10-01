#!/usr/bin/env python3
"""Build a deterministic ChIP-seq peak-calling plan from ready IP/Input analyses.

The upstream analysis plan determines which IP/Input analyses are biologically
eligible and unambiguously paired. This layer determines how each ready
analysis should be peak-called.

Unsupported studies, targets, layouts, controls or fragment-size policies
fail closed.
"""

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


def sha256(path):
    path = Path(path)
    h = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)

    return h.hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path):
    path = Path(path)

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []

    if not fields:
        raise ValueError(f"Missing TSV header: {path}")

    return rows, fields


def atomic_write(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(path)


def tsv_text(fields, rows):
    from io import StringIO

    buffer = StringIO(newline="")

    writer = csv.DictWriter(
        buffer,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )

    writer.writeheader()
    writer.writerows(rows)

    return buffer.getvalue()


def require_nonempty(row, field, label):
    value = (row.get(field) or "").strip()

    if not value:
        raise ValueError(f"Missing {field} in {label}")

    return value


SAFE_COMPONENT = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")
SUPPORTED_GENOME_SIZE_POLICY = (
    "derive_nuclear_reference_span_from_fai_"
    "excluding_mitochondrial_contig"
)


def require_safe_identifier(value, label, pattern=SAFE_COMPONENT):
    if (
        not isinstance(value, str)
        or not pattern.fullmatch(value)
        or value in {".", ".."}
    ):
        raise ValueError(f"Unsafe {label} for path use: {value!r}")
    return value


def require_boolean(mapping, field, label):
    value = mapping.get(field)
    if type(value) is not bool:
        raise ValueError(f"{label}.{field} must be a boolean")
    return value


def require_policy_text(mapping, field, label):
    value = mapping.get(field)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Missing {label}.{field}")
    return value


def validate_policy(policy):
    if policy.get("schema_version") != 1:
        raise ValueError("Unsupported ChIP peak-calling policy schema")

    if policy.get("scope") != "dynamic_chipseq_peak_calling_policy":
        raise ValueError("Unexpected ChIP peak-calling policy scope")

    defaults = policy.get("technical_defaults")

    if not isinstance(defaults, dict):
        raise ValueError(
            "Missing technical_defaults in peak-calling policy"
        )

    layouts = defaults.get("supported_library_layouts")

    if not isinstance(layouts, list) or not layouts:
        raise ValueError(
            "No supported library layouts configured"
        )

    qvalue = float(defaults.get("qvalue", 0))

    if not 0 < qvalue < 1:
        raise ValueError(
            "qvalue must be between 0 and 1"
        )

    broad_cutoff = float(defaults.get("broad_cutoff", 0))

    if not 0 < broad_cutoff < 1:
        raise ValueError(
            "broad_cutoff must be between 0 and 1"
        )

    if defaults.get("keep_dup") != "all":
        raise ValueError(
            "Upstream-deduplicated BAMs require keep_dup=all"
        )

    if defaults.get("format") != "BAM":
        raise ValueError(
            "Current dynamic ChIP-seq policy supports BAM input"
        )

    if defaults.get("scale_to") not in {"small", "large"}:
        raise ValueError(
            "scale_to must be small or large"
        )

    for boolean_field in (
        "spmr",
        "store_bdg",
        "call_summits",
        "cutoff_analysis",
    ):
        require_boolean(defaults, boolean_field, "technical_defaults")

    if defaults.get("spmr") is not True:
        raise ValueError(
            "Current harmonized policy requires SPMR"
        )

    if defaults.get("store_bdg") is not True:
        raise ValueError(
            "Current harmonized policy requires bedGraph outputs"
        )

    if defaults.get("mitochondrial_accession") != "NC_007936.1":
        raise ValueError("Unsupported mitochondrial accession")

    if (
        defaults.get("effective_genome_size_policy")
        != SUPPORTED_GENOME_SIZE_POLICY
    ):
        raise ValueError("Unsupported effective genome-size policy")

    require_policy_text(defaults, "parameter_basis", "technical_defaults")

    studies = policy.get("studies")

    if not isinstance(studies, dict) or not studies:
        raise ValueError(
            "No study-specific peak-calling policies configured"
        )

    for study, study_policy in studies.items():
        require_safe_identifier(study, "study identifier")
        require_policy_text(study_policy, "evidence_type", study)
        require_policy_text(study_policy, "source_doi", study)
        require_policy_text(study_policy, "source_url", study)
        targets = study_policy.get("targets")
        fragment = study_policy.get("fragment_size_policy")

        if not isinstance(targets, dict) or not targets:
            raise ValueError(
                f"No target policies configured for {study}"
            )

        if not isinstance(fragment, dict):
            raise ValueError(
                f"Missing fragment_size_policy for {study}"
            )

        require_policy_text(fragment, "basis", f"{study}.fragment_size_policy")

        mode = fragment.get("mode")

        if mode == "fixed":
            if type(fragment.get("extsize")) is not int or fragment["extsize"] <= 0:
                raise ValueError(
                    f"Fixed fragment policy needs positive extsize: "
                    f"{study}"
                )

        elif mode == "phantompeakqualtools":
            require_boolean(
                fragment,
                "requires_external_estimation",
                f"{study}.fragment_size_policy",
            )
            if fragment["requires_external_estimation"] is not True:
                raise ValueError(
                    "PhantomPeakQualTools policy must require "
                    f"external estimation: {study}"
                )

        else:
            raise ValueError(
                f"Unsupported fragment-size policy for {study}: "
                f"{mode}"
            )

        for target, target_policy in targets.items():
            require_safe_identifier(target, f"target identifier for {study}")
            if not isinstance(target_policy, dict):
                raise ValueError(f"Invalid target policy for {study}/{target}")
            peak_mode = target_policy.get("peak_mode")

            if peak_mode not in {"narrow", "broad"}:
                raise ValueError(
                    f"Unsupported peak mode for "
                    f"{study}/{target}: {peak_mode}"
                )

    return defaults, studies


def build_plan(analysis_plan_path, policy_path):
    rows, fields = read_tsv(analysis_plan_path)

    policy = load_json(policy_path)
    defaults, studies = validate_policy(policy)

    required_fields = {
        "analysis_id",
        "study_accession",
        "ip_run_accession",
        "control_run_accession",
        "declared_target",
        "library_layout",
        "analysis_status",
    }

    missing = sorted(required_fields - set(fields))

    if missing:
        raise ValueError(
            f"Analysis plan missing required columns: {missing}"
        )

    supported_layouts = set(
        defaults["supported_library_layouts"]
    )

    seen_analysis_ids = set()
    seen_ip_runs = set()
    peak_plan = []

    for row in rows:
        if row.get("analysis_status") != "ready":
            continue

        analysis_id = require_nonempty(
            row,
            "analysis_id",
            "ready analysis",
        )
        require_safe_identifier(analysis_id, "analysis_id")

        if analysis_id in seen_analysis_ids:
            raise ValueError(
                f"Duplicate ready analysis_id: {analysis_id}"
            )

        seen_analysis_ids.add(analysis_id)

        study = require_nonempty(
            row,
            "study_accession",
            analysis_id,
        )
        require_safe_identifier(study, "study_accession")

        ip_run = require_nonempty(
            row,
            "ip_run_accession",
            analysis_id,
        )
        require_safe_identifier(ip_run, "IP run accession")

        if ip_run in seen_ip_runs:
            raise ValueError(f"Duplicate ready IP analysis: {ip_run}")
        seen_ip_runs.add(ip_run)

        control_run = require_nonempty(
            row,
            "control_run_accession",
            analysis_id,
        )
        require_safe_identifier(control_run, "Input run accession")

        target = require_nonempty(
            row,
            "declared_target",
            analysis_id,
        )
        require_safe_identifier(target, "declared target")

        layout = require_nonempty(
            row,
            "library_layout",
            analysis_id,
        )

        if ip_run == control_run:
            raise ValueError(
                f"IP and Input are identical for {analysis_id}"
            )

        if layout not in supported_layouts:
            raise ValueError(
                "Unsupported library layout for ready analysis "
                f"{analysis_id}: {layout}"
            )

        if study not in studies:
            raise ValueError(
                "No peak-calling study policy for ready analysis "
                f"{analysis_id}: {study}"
            )

        study_policy = studies[study]

        target_policy = study_policy["targets"].get(target)

        if target_policy is None:
            raise ValueError(
                "No peak-calling target policy for ready analysis "
                f"{analysis_id}: {study}/{target}"
            )

        peak_mode = target_policy["peak_mode"]

        fragment = study_policy["fragment_size_policy"]
        fragment_mode = fragment["mode"]

        if fragment_mode == "fixed":
            fragment_size_bp = str(
                int(fragment["extsize"])
            )

            fragment_size_status = "resolved"
            execution_status = "ready_for_peak_calling"
            macs3_model = "fixed"

        elif fragment_mode == "phantompeakqualtools":
            fragment_size_bp = ""

            fragment_size_status = "requires_estimation"
            execution_status = "requires_fragment_estimation"
            macs3_model = "external_fixed_after_estimation"

        else:
            raise ValueError(
                "Unsupported fragment-size mode for ready "
                f"analysis {analysis_id}: {fragment_mode}"
            )

        peak_plan.append({
            "analysis_id": analysis_id,
            "study_accession": study,
            "ip_run_accession": ip_run,
            "control_run_accession": control_run,
            "declared_target": target,
            "library_layout": layout,
            "peak_mode": peak_mode,
            "format": defaults["format"],
            "qvalue": str(defaults["qvalue"]),
            "broad_cutoff": (
                str(defaults["broad_cutoff"])
                if peak_mode == "broad"
                else ""
            ),
            "keep_dup": defaults["keep_dup"],
            "scale_to": defaults["scale_to"],
            "spmr": str(defaults["spmr"]).lower(),
            "store_bdg": str(
                defaults["store_bdg"]
            ).lower(),
            "call_summits": str(
                defaults["call_summits"]
            ).lower(),
            "cutoff_analysis": str(
                defaults["cutoff_analysis"]
            ).lower(),
            "mitochondrial_accession": (
                defaults["mitochondrial_accession"]
            ),
            "effective_genome_size_policy": (
                defaults[
                    "effective_genome_size_policy"
                ]
            ),
            "fragment_size_policy": fragment_mode,
            "fragment_size_bp": fragment_size_bp,
            "fragment_size_status": fragment_size_status,
            "macs3_model": macs3_model,
            "fragment_size_basis": fragment["basis"],
            "execution_status": execution_status,
            "evidence_type": study_policy[
                "evidence_type"
            ],
            "source_doi": study_policy["source_doi"],
            "source_url": study_policy["source_url"],
            "evidence_note": target_policy.get(
                "evidence_note",
                "",
            ),
        })

    return peak_plan


FIELDS = [
    "analysis_id",
    "study_accession",
    "ip_run_accession",
    "control_run_accession",
    "declared_target",
    "library_layout",
    "peak_mode",
    "format",
    "qvalue",
    "broad_cutoff",
    "keep_dup",
    "scale_to",
    "spmr",
    "store_bdg",
    "call_summits",
    "cutoff_analysis",
    "mitochondrial_accession",
    "effective_genome_size_policy",
    "fragment_size_policy",
    "fragment_size_bp",
    "fragment_size_status",
    "macs3_model",
    "fragment_size_basis",
    "execution_status",
    "evidence_type",
    "source_doi",
    "evidence_note",
    "source_url",
]


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--analysis-plan",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--policy",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
    )

    args = parser.parse_args()

    plan = build_plan(
        args.analysis_plan,
        args.policy,
    )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    plan_path = (
        args.outdir /
        "chipseq_peak_calling_plan.tsv"
    )

    atomic_write(
        plan_path,
        tsv_text(FIELDS, plan),
    )

    status_counts = Counter(
        row["execution_status"]
        for row in plan
    )

    study_counts = Counter(
        row["study_accession"]
        for row in plan
    )

    target_counts = Counter(
        row["declared_target"]
        for row in plan
    )

    peak_mode_counts = Counter(
        row["peak_mode"]
        for row in plan
    )

    fragment_policy_counts = Counter(
        row["fragment_size_policy"]
        for row in plan
    )

    summary = {
        "schema_version": 1,
        "scope": "dynamic_chipseq_peak_calling_plan",
        "peak_analyses": len(plan),
        "execution_status_counts": dict(
            sorted(status_counts.items())
        ),
        "study_counts": dict(
            sorted(study_counts.items())
        ),
        "target_counts": dict(
            sorted(target_counts.items())
        ),
        "peak_mode_counts": dict(
            sorted(peak_mode_counts.items())
        ),
        "fragment_policy_counts": dict(
            sorted(fragment_policy_counts.items())
        ),
        "input_sha256": {
            "analysis_plan": sha256(
                args.analysis_plan
            ),
            "policy": sha256(
                args.policy
            ),
        },
        "peak_calling_plan_sha256": sha256(
            plan_path
        ),
    }

    summary_path = (
        args.outdir /
        "chipseq_peak_calling_plan_summary.json"
    )

    atomic_write(
        summary_path,
        json.dumps(
            summary,
            indent=2,
            sort_keys=True,
        ) + "\n",
    )

    print(
        "[OK] Dynamic ChIP-seq peak-calling plan "
        f"built: {len(plan)} analyses"
    )

    for key, count in sorted(
        status_counts.items()
    ):
        print(f"{key}: {count}")


if __name__ == "__main__":
    main()
