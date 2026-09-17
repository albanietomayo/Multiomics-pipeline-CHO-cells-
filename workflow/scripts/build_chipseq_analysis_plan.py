#!/usr/bin/env python3
"""Build a deterministic, fail-closed ChIP-seq IP/Input analysis plan.

This layer converts the existing ChIP metadata, annotation, condition,
control-candidate and experimental-eligibility evidence into a downstream
analysis plan.

It does not infer biological relationships when evidence is ambiguous.
Only a unique, exact-condition, eligible Input may be promoted automatically.
"""

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


CONDITION_FIELDS = (
    "study_accession",
    "condition_rule",
    "label_context",
    "time_value",
    "time_unit",
)


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


def index_rows(rows, key, label):
    result = {}

    for row in rows:
        value = (row.get(key) or "").strip()

        if not value:
            raise ValueError(f"Missing {key} in {label}")

        if value in result:
            raise ValueError(f"Duplicate {key} in {label}: {value}")

        result[value] = row

    return result


def parse_json_list(value, label):
    try:
        parsed = json.loads(value or "[]")
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON list in {label}: {value}") from exc

    if not isinstance(parsed, list):
        raise ValueError(f"Expected JSON list in {label}")

    if len(parsed) != len(set(parsed)):
        raise ValueError(f"Duplicate accessions in {label}")

    return parsed


def condition_key(row, fields):
    return tuple((row.get(field) or "").strip() for field in fields)


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


def build_plan(
    runs_path,
    annotations_path,
    conditions_path,
    candidates_path,
    eligibility_path,
    policy_path,
):
    policy = load_json(policy_path)

    if policy.get("schema_version") != 1:
        raise ValueError("Unsupported ChIP analysis policy schema")

    if policy.get("scope") != "incremental_chipseq_analysis_planning":
        raise ValueError("Unexpected ChIP analysis policy scope")

    match_fields = tuple(policy.get("condition_match_fields", []))

    if match_fields != CONDITION_FIELDS:
        raise ValueError(
            "Condition match fields differ from the validated exact-match policy"
        )

    supported_layouts = set(policy.get("supported_library_layouts", []))

    if not supported_layouts:
        raise ValueError("No supported ChIP-seq library layouts configured")

    runs_rows, _ = read_tsv(runs_path)
    annotation_rows, _ = read_tsv(annotations_path)
    condition_rows, _ = read_tsv(conditions_path)
    candidate_rows, _ = read_tsv(candidates_path)

    eligibility = load_json(eligibility_path)

    if eligibility.get("schema_version") != 1:
        raise ValueError("Unsupported experimental eligibility schema")

    eligibility_rows = eligibility.get("runs")

    if not isinstance(eligibility_rows, list):
        raise ValueError("Experimental eligibility JSON has no runs list")

    runs = index_rows(runs_rows, "run_accession", "run inventory")
    annotations = index_rows(
        annotation_rows,
        "run_accession",
        "target annotations",
    )
    conditions = index_rows(
        condition_rows,
        "run_accession",
        "condition table",
    )
    decisions = index_rows(
        eligibility_rows,
        "run_accession",
        "experimental eligibility",
    )

    run_set = set(runs)

    for label, collection in (
        ("target annotations", annotations),
        ("conditions", conditions),
        ("experimental eligibility", decisions),
    ):
        if set(collection) != run_set:
            missing = sorted(run_set - set(collection))
            extra = sorted(set(collection) - run_set)
            raise ValueError(
                f"{label} does not match active ChIP inventory; "
                f"missing={missing}, extra={extra}"
            )

    ip_runs = {
        run
        for run, row in annotations.items()
        if row.get("library_role") == "ip"
    }

    input_runs = {
        run
        for run, row in annotations.items()
        if row.get("library_role") == "input"
    }

    candidate_by_ip = index_rows(
        candidate_rows,
        "ip_run_accession",
        "control-candidate table",
    )

    if set(candidate_by_ip) != ip_runs:
        missing = sorted(ip_runs - set(candidate_by_ip))
        extra = sorted(set(candidate_by_ip) - ip_runs)
        raise ValueError(
            "Control-candidate table does not cover exactly the IP cohort; "
            f"missing={missing}, extra={extra}"
        )

    plan = []

    for ip_run in sorted(ip_runs):
        inventory = runs[ip_run]
        annotation = annotations[ip_run]
        condition = conditions[ip_run]
        decision = decisions[ip_run]
        candidate = candidate_by_ip[ip_run]

        study = annotation["study_accession"]
        target = annotation["declared_target"]
        layout = inventory["library_layout"]

        if not target:
            raise ValueError(f"IP has no declared target: {ip_run}")

        for source_name, source_row in (
            ("inventory", inventory),
            ("condition", condition),
            ("eligibility", decision),
            ("candidate", candidate),
        ):
            if source_row.get("study_accession") != study:
                raise ValueError(
                    f"Study mismatch for {ip_run} in {source_name}"
                )

        if condition.get("library_role") != "ip":
            raise ValueError(f"Condition role mismatch for IP: {ip_run}")

        if decision.get("library_role") != "ip":
            raise ValueError(f"Eligibility role mismatch for IP: {ip_run}")

        if decision.get("declared_target") != target:
            raise ValueError(
                f"Eligibility target mismatch for IP: {ip_run}"
            )

        candidate_runs = parse_json_list(
            candidate.get("candidate_input_runs", "[]"),
            f"{ip_run}/candidate_input_runs",
        )

        issue_runs = parse_json_list(
            candidate.get("candidates_with_documented_issue", "[]"),
            f"{ip_run}/candidates_with_documented_issue",
        )

        try:
            declared_count = int(candidate.get("candidate_count", ""))
        except ValueError as exc:
            raise ValueError(
                f"Invalid candidate_count for {ip_run}"
            ) from exc

        if declared_count != len(candidate_runs):
            raise ValueError(
                f"candidate_count disagrees with candidate_input_runs: {ip_run}"
            )

        candidate_status = candidate.get("candidate_status", "")
        eligibility_status = decision.get(
            "experimental_eligibility",
            "",
        )

        analysis_status = ""
        reason = ""
        control_run = ""
        pairing_rule = ""

        if eligibility_status == "excluded":
            analysis_status = "excluded"
            reason = "ip_excluded_by_experimental_eligibility"

        elif eligibility_status == "review_required":
            analysis_status = "review_required"
            reason = "ip_experimental_eligibility_requires_review"

        elif eligibility_status != "retained_by_rules":
            raise ValueError(
                f"Unknown eligibility state for {ip_run}: "
                f"{eligibility_status}"
            )

        elif condition.get("condition_status") != "resolved_from_labels":
            analysis_status = "review_required"
            reason = "ip_condition_not_resolved"

        elif layout not in supported_layouts:
            analysis_status = "review_required"
            reason = "library_layout_not_yet_supported"

        elif candidate_status == "no_candidate":
            if candidate_runs:
                raise ValueError(
                    f"no_candidate row contains candidates: {ip_run}"
                )

            analysis_status = "blocked_no_control"
            reason = "no_compatible_input_control"

        elif candidate_status == "multiple_candidates_needs_review":
            if len(candidate_runs) < 2:
                raise ValueError(
                    f"multiple-candidate status inconsistent for {ip_run}"
                )

            analysis_status = "review_required"
            reason = "multiple_compatible_input_controls"

        elif candidate_status == "single_candidate_needs_review":
            if len(candidate_runs) != 1:
                raise ValueError(
                    f"single-candidate status inconsistent for {ip_run}"
                )

            proposed_control = candidate_runs[0]

            if proposed_control not in input_runs:
                raise ValueError(
                    f"Candidate is not an annotated Input: "
                    f"{ip_run} -> {proposed_control}"
                )

            control_annotation = annotations[proposed_control]
            control_condition = conditions[proposed_control]
            control_decision = decisions[proposed_control]
            control_inventory = runs[proposed_control]

            if proposed_control in issue_runs:
                analysis_status = "review_required"
                reason = "unique_input_has_documented_issue"

            elif control_annotation.get("study_accession") != study:
                raise ValueError(
                    f"Cross-study candidate pairing detected: "
                    f"{ip_run} -> {proposed_control}"
                )

            elif control_condition.get("library_role") != "input":
                raise ValueError(
                    f"Candidate condition role is not Input: "
                    f"{proposed_control}"
                )

            elif control_decision.get(
                "experimental_eligibility"
            ) != "retained_by_rules":
                analysis_status = "review_required"
                reason = "unique_input_not_retained_by_eligibility"

            elif control_inventory.get(
                "library_layout"
            ) not in supported_layouts:
                analysis_status = "review_required"
                reason = "control_library_layout_not_yet_supported"

            elif condition_key(
                condition,
                match_fields,
            ) != condition_key(
                control_condition,
                match_fields,
            ):
                raise ValueError(
                    f"Candidate does not match exact condition key: "
                    f"{ip_run} -> {proposed_control}"
                )

            else:
                analysis_status = "ready"
                reason = "unique_exact_condition_input"
                control_run = proposed_control
                pairing_rule = "unique_exact_condition_input"

        else:
            raise ValueError(
                f"Unknown candidate status for {ip_run}: "
                f"{candidate_status}"
            )

        plan.append({
            "analysis_id": ip_run,
            "study_accession": study,
            "ip_run_accession": ip_run,
            "control_run_accession": control_run,
            "declared_target": target,
            "library_layout": layout,
            "instrument_platform": inventory.get(
                "instrument_platform",
                "",
            ),
            "condition_status": condition.get(
                "condition_status",
                "",
            ),
            "condition_rule": condition.get(
                "condition_rule",
                "",
            ),
            "label_context": condition.get(
                "label_context",
                "",
            ),
            "time_value": condition.get(
                "time_value",
                "",
            ),
            "time_unit": condition.get(
                "time_unit",
                "",
            ),
            "candidate_status": candidate_status,
            "candidate_count": str(declared_count),
            "candidate_input_runs": json.dumps(
                candidate_runs,
                separators=(",", ":"),
            ),
            "pairing_rule": pairing_rule,
            "analysis_status": analysis_status,
            "status_reason": reason,
        })

    return plan


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--runs", required=True, type=Path)
    parser.add_argument("--annotations", required=True, type=Path)
    parser.add_argument("--conditions", required=True, type=Path)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--eligibility", required=True, type=Path)
    parser.add_argument("--policy", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)

    args = parser.parse_args()

    plan = build_plan(
        args.runs,
        args.annotations,
        args.conditions,
        args.candidates,
        args.eligibility,
        args.policy,
    )

    fields = [
        "analysis_id",
        "study_accession",
        "ip_run_accession",
        "control_run_accession",
        "declared_target",
        "library_layout",
        "instrument_platform",
        "condition_status",
        "condition_rule",
        "label_context",
        "time_value",
        "time_unit",
        "candidate_status",
        "candidate_count",
        "candidate_input_runs",
        "pairing_rule",
        "analysis_status",
        "status_reason",
    ]

    args.outdir.mkdir(parents=True, exist_ok=True)

    plan_path = args.outdir / "chipseq_analysis_plan.tsv"

    atomic_write(
        plan_path,
        tsv_text(fields, plan),
    )

    status_counts = Counter(
        row["analysis_status"]
        for row in plan
    )

    ready_by_target = Counter(
        row["declared_target"]
        for row in plan
        if row["analysis_status"] == "ready"
    )

    input_paths = {
        "runs": args.runs,
        "annotations": args.annotations,
        "conditions": args.conditions,
        "candidates": args.candidates,
        "eligibility": args.eligibility,
        "policy": args.policy,
    }

    summary = {
        "schema_version": 1,
        "scope": "incremental_chipseq_analysis_planning",
        "ip_analyses": len(plan),
        "status_counts": dict(sorted(status_counts.items())),
        "ready_by_target": dict(sorted(ready_by_target.items())),
        "input_sha256": {
            key: sha256(path)
            for key, path in sorted(input_paths.items())
        },
        "analysis_plan_sha256": sha256(plan_path),
    }

    summary_path = args.outdir / "analysis_plan_summary.json"

    atomic_write(
        summary_path,
        json.dumps(
            summary,
            indent=2,
            sort_keys=True,
        ) + "\n",
    )

    print(
        "[OK] Incremental ChIP-seq analysis plan built: "
        f"{len(plan)} IP analyses"
    )

    for status, count in sorted(status_counts.items()):
        print(f"{status}: {count}")


if __name__ == "__main__":
    main()
