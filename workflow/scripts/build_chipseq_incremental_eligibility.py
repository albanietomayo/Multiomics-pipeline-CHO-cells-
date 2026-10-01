#!/usr/bin/env python3
"""Build fail-closed incremental ChIP-seq protocol eligibility decisions.

Previously reviewed runs retain their curated decisions.

Previously unseen runs may be retained automatically only when they fall
inside a protocol envelope already represented by retained baseline runs:
same reviewed study, library role, histone target and resolved experimental
condition.

New studies, conditions, targets, ambiguous labels or documented issues
are routed to review_required rather than being silently approved.
"""

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


SIGNATURE_FIELDS = (
    "study_accession",
    "library_role",
    "declared_target",
    "condition_rule",
    "label_context",
    "time_value",
    "time_unit",
)


def sha256(path):
    path = Path(path)
    h = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


def load_json(path):
    return json.loads(
        Path(path).read_text(encoding="utf-8")
    )


def read_tsv(path):
    with Path(path).open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )
        rows = list(reader)
        fields = reader.fieldnames or []

    if not fields:
        raise ValueError(
            f"Missing TSV header: {path}"
        )

    return rows, fields


def index_rows(rows, key, label):
    result = {}

    for row in rows:
        value = str(
            row.get(key, "")
        ).strip()

        if not value:
            raise ValueError(
                f"Missing {key} in {label}"
            )

        if value in result:
            raise ValueError(
                f"Duplicate {key} in {label}: "
                f"{value}"
            )

        result[value] = row

    return result


def is_true(value):
    return str(value).strip().lower() in {
        "1",
        "true",
        "yes",
    }


def signature(row):
    return tuple(
        str(row.get(field, "")).strip()
        for field in SIGNATURE_FIELDS
    )


def condition_identity(row):
    return {
        key: str(row.get(key, "")).strip()
        for key in (
            "study_accession",
            "library_role",
            "declared_target",
            "condition_status",
            "condition_rule",
            "label_context",
            "time_value",
            "time_unit",
            "documented_issue",
        )
    }


def build(
    runs_path,
    annotations_path,
    conditions_path,
    baseline_review_path,
    baseline_conditions_path,
    study_evidence_path,
    policy_path,
):
    policy = load_json(policy_path)

    if policy.get("schema_version") != 1:
        raise ValueError(
            "Unsupported incremental eligibility policy schema"
        )

    if (
        policy.get("scope")
        != "incremental_experimental_protocol_eligibility"
    ):
        raise ValueError(
            "Unexpected incremental eligibility policy scope"
        )

    baseline = load_json(
        baseline_review_path
    )

    if baseline.get("schema_version") != 1:
        raise ValueError(
            "Unsupported baseline eligibility schema"
        )

    if (
        baseline.get("scope")
        != "experimental_protocol_eligibility_only"
    ):
        raise ValueError(
            "Unexpected baseline eligibility scope"
        )

    study_evidence = load_json(
        study_evidence_path
    )

    if study_evidence.get("schema_version") != 1:
        raise ValueError(
            "Unsupported study evidence schema"
        )

    evidence_studies = study_evidence.get(
        "studies",
        {},
    )

    if not isinstance(
        evidence_studies,
        dict,
    ):
        raise ValueError(
            "study_evidence/studies must be a dictionary"
        )

    automatic_studies = policy.get(
        "automatic_studies",
        {},
    )

    if not isinstance(
        automatic_studies,
        dict,
    ):
        raise ValueError(
            "automatic_studies must be a dictionary"
        )

    runs_rows, _ = read_tsv(
        runs_path
    )

    annotation_rows, _ = read_tsv(
        annotations_path
    )

    condition_rows, _ = read_tsv(
        conditions_path
    )

    baseline_condition_rows, _ = read_tsv(
        baseline_conditions_path
    )

    runs = index_rows(
        runs_rows,
        "run_accession",
        "active run inventory",
    )

    annotations = index_rows(
        annotation_rows,
        "run_accession",
        "active target annotations",
    )

    conditions = index_rows(
        condition_rows,
        "run_accession",
        "active conditions",
    )

    frozen_conditions = index_rows(
        baseline_condition_rows,
        "run_accession",
        "baseline conditions",
    )

    baseline_rows = baseline.get(
        "runs",
    )

    if not isinstance(
        baseline_rows,
        list,
    ):
        raise ValueError(
            "Baseline eligibility has no runs list"
        )

    baseline_decisions = index_rows(
        baseline_rows,
        "run_accession",
        "baseline eligibility",
    )

    active_set = set(runs)

    for label, collection in (
        (
            "active annotations",
            annotations,
        ),
        (
            "active conditions",
            conditions,
        ),
    ):
        if set(collection) != active_set:
            missing = sorted(
                active_set - set(collection)
            )
            extra = sorted(
                set(collection) - active_set
            )

            raise ValueError(
                f"{label} does not match active run inventory; "
                f"missing={missing}, extra={extra}"
            )

    baseline_set = set(
        baseline_decisions
    )

    if not baseline_set.issubset(
        active_set
    ):
        missing = sorted(
            baseline_set - active_set
        )

        raise ValueError(
            "Previously reviewed runs disappeared from "
            f"the active ChIP catalogue: {missing}"
        )

    if set(
        frozen_conditions
    ) != baseline_set:
        raise ValueError(
            "Frozen baseline conditions do not cover "
            "exactly the baseline eligibility cohort"
        )

    reviewed_envelope = set()

    for run in sorted(
        baseline_set
    ):
        decision = baseline_decisions[
            run
        ]

        frozen = frozen_conditions[
            run
        ]

        active_annotation = annotations[
            run
        ]

        active_condition = conditions[
            run
        ]

        for key in (
            "study_accession",
            "library_role",
            "declared_target",
        ):
            expected = str(
                decision.get(
                    key,
                    "",
                )
            ).strip()

            frozen_value = str(
                frozen.get(
                    key,
                    "",
                )
            ).strip()

            active_annotation_value = str(
                active_annotation.get(
                    key,
                    "",
                )
            ).strip()

            active_condition_value = str(
                active_condition.get(
                    key,
                    "",
                )
            ).strip()

            if not (
                expected
                == frozen_value
                == active_annotation_value
                == active_condition_value
            ):
                raise ValueError(
                    "Baseline metadata changed for "
                    f"{run}/{key}"
                )

        if (
            condition_identity(
                active_condition
            )
            != condition_identity(
                frozen
            )
        ):
            raise ValueError(
                "Baseline experimental condition changed: "
                f"{run}"
            )

        if (
            decision.get(
                "experimental_eligibility"
            )
            == "retained_by_rules"
        ):
            reviewed_envelope.add(
                signature(
                    frozen
                )
            )

    output = []

    for run in sorted(
        active_set
    ):
        inventory = runs[
            run
        ]

        annotation = annotations[
            run
        ]

        condition = conditions[
            run
        ]

        study = str(
            annotation.get(
                "study_accession",
                "",
            )
        ).strip()

        role = str(
            annotation.get(
                "library_role",
                "",
            )
        ).strip()

        target = str(
            annotation.get(
                "declared_target",
                "",
            )
        ).strip()

        if (
            str(
                inventory.get(
                    "study_accession",
                    "",
                )
            ).strip()
            != study
        ):
            raise ValueError(
                f"Inventory/annotation study mismatch: {run}"
            )

        if (
            str(
                condition.get(
                    "study_accession",
                    "",
                )
            ).strip()
            != study
        ):
            raise ValueError(
                f"Condition/annotation study mismatch: {run}"
            )

        if (
            str(
                condition.get(
                    "library_role",
                    "",
                )
            ).strip()
            != role
        ):
            raise ValueError(
                f"Condition/annotation role mismatch: {run}"
            )

        if (
            str(
                condition.get(
                    "declared_target",
                    "",
                )
            ).strip()
            != target
        ):
            raise ValueError(
                f"Condition/annotation target mismatch: {run}"
            )

        if run in baseline_decisions:
            original = dict(
                baseline_decisions[
                    run
                ]
            )

            original[
                "decision_origin"
            ] = (
                "baseline_curated_review"
            )

            output.append(
                original
            )
            continue

        state = (
            policy.get(
                "default_new_run_status",
                "review_required",
            )
        )

        if state != "review_required":
            raise ValueError(
                "Fail-closed default must be review_required"
            )

        decision = {
            "run_accession": run,
            "study_accession": study,
            "library_role": role,
            "declared_target": target,
            "experimental_eligibility": "review_required",
            "rule_id": "incremental_fail_closed_v1",
            "reason": "",
            "evidence_scope": (
                "incremental_metadata_and_reviewed_protocol_envelope"
            ),
            "decision_origin": (
                "incremental_fail_closed_review"
            ),
        }

        label_status = str(
            annotation.get(
                "label_status",
                "",
            )
        ).strip()

        condition_status = str(
            condition.get(
                "condition_status",
                "",
            )
        ).strip()

        if (
            policy.get(
                "require_resolved_labels"
            )
            and label_status
            != "resolved_from_labels"
        ):
            decision[
                "reason"
            ] = "label_annotation_unresolved"

        elif (
            policy.get(
                "require_resolved_conditions"
            )
            and condition_status
            != "resolved_from_labels"
        ):
            decision[
                "reason"
            ] = "experimental_condition_unresolved"

        elif role not in {
            "ip",
            "input",
        }:
            decision[
                "reason"
            ] = "library_role_not_supported"

        elif (
            role == "ip"
            and not target
        ):
            decision[
                "reason"
            ] = "ip_target_unresolved"

        elif (
            role == "input"
            and target
        ):
            decision[
                "reason"
            ] = "input_has_unexpected_target"

        elif (
            policy.get(
                "reject_documented_issues_for_new_runs"
            )
            and is_true(
                condition.get(
                    "documented_issue",
                    "",
                )
            )
        ):
            decision[
                "reason"
            ] = "documented_issue_requires_review"

        elif study not in automatic_studies:
            decision[
                "reason"
            ] = "new_study_requires_protocol_review"

        elif not automatic_studies[
            study
        ].get(
            "automatic_new_runs",
            False,
        ):
            decision[
                "reason"
            ] = (
                "study_not_enabled_for_automatic_incremental_eligibility"
            )

        elif study not in evidence_studies:
            decision[
                "reason"
            ] = "study_evidence_missing"

        elif (
            str(
                evidence_studies[
                    study
                ].get(
                    "evidence_type",
                    "",
                )
            )
            != str(
                automatic_studies[
                    study
                ].get(
                    "required_evidence_type",
                    "",
                )
            )
        ):
            decision[
                "reason"
            ] = "study_evidence_type_changed"

        elif signature(
            condition
        ) not in reviewed_envelope:
            decision[
                "reason"
            ] = "run_outside_reviewed_protocol_envelope"

        else:
            decision[
                "experimental_eligibility"
            ] = "retained_by_rules"

            decision[
                "rule_id"
            ] = (
                "incremental_reviewed_protocol_envelope_v1"
            )

            decision[
                "reason"
            ] = (
                "matches_retained_baseline_study_role_target_condition_envelope"
            )

            decision[
                "evidence_scope"
            ] = (
                "baseline_review_plus_study_evidence_and_run_metadata"
            )

            decision[
                "decision_origin"
            ] = (
                "incremental_reviewed_envelope"
            )

        output.append(
            decision
        )

    return output


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--runs",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--annotations",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--conditions",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--baseline-review",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--baseline-conditions",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--study-evidence",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--policy",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )

    args = parser.parse_args()

    decisions = build(
        args.runs,
        args.annotations,
        args.conditions,
        args.baseline_review,
        args.baseline_conditions,
        args.study_evidence,
        args.policy,
    )

    baseline = load_json(
        args.baseline_review
    )

    baseline_runs = {
        row["run_accession"]
        for row in baseline["runs"]
    }

    eligibility_counts = Counter(
        row[
            "experimental_eligibility"
        ]
        for row in decisions
    )

    origin_counts = Counter(
        row[
            "decision_origin"
        ]
        for row in decisions
    )

    result = {
        "schema_version": 1,
        "scope": "experimental_protocol_eligibility_only",
        "mode": "incremental_baseline_plus_reviewed_protocol_envelope",
        "baseline_run_count": len(
            baseline_runs
        ),
        "active_run_count": len(
            decisions
        ),
        "new_run_count": (
            len(decisions)
            - len(baseline_runs)
        ),
        "eligibility_counts": dict(
            sorted(
                eligibility_counts.items()
            )
        ),
        "decision_origin_counts": dict(
            sorted(
                origin_counts.items()
            )
        ),
        "input_sha256": {
            "runs": sha256(
                args.runs
            ),
            "annotations": sha256(
                args.annotations
            ),
            "conditions": sha256(
                args.conditions
            ),
            "baseline_review": sha256(
                args.baseline_review
            ),
            "baseline_conditions": sha256(
                args.baseline_conditions
            ),
            "study_evidence": sha256(
                args.study_evidence
            ),
            "policy": sha256(
                args.policy
            ),
        },
        "runs": decisions,
    }

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = args.output.with_name(
        args.output.name + ".tmp"
    )

    temporary.write_text(
        json.dumps(
            result,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    temporary.replace(
        args.output
    )

    print(
        "[OK] Incremental ChIP-seq eligibility built"
    )

    print(
        f"active_runs={result['active_run_count']}"
    )

    print(
        f"baseline_runs={result['baseline_run_count']}"
    )

    print(
        f"new_runs={result['new_run_count']}"
    )

    for state, count in sorted(
        eligibility_counts.items()
    ):
        print(
            f"{state}: {count}"
        )

    for origin, count in sorted(
        origin_counts.items()
    ):
        print(
            f"{origin}: {count}"
        )


if __name__ == "__main__":
    main()
