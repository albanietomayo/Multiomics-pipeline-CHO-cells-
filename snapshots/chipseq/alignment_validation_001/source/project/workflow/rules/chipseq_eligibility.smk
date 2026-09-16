# Shared by the metadata and preprocessing entry points.
rule chipseq_experimental_eligibility:
    input:
        decisions="config/chipseq_experimental_eligibility.json",
        samples=CHIP_ELIGIBILITY_SAMPLES,
        conditions=CHIP_ELIGIBILITY_CONDITIONS,
        pilot=([CHIP_ELIGIBILITY_PILOT] if CHIP_ELIGIBILITY_PILOT else []),
        evidence=[
            "snapshots/chipseq/control_validation_001/chipseq_conditions.tsv",
            "snapshots/chipseq/control_validation_001/control_candidate_summary.json",
            "snapshots/chipseq/control_validation_001/request_plan_provenance.json",
        ],
        code="workflow/scripts/validate_chipseq_eligibility.py"
    output:
        CHIP_ELIGIBILITY_REPORT
    params:
        pilot_arg=("--pilot config/chipseq_pilot.json" if CHIP_ELIGIBILITY_PILOT else "")
    shell:
        "python3 {input.code:q} --decisions {input.decisions:q} "
        "--samples {input.samples:q} --conditions {input.conditions:q} "
        "{params.pilot_arg} --report {output:q}"
