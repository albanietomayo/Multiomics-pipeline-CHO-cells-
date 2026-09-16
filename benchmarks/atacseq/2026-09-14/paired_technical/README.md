# Historical paired-end technical validation

Base commit: 354cf261d121df32eda46710ecd1383de404bc20, branch atacseq-paired-filter.
The included three working files and patch were supplied from Vera on 2026-09-16.
Job 10280253 completed on 2026-09-14; the earlier two jobs failed.
SRR29929613 is a CHO multiome library and is excluded from the bulk atlas.
These files preserve historical execution details; they are not the current entry point.
The historical worker has a hard-coded excluded run, user paths, and quota checks.
Do not submit it for bulk production. Paired-end technical support does not establish
validation on eligible bulk CHO paired-end libraries. No eligible PE candidate was
identified in the audited catalogue. Raw PE validation reports remain in Vera at:
/cephyr/users/mayoa/Vera/atacseq_paired_validation_results/SRR29929613/job_10280253
