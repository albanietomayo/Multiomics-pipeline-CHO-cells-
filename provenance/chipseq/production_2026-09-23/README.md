# Frozen ChIP-seq production evidence: 2026-09-23

These files document historical production; they are not current execution
configuration. The scientific cohort has 18 PASS analyses, eight narrow and ten
broad, using 18 IPs and four shared Inputs. Historical paths in parameter and
reference records document their original locations and are not required by the
portable workflow. Existing public summaries remain under
`benchmark/derived/chipseq_source_metadata/`.

`analysis_sources.tsv` maps each analysis to its recorded production commit,
source manifest, production manifest and installed package-record hashes.
`source_versions.tsv` maps original paths and commits to exact source hashes.
There are 27 distinct source blobs in the recorded production source sets.
One current source version is not claimed to have produced all 18 analyses.

| Recorded commit | Analyses |
| --- | --- |
| 48d80c961acfd188033db9e1cf2972222fbc48a5 | SRR20770294, SRR20770296, SRR20770297 |
| fb62508d51efa41ab73067ba7ca57966516cd5f9 | ERR868151–ERR868154, SRR20770291–SRR20770293 |
| 290c35dcec4b76c258b2788a4adbc6940654d04f | ERR868155, ERR868156, ERR868171–ERR868176 |

Changes between these source sets concern download resumption, submission,
locking, protected cleanup, ledger updates, documentation and related tests.
The scientific configs, core alignment/filtering/peak/reference helpers and
combined environment specification are unchanged. The portable implementation
uses reviewed scientific helpers and adapts the worker's execution interface.

## Installed environment records

| Record | Analyses |
| --- | --- |
| environment_a.explicit.txt | ERR868151–ERR868156, ERR868171–ERR868175 |
| environment_b.explicit.txt | ERR868176 |
| environment_c.explicit.txt | SRR20770291–SRR20770293 |
| environment_d.explicit.txt | SRR20770294, SRR20770296, SRR20770297 |

The direct Table A9 versions match in all four records. Some transitive package
builds differ. These are original Linux execution records, not cross-platform
locks or claims about a new environment solve.

## Checksums and intended release assets

`frozen_inputs.sha256` protects the imported scientific inputs and environments.
`frozen_outputs.sha256` records existing public benchmark/integration products
before the patch. `historical_artifacts.tsv` lists external production artifacts
and original hashes. Their large files are not included. `frozen_contract.json`
preserves historical parameters, counts, limitations and release identifiers.
`SHA256SUMS` protects this compact provenance package.

`source_blobs.tar.gz` and `run_records.tar.gz` have not been created or imported.
`source_versions.tsv` defines the former's intended SHA256-addressed members;
`release_asset_contents.tsv` defines the latter's intended members, byte sizes
and hashes. `archived_record_location` is an intended future asset location, not
a claim that an asset is already available. No container hash is invented.
Historical worker and submitter source blobs are evidence, not portable launchers.

Raw ENA XML and its archive are prohibited from Git, GitHub release assets and
Zenodo. Sanitized scientific JSON replaces the runtime XML dependency. The
historical archive is retained privately only; its original hash
is `5b01f03bacf45548c004c70c28e46222e6313a0f457a3549f2165f174c0c54cf`.
Its metadata inventory and checksum record are in
`snapshots/chipseq/xml_validation_001/`. That original checksum list includes the
external archive; absence from the checkout is intentional. Do not run that
list as if every member were bundled.

No authentication credentials are part of this package. Historical execution
paths are preserved rather than silently rewritten. Future asset publication,
full artifact verification and any production rerun remain separate actions.
