# DotMatch Open-Core Boundary

DotMatch Core is the public trust layer:

- exact known-target assignment engine;
- `k=1` Levenshtein and Hamming modes;
- unique / ambiguous / none / invalid semantics;
- target audit;
- FASTQ/FASTQ.gz count CLI;
- basic QC and HTML reports;
- public schemas;
- public benchmarks and validation harnesses.

These pieces should remain inspectable, installable, benchmarkable, and citable. The scientific claim depends on users being able to reproduce the assignment semantics.

DotMatch Pro / Workbench is the commercial operations layer:

- hosted or private run workbench;
- team workspaces and permissions;
- run history, provenance, and audit logs;
- LIMS / ELN / Benchling-style integrations;
- S3/GCS/Azure/HPC/Nextflow orchestration;
- enterprise report publishing;
- advanced root-cause diagnosis;
- quality-aware confidence and misassignment-risk models;
- support, SLAs, and validated deployment bundles.

The boundary is file-contract based. DotMatch Core emits stable TSV/JSON/HTML artifacts; a proprietary workbench consumes those artifacts and adds operational value without forking the scientific engine.

Recommended core license: Apache-2.0. It is permissive, enterprise-friendly, and includes an explicit patent grant for the open core.

Do not hide:

- assignment status semantics;
- ambiguity policy;
- edit provenance;
- target audit outputs;
- validation/oracle behavior;
- benchmark commands and raw CSVs.

Keep private initially:

- quality-aware confidence scoring;
- run-level root-cause diagnosis heuristics;
- private benchmark corpus and instrument signatures;
- hosted workbench;
- customer integrations;
- workflow monitoring and support automation.
