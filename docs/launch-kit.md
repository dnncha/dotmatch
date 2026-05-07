# GitHub Launch Kit

Use this page when turning the local repository into the public GitHub project, release, registry submissions, and launch posts. Keep wording aligned with `docs/scientific-claims.md`.

## Repository Setup

Repository:

```text
Dnncha/dotmatch
```

Description:

```text
Fast exact short-DNA known-target assignment for CRISPR guides, barcodes, primers, panels, and whitelists.
```

Topics:

```text
bioinformatics
computational-biology
crispr
fastq
barcode-demultiplexing
edit-distance
sequence-analysis
genomics
c
python
```

Enable:

- Issues;
- Discussions;
- Dependabot alerts;
- GitHub Pages with GitHub Actions as the source;
- Code scanning with the committed CodeQL workflow;
- private vulnerability reporting;
- branch protection after the first CI pass on `main`.

Set the repository website to the Pages URL after the first `pages` workflow succeeds. For `Dnncha/dotmatch`, that should be:

```text
https://dnncha.github.io/dotmatch
```

Use `public/dotmatch-og.png` as the GitHub repository social preview image.

Do not tag `v0.1.0` until the pre-tag checks in `docs/release-process.md` pass on the pushed repository.

## First Release Notes

Title:

```text
DotMatch v0.1.0: exact known-target short-DNA assignment for real FASTQ workflows
```

Body:

```markdown
DotMatch v0.1.0 is the first public release of the open core: a native C and Python toolkit for exact known-target short-DNA assignment in CRISPR guide-counting, barcode, primer, panel, and whitelist-style workflows.

Highlights:

- deterministic `unique`, `ambiguous`, `none`, and `invalid` assignment semantics;
- `k=0` exact lookup and `k=1` Hamming/Levenshtein assignment lanes;
- native FASTQ/FASTQ.gz counting, demultiplexing, audit, validation, and QC report commands;
- MAGeCK-compatible CRISPR count output;
- Python package builds with bundled native library for local/GitHub artifact testing;
- reproducible benchmark reports, raw CSV evidence, and claim gates;
- Apache-2.0 open-core boundary with public schemas and validation harnesses.

Current supported claim:

DotMatch supports narrow known-target short-DNA assignment and checked CRISPR guide-counting claims where `make public-crispr-claim-gate` and `make crispr-sota-gate` pass.

Blocked claims:

- barcode demultiplexing SOTA remains blocked until `make barcode-sota-gate` passes on claim-grade real-data evidence;
- raw BCL/CBCL demultiplexing SOTA remains blocked until `make bcl-sota-gate` passes;
- DotMatch is not a genome aligner and does not claim SAM/BAM/CIGAR/reference-mapping semantics.

Reproducibility entry points:

- `docs/scientific-claims.md`
- `docs/benchmarks/`
- `docs/methods-and-citation.md`
- `CHANGELOG.md`
```

## Short Public Pitch

```text
DotMatch is an Apache-2.0 open-core toolkit for exact known-target short-DNA assignment. It targets CRISPR guides, barcodes, primers, panels, and whitelist-like sequences where deterministic unique/ambiguous/no-match semantics matter more than broad genome alignment.
```

## Social Post

```text
I am open-sourcing DotMatch: a fast exact known-target short-DNA assignment engine for CRISPR guide-counting, barcode-style workflows, primers, panels, and whitelists.

It is not a genome aligner. The focus is deterministic `unique` / `ambiguous` / `no-match` assignment, real FASTQ workflows, public benchmarks, and claim gates that fail when evidence is missing.

Repo: https://github.com/Dnncha/dotmatch
```

## Scientific Registry Metadata

Name:

```text
DotMatch
```

Description:

```text
DotMatch performs exact known-target short-DNA assignment for FASTQ workflows, including CRISPR guide counting, barcode-style demultiplexing, primer/panel matching, and whitelist-like sequence assignment.
```

Homepage and source:

```text
https://github.com/Dnncha/dotmatch
```

License:

```text
Apache-2.0
```

Keywords:

```text
bioinformatics; computational biology; CRISPR; FASTQ; barcode demultiplexing; edit distance; known-target assignment; sequence analysis
```

Languages:

```text
C; Python; TypeScript
```

Primary functions:

```text
short sequence assignment; guide counting; barcode demultiplexing; edit-distance matching; FASTQ quality control; target-library audit
```

## Citation And Archive

Before DOI:

```text
O'Toole D. DotMatch: Streaming Exact One-Edit Barcode and Guide Assignment Without Exhaustive Scanning. Software release v0.1.0. https://github.com/Dnncha/dotmatch
```

After Zenodo creates the DOI:

- add the DOI to `CITATION.cff`;
- add the DOI badge to `README.md`;
- update `docs/methods-and-citation.md`;
- rerun `make publication-ready`.

## Open-Core Messaging

Say:

```text
DotMatch Core is the public scientific trust layer: assignment engine, CLI, Python bindings, audit, schemas, validation harnesses, and benchmarks.
```

Say:

```text
Future proprietary work can add hosted workbench, team workflows, LIMS/ELN integrations, run history, support, deployment bundles, and operational automation without hiding the core assignment semantics.
```

Do not say:

```text
The open repository is a teaser, demo, or limited proof of concept.
```

The public core must remain useful and citable on its own.

## Do-Not-Claim List

Avoid these phrases until the relevant gates and artifacts exist:

- state-of-the-art barcode demultiplexing;
- production Illumina demultiplexing replacement;
- universal guide-counter replacement;
- general aligner;
- Edlib replacement;
- genome-scale mapper;
- supports SAM/BAM/CIGAR output.

Preferred phrase:

```text
exact known-target short-DNA assignment with deterministic ambiguity handling and reproducible claim gates
```
