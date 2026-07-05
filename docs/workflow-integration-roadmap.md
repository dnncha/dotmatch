# DotMatch Workflow Integration Roadmap

This roadmap builds on the workflow integration kit in
`docs/workflow-integration-kit.md`. It turns reviewer-facing work into concrete
assets for maintainers, workflow reviewers, core facilities, CRISPR teams, and
assay developers without expanding DotMatch's validated scope.

Use `docs/workflow-integration-plan.json` as the machine-readable checklist for
these ten assets. Reviewer readiness is checked by
`docs/reviewer-readiness.json` and `scripts/check_reviewer_readiness_assets.py`.
Public use records follow the approval requirements in `docs/adopters/README.md`.

## 1. Reviewer Decision Tree

Use this when someone asks whether DotMatch is relevant.

| Question | If yes | If no |
| --- | --- | --- |
| Do you already know the expected short target sequences? | Continue. | DotMatch is probably not the right first tool. |
| Is the read window fixed, scaffolded, or inferable? | Use count, demux, assay, barcode, or panel docs. | Use extraction or preprocessing before DotMatch. |
| Do ambiguous, unmatched, or invalid reads affect interpretation? | Lead with reliability and output artifacts. | Use DotMatch only if ordinary count outputs still help. |
| Do you need genome coordinates, CIGAR strings, UMI deduplication, or hit statistics? | Use established downstream tools. | DotMatch can own the assignment layer. |
| Do you need a public pipeline integration? | Start with workflow submissions. | Start with local install and methods text. |

## 2. Persona One-Pagers

### Core Facilities

- Promise: visible assignment QC for known guides, barcodes, and panels.
- Proof path: homepage, barcode troubleshooting, evidence gallery, packaging.
- First command: `dotmatch barcode autopsy ...`
- Review artifact: HTML report plus TSV/JSON outputs for lab handoff.

### CRISPR Screen Teams

- Promise: guide counts with explicit ambiguity and unmatched diagnostics.
- Proof path: CRISPR tutorial, public CRISPR evidence, methods citation text.
- First command: `dotmatch crispr-count ...`
- Review artifact: count matrix, sample QC, top unmatched, MAGeCK-compatible table.

### Workflow Maintainers

- Promise: stable TSV, JSON, FASTQ, and HTML artifacts for wrappers.
- Proof path: workflow submission pack, schemas, MultiQC parser, release gates.
- First command: `make workflow-examples-ready`
- Review artifact: wrapper fixture, command log, expected output contract.

### Assay Developers

- Promise: barcode panel design and correction-safety review before sequencing.
- Proof path: barcode panel design docs, assay evidence, panel report.
- First command: `dotmatch panel design ...`
- Review artifact: design report, collision tables, plate layout, lab README.

## 3. Integration Target Tracker

Track external integration work separately from local examples.

| Target | Why it matters | Ready asset | Public state to record |
| --- | --- | --- | --- |
| nf-core modules | High-trust workflow reuse and container automation | `docs/workflow-submissions.md` | Merged PR or released module page |
| MultiQC module | Makes DotMatch visible in existing pipeline reports | `python/dotmatch/multiqc.py` | Released plugin or upstream integration |
| Galaxy/IUC | Reaches core facilities and wet-lab teams | Galaxy wrapper examples | IUC acceptance or ToolShed publication |
| Snakemake wrapper | Easy lab workflow reuse | Snakemake example workflow | Public wrapper or external lab pipeline |
| bio.tools entry | Searchable bioinformatics registry presence | Homepage and metadata | Accepted bio.tools record |

## 4. Reviewer Packet

Send this packet when a maintainer, reviewer, procurement evaluator, or PI asks
what is real today.

- Positioning: homepage.
- Bioinformatics evaluation: `docs/bioinformatics-evaluation.md`.
- External review packet: `docs/external-review-packet.md`.
- Install proof: PyPI, Bioconda, packaging notes.
- Output contract: schemas and command reference.
- Validated scope: scientific claims and trust/scope docs.
- Evidence: evidence gallery and benchmark pages.
- Citation: methods and citation template.
- Public use records: adopter notes and workflow status JSON.
- Integration status: `docs/integration-targets.json`.

## 5. Conference Abstracts

### Short Abstract

DotMatch is a deterministic known-target sequencing assignment toolkit for
fixed read windows such as CRISPR guides, inline barcodes, feature tags,
primers, and panel targets. It reports unique, ambiguous, unmatched, and invalid
read outcomes so assignment failures remain visible in workflow artifacts.

### Methods Abstract

Known-target sequencing workflows often collapse assignment decisions into count
tables, making ambiguous reads, shifted windows, unsafe correction, and recurring
unmatched sequences hard to inspect. DotMatch separates the assignment layer
from downstream interpretation: it compares configured read windows with known
short DNA targets, records explicit read outcomes, and writes TSV, JSON, FASTQ,
and HTML artifacts for workflow review. Public claims are scoped to checked
repository evidence and release gates.

### Core Facility Abstract

DotMatch helps sequencing cores review known-target assays before results leave
the facility. It supports guide counting, inline barcode demultiplexing, barcode
panel design, and assignment autopsy reports while preserving ambiguity,
unmatched reads, and invalid extraction windows as visible QC signals.

## 6. Technical Communication Pack

### Technical Thread

```text
DotMatch focuses on one layer: assigning fixed read windows to known short DNA
targets.

Why that matters:
1. unique reads can be counted
2. ambiguous reads stay out of forced calls
3. unmatched reads remain reviewable
4. invalid windows are visible QC failures

Docs: https://dotmatch.readthedocs.io/
Homepage: https://dnncha.github.io/dotmatch
```

### Forum Prompt

```text
I am looking for feedback from teams that run known-target sequencing assays:
CRISPR guide counting, inline barcode demultiplexing, feature tags, primers, or
panel starts. DotMatch is scoped to assignment reliability, not downstream
screen statistics or genome alignment. Which workflow wrapper would make review
easiest for your lab: nf-core, MultiQC, Galaxy, Snakemake, or something else?
```

### Release Follow-Up

```text
The latest DotMatch release gates package installability, docs, scientific
validated scope, workflow examples, and public evidence checks before release
tagging. The project is looking for reviewed workflow integrations and scoped
technical feedback.
```

## 7. Maintainer Issue Templates

### nf-core / Workflow Module Opening

```text
Title: Add DotMatch known-target assignment module

DotMatch assigns fixed read windows to known short DNA targets and writes TSV,
JSON, FASTQ, and HTML outputs. This module proposal is scoped to assignment
artifacts and explicit unique/ambiguous/unmatched/invalid outcomes.

Review asks:
- command shape and metadata
- output contract
- container pinning
- fixture coverage
- MultiQC compatibility
```

### MultiQC Opening

```text
Title: Parse DotMatch assignment QC outputs

DotMatch writes sample QC, summaries, top-unmatched rows, and panel-safety
outputs for known-target sequencing assignments. A MultiQC module should expose
assignment rate, ambiguity rate, unmatched rate, invalid windows, and panel
safety status without implying downstream biological pass/fail calls.
```

## 8. Evaluation Scorecard

Use `docs/pilot-program.md` for intake fields, review steps, and public-use
record requirements.

| Dimension | Score | Notes |
| --- | --- | --- |
| Install worked from released channel | 0-2 | PyPI, Bioconda, source, or container |
| Input mapping was understandable | 0-2 | Target table, sample table, read window |
| Assignment failures were clearer | 0-2 | Ambiguous, unmatched, invalid, unsafe correction |
| Outputs fit existing workflow | 0-2 | TSV, JSON, FASTQ, HTML, MultiQC, notebook |
| Citation and methods text was usable | 0-2 | Version, command, ambiguity policy |
| Public use record approved | 0-2 | Approved wording and public URL, if applicable |

Interpretation:

- 0-4: workflow fit is poor or incomplete.
- 5-8: useful evaluation; record blockers and repeat after fixes.
- 9-12: strong candidate for a documented public workflow example or approved
  use record.

## 9. Integration Tracking Metrics

Track integration status separately from scientific evidence.

| KPI | Source | Cadence |
| --- | --- | --- |
| Homepage visits to install clicks | site analytics if enabled | monthly |
| Docs visits to tutorial starts | docs analytics if enabled | monthly |
| External workflow PRs opened | GitHub URLs | weekly during push |
| External workflow PRs merged | accepted public records | release cycle |
| Public use records approved | `docs/adopters/` | release cycle |
| Citation artifacts generated | release or assay outputs | release cycle |
| Distribution channel health | `make distribution-channels` | release cycle |

Do not combine these metrics with performance or correctness claims. Workflow
interest can increase before external integration is accepted.

## 10. Release Communications Calendar

Use this around each release or major integration push.

| Time | Action | Evidence link |
| --- | --- | --- |
| T-7 days | Confirm validated scope and release notes | `docs/scientific-claims.md` |
| T-5 days | Prepare maintainer issue or PR drafts | `docs/workflow-submissions.md` |
| T-3 days | Prepare short social and forum posts | this page |
| Tag day | Announce only after release workflow artifacts are visible | release URL |
| T+1 day | Verify PyPI, Bioconda, containers, Zenodo as applicable | `docs/distribution-release.json` |
| T+7 days | Follow up with maintainers and pilot contacts | public URLs only |
| T+30 days | Update integration metrics snapshot | public records only |

## Completion Rule

These assets are complete only when the homepage, docs index, integration kit,
structured plan, and site guard all reference the same ten items. If the JSON
plan and markdown playbook diverge, treat the playbook as not release-ready.
