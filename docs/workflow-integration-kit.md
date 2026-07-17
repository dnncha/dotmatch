# DotMatch Workflow Integration Kit

This page collects the materials needed to review DotMatch for workflow
integration in nf-core, MultiQC, Galaxy/IUC, Snakemake, or an institutional
pipeline. Public performance, correctness, packaging, and workflow status
statements should remain tied to `docs/scientific-claims.md`,
`docs/packaging.md`, and `docs/workflow-adoption.json`.

## Priority Integration Work

| Work item | Relevant reviewers | Asset in this repository | Done when |
| --- | --- | --- | --- |
| Homepage routes | Core facility leads, CRISPR screen teams, assay developers, workflow maintainers | Public homepage section: "Routes into industry workflows" | Each reviewer has a clear use case, proof path, and next-click destination. |
| Workflow distribution handoff | nf-core, MultiQC, Galaxy/IUC, Snakemake, institutional pipeline owners | `docs/workflow-submissions.md` | At least one external PR, package, or reviewed wrapper is public and tracked in `docs/workflow-adoption.json`. |
| Methods and citation artifacts | PIs, methods writers, bioinformatics leads, paper authors | `docs/methods-and-citation.md`, `docs/citation-flywheel.md` | Runs can produce copyable methods and citation artifacts, and reviewers know how to cite the exact release. |
| Reviewer packet | Technical reviewers, core facilities, procurement reviewers | Evidence gallery, benchmark docs, scientific scope notes | Public statements link to scoped evidence instead of broad claims. |
| Public use records | Labs, cores, CROs, biotech, pharma teams, workflow projects | `docs/adopters/README.md` | Independent use is listed only after the named party approves the wording and public URL. |
| Reviewer readiness gate | Maintainers and reviewers | `docs/reviewer-readiness.json`, `scripts/check_reviewer_readiness_assets.py` | `make reviewer-readiness-ready` passes before reviewer, integration, registry, or public-use materials are merged. |

## Positioning

Short description:

> DotMatch is a deterministic known-target sequencing assignment toolkit. It
> counts or demultiplexes fixed read windows against expected short DNA targets
> while keeping unique, ambiguous, unmatched, and invalid outcomes visible.

One-line audience variants:

- Core facilities: "DotMatch makes barcode and guide assignment failures visible before results leave the core."
- CRISPR screen teams: "DotMatch produces guide-counting artifacts while preserving ambiguous and unmatched reads for review."
- Workflow maintainers: "DotMatch writes stable TSV, JSON, FASTQ, and HTML outputs that can be wrapped in nf-core, Galaxy, Snakemake, and MultiQC."
- Assay developers: "DotMatch designs and audits barcode panels so unsafe correction rules are caught before sequencing."

Avoid these unsupported shortcuts unless the linked evidence explicitly supports
them for the exact release and setting:

- "replacement for genome aligners";
- "replacement for basecallers or full BCL conversion";
- "screen-analysis package";
- "guaranteed production demultiplexing replacement";
- broad speed claims without benchmark scope, hardware, command, and comparator.

## Launch Checklist

Use this checklist when announcing a release, opening an integration PR, or
asking an external maintainer to evaluate DotMatch.

- Link the homepage first for positioning.
- Link `docs/getting-started.md` or a specific tutorial for a runnable path.
- Link `docs/scientific-claims.md` before making scope-sensitive claims.
- Link `docs/workflow-submissions.md` for pipeline maintainers.
- Link `docs/methods-and-citation.md` for citation and methods text.
- Include PyPI, Bioconda, and repository links only after release smoke tests pass.
- Record any accepted external integration in `docs/workflow-adoption.json`.
- Add public use records to `docs/adopters/` only with approved wording and a
  public URL.

## Copy-Paste Outreach

### Repository announcement

```text
DotMatch is a deterministic known-target sequencing assignment toolkit for
CRISPR guides, inline barcodes, feature tags, primers, and panel targets. It
reports every read as unique, ambiguous, unmatched, or invalid, so assignment
failures stay visible in TSV, JSON, FASTQ, and HTML outputs.

Homepage: https://dnncha.github.io/dotmatch
Docs: https://dotmatch.readthedocs.io/
Repository: https://github.com/dnncha/dotmatch
```

### Workflow maintainer email

```text
Subject: DotMatch handoff for known-target sequencing workflows

Hi <name>,

DotMatch may be a fit for workflows that count CRISPR guides, split inline
barcodes, or audit known-target read windows. The useful distinction is that it
keeps unique, ambiguous, unmatched, and invalid read outcomes explicit instead
of collapsing them into a count-only result.

I put together a workflow submission pack with expected outputs, review notes,
and public workflow-status rules:
https://github.com/dnncha/dotmatch/blob/main/docs/workflow-submissions.md

If this overlaps with your pipeline, I would value a review of the output
contract and wrapper shape before opening or expanding an integration PR.
```

### Core facility pilot email

```text
Subject: Pilot request: visible assignment QC for barcodes and guide counts

Hi <name>,

I am looking for feedback from core facilities that run known-target sequencing
assays: CRISPR guide counting, inline barcodes, feature tags, or panel starts.
DotMatch focuses on the assignment reliability layer, especially ambiguous
reads, unsafe correction, shifted windows, and recurring unmatched sequences.

The homepage and evidence boundary are here:
https://dnncha.github.io/dotmatch
https://github.com/dnncha/dotmatch/blob/main/docs/scientific-claims.md

If a small public or anonymized evaluation is possible, I can help scope the
commands and outputs so the review does not require sensitive sample data.
```

### Short social post

```text
DotMatch is a known-target sequencing assignment toolkit for CRISPR guides,
inline barcodes, feature tags, primers, and panel targets. It keeps unique,
ambiguous, unmatched, and invalid read outcomes visible in workflow-friendly
artifacts.

Homepage: https://dnncha.github.io/dotmatch
Docs: https://dotmatch.readthedocs.io/
```

## Where To Put Effort First

1. Open or prepare the official nf-core module PRs after local workflow gates
   pass. This gives pipeline authors concrete module code and test fixtures.
2. Package or upstream the MultiQC module so DotMatch outputs are automatically
   visible in existing reports.
3. Use the [DotMatch Evaluation Protocol](pilot-program.md) with two to five
   core facilities or CRISPR-screen teams.
4. Convert strong evaluations into methods-focused examples only when the
   reviewing team approves the public wording and URL.
5. Keep the homepage, docs index, citation page, and public use records aligned
   so every external mention points to the same source of truth.

## Workflow Integration Roadmap

The next layer is tracked in
[DotMatch Workflow Integration Roadmap](workflow-integration-roadmap.md)
and mirrored in `docs/workflow-integration-plan.json`.

1. Evaluator decision tree.
2. Persona one-pagers.
3. Integration target tracker.
4. Reviewer evidence packet.
5. Conference abstracts.
6. Social and forum pack.
7. Maintainer issue templates.
8. Evaluation scorecard.
9. Integration tracking metrics.
10. Release publication checklist.

## Tracking Rules

- Unmerged PRs can be listed as work in progress, not as accepted integration.
- External use records require an approved public URL.
- Names, organization names, logos, and quotes require approved wording before
  they appear in the repository or public site.
- If a statement would influence purchase, publication, or pipeline replacement,
  it needs a checked artifact and a scoped wording review.
