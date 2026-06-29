# DotMatch Industry Exposure Kit

This page turns the adoption roadmap into concrete marketing and distribution
work. It is about reaching the right people without overstating claims. Public
performance, correctness, packaging, and adoption statements still need the
evidence gates in `docs/scientific-claims.md`, `docs/packaging.md`, and
`docs/workflow-adoption.json`.

## The Big 5 Wins

| Win | Relevant people | Asset in this repository | Done when |
| --- | --- | --- | --- |
| Audience-specific homepage routes | Core facility leads, CRISPR screen teams, assay developers, workflow maintainers | Public homepage section: "Routes into industry workflows" | Each audience has a clear use case, proof path, and next-click destination. |
| Workflow distribution handoff | nf-core, MultiQC, Galaxy/IUC, Snakemake, institutional pipeline owners | `docs/workflow-submissions.md` | At least one external PR, package, or reviewed wrapper is public and tracked in `docs/workflow-adoption.json`. |
| Citation and methods flywheel | PIs, methods writers, bioinformatics leads, paper authors | `docs/methods-and-citation.md`, `docs/citation-flywheel.md` | Runs can produce copyable methods and citation artifacts, and external users know how to cite the exact release. |
| Evidence-first launch packet | Technical evaluators, skeptical industry users, procurement reviewers | Evidence gallery, benchmark docs, scientific claim boundary | Every public announcement links to scoped evidence instead of broad claims. |
| Public adopter record | Labs, cores, CROs, biotech, pharma teams, workflow projects | `docs/adopters/README.md` | Quote-approved independent use is listed only after the external artifact is public. |

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
- Link `docs/scientific-claims.md` before making evidence-sensitive claims.
- Link `docs/workflow-submissions.md` for pipeline maintainers.
- Link `docs/methods-and-citation.md` for citation and methods text.
- Include PyPI, Bioconda, and repository links only after release smoke tests pass.
- Record any accepted external integration in `docs/workflow-adoption.json`.
- Add public pilots to `docs/adopters/` only with approval and a public URL.

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
and adoption-record rules:
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

If a small public or anonymized pilot is possible, I can help scope it so the
output is useful without exposing private sample data.
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
   pass. This reaches pipeline authors and inherited container/test automation.
2. Package or upstream the MultiQC module so DotMatch outputs are automatically
   visible in existing reports.
3. Ask two to five core facilities or CRISPR-screen teams for scoped pilot
   feedback. Do not count the pilot as adoption until they approve a public
   record.
4. Convert the best pilot into a methods-focused example, not a broad marketing
   claim.
5. Keep the homepage, docs index, citation page, and adopter records aligned so
   every external mention points to the same source of truth.

## Tracking Rules

- Private feedback can inform the roadmap, but it is not adoption evidence.
- Unmerged PRs can be listed as outreach activity, but not as external adoption.
- External use counts only when a public artifact can be linked.
- Quote-approved user names, organization names, and logos need explicit
  approval before they appear in the repository or public site.
- If a claim would influence purchase, publication, or pipeline replacement, it
  needs a checked artifact and a scoped wording review.
