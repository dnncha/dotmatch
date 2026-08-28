# External Review Packet

This packet is for maintainers, PIs, core-facility leads, and technical
reviewers who need a concise answer to: what is DotMatch today, what can be
reviewed, and what should not be claimed yet?

## Review Summary

DotMatch is a local command-line and Python package for deterministic
known-target assignment of short read windows. It records each read as
`unique`, `ambiguous`, `none`, or `invalid`, and writes workflow-readable TSV,
JSON, FASTQ, and HTML artifacts.

The package is strongest today where the assay has:

- a known target table;
- a fixed, scaffolded, or inferable read window;
- a need to preserve ambiguity, unmatched reads, invalid windows, and unsafe
  correction states;
- downstream tools that can consume count tables, split FASTQs, QC tables, or
  reports.

It should not be reviewed as a genome aligner, basecaller, adapter trimmer,
variant caller, UMI/cell quantifier, or screen-level statistics package.

## What To Check First

| Review area | Start here | Notes |
| --- | --- | --- |
| Installation | `docs/bioinformatics-evaluation.md` | PyPI and Bioconda install checks are verified after public channel propagation; container runtime proof needs an OCI runtime. |
| Validated scope | `docs/scientific-claims.md` | Public claims must stay inside checked lanes and gates. |
| Output contracts | `docs/schemas.md` | TSV/JSON/HTML artifacts are intentionally plain for workflow systems. |
| Workflow handoff | `docs/workflow-submissions.md` | Links exact open nf-core, MultiQC, Galaxy, and Snakemake submissions and their remaining gates. |
| Workflow status | `docs/workflow-adoption.json` | External workflow integration is listed only after acceptance outside this repository. |
| Ecosystem status | `docs/ecosystem-status.md` | Separates local, submitted, accepted, released, and installable states. |
| Registry status | `docs/registries/biotools.yml` | Draft metadata only; not an accepted bio.tools record. |
| Citation | `docs/methods-and-citation.md` | Use release-specific citation text and generated methods artifacts. |

## Minimum Review Commands

```bash
python3 -m pip install dotmatch==0.2.2
dotmatch --version
dotmatch dist ACGT AGGT
```

For source review:

```bash
make test
make cli-test
make python-test
make workflow-examples-ready
make reviewer-readiness-ready
make repository-ready
```

For a claim-sensitive review, run the assay-specific evidence gate named beside
the claim. A generic test pass is not enough to support a new scientific or
comparative statement.

## Reviewer Questions

Ask these before recommending DotMatch for a workflow:

1. Does the assay have known short targets?
2. Is the read window fixed, scaffolded, or inferable?
3. Are ambiguous, unmatched, or invalid reads important to interpretation?
4. Are TSV, JSON, FASTQ, HTML, and methods/citation artifacts sufficient for the
   receiving workflow?
5. Does the claim being made have a matching raw artifact, report, and gate?
6. Is any public adoption record accepted outside this repository?

## Known Blockers To Keep Visible

- `docs/workflow-adoption.json` is still `not_ready`; local examples are not
  accepted external workflow integrations.
- The nf-core, Galaxy IUC, Snakemake wrappers, and MultiQC submissions are open;
  none is an accepted workflow integration unless it is reviewed and merged.
- The Spack package recipe is merged in `spack-packages` `develop`; no dated
  Spack release containing it is claimed here.
- The MultiQC plugin discovery fix is source-only until a later DotMatch
  release includes it.
- The bio.tools metadata is prepared locally, but no public registry record
  exists yet.
- A local BioContainers runtime check still needs a Docker or OCI host; the
  public tag and manifest digest are available.
- Broader BCL/CBCL, Cell Ranger-style quantification, adapter trimming,
  variant calling, and screen-level statistical claims remain outside the
  current evidence boundary.

## Public Wording Rule

Use the narrowest accurate wording. If a sentence would influence purchase,
publication, workflow replacement, or regulated use, it needs a checked artifact
and a scoped wording review.
