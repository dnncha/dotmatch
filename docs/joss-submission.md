# JOSS Submission Notes

This page tracks the short software-paper path for DotMatch. The draft lives in
`paper/paper.md` with references in `paper/paper.bib`.

Public tracking issue: <https://github.com/dnncha/dotmatch/issues/43>.

## Target

Primary target: Journal of Open Source Software (JOSS).

Rationale:

- JOSS papers are short software papers designed to give research software a
  normal scholarly citation.
- DotMatch already has a public repository, OSI-approved license, tests,
  release metadata, and a Zenodo DOI.
- A JOSS paper is a better near-term fit than a bioinformatics application note
  while broader package-channel and workflow adoption evidence is still being
  completed.

## Before Submission

- Confirm the author affiliation and optional ORCID in `paper/paper.md`.
- Decide whether to keep the current title or use the longer release title.
- Review the AI usage disclosure for accuracy.
- Make sure the public README points to the paper once submitted or accepted.
- Keep the Zenodo DOI badge and BibTeX block in the README aligned with the
  current archived release.
- Run the `JOSS draft PDF` GitHub Actions workflow and download its
  `joss-paper` artifact.
- Confirm `make joss-paper-ready`, `make citation-metadata-ready`, and
  `python3 -m pytest python/tests` pass after any paper-related metadata edits.

## Submission Packet

Use the JOSS submission site: <https://joss.theoj.org/>.

Suggested form values:

| Field | Value |
| --- | --- |
| Journal | Journal of Open Source Software |
| Software repository | <https://github.com/dnncha/dotmatch> |
| Repository branch | `main` |
| Paper path | `paper/paper.md` |
| Title | DotMatch: deterministic known-target short-DNA assignment for sequencing workflows |
| Software version | `v0.1.7` for the current archived software release |
| Software archive DOI | `10.5281/zenodo.20541629` |
| Concept DOI | `10.5281/zenodo.20541628` |
| License | Apache-2.0, plain-text `LICENSE` file |
| Primary language/runtime | C command-line tool with Python package wrappers |
| Issue tracker | <https://github.com/dnncha/dotmatch/issues> |

Suggested short description:

> DotMatch is an open-source command-line and Python package for deterministic
> assignment of fixed short-DNA read windows to known target tables, including
> CRISPR guide counting, barcode-style demultiplexing, feature-barcode
> assignment, and panel-design safety checks. It reports unique, ambiguous,
> unmatched, and invalid reads explicitly and includes public-data evidence
> gates for its documented claims.

Important submission note:

The current Zenodo DOI archives the `v0.1.7` software release. If reviewers ask
for paper, documentation, or software changes, make a new tagged release after
the accepted review state and archive that final release on Zenodo before
posting the final version and DOI in the JOSS review issue.

## Review-Readiness Map

JOSS review checklist item | DotMatch evidence
--- | ---
Repository is public and cloneable | GitHub repository: <https://github.com/dnncha/dotmatch>
Issue tracker is public | GitHub issues are enabled: <https://github.com/dnncha/dotmatch/issues>
OSI-approved license | `LICENSE` contains Apache License 2.0 text
Author contribution | Commit history under `dnncha/dotmatch`; confirm final author list before submission
Research application | `paper/paper.md`, `README.md`, `docs/scientific-claims.md`, and public benchmark reports
Installation instructions | `README.md#installation`, `docs/getting-started.md`, `docs/packaging.md`
Example usage | `README.md`, `examples/crispr_guides/`, `examples/workflows/`, and evidence-gallery scenarios
Functionality documentation | `README.md`, `docs/assayspec.md`, command help, and workflow docs
Automated tests | `make test`, `make cli-test`, `make python-test`, `make python-package-test`, GitHub Actions `ci`
Claim verification | `docs/scientific-claims.md`, `make assay-evidence-ready`, `make release-ready`
Paper required sections | `make joss-paper-ready` checks metadata, sections, bibliography, DOI, and word count
Paper PDF compilation | GitHub Actions `JOSS draft PDF` workflow builds `paper/paper.pdf`
AI usage disclosure | `paper/paper.md#ai-usage-disclosure`

## Current Automated Evidence

The current `main` branch has passing GitHub Actions for:

- `ci`;
- `codeql`;
- `pages`;
- `JOSS draft PDF`.

Before submitting, rerun the draft PDF workflow and download the `joss-paper`
artifact. The expected artifact is `paper.pdf`.

## Evidence Boundaries For The Paper

The paper should not claim:

- genome alignment;
- production Illumina BCL conversion;
- adapter trimming;
- UMI/cell quantification;
- amplicon consensus generation or variant calling;
- downstream CRISPR phenotype statistics;
- broad superiority over unrelated aligners or demultiplexers.

Claims should stay tied to the checked evidence in `docs/scientific-claims.md`
and public benchmark reports under `docs/benchmarks/`.
