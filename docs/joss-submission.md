# JOSS Submission Notes

This page tracks the short software-paper path for DotMatch. The draft lives in
`paper/paper.md` with references in `paper/paper.bib`.

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
- Run a JOSS PDF build with the Open Journals/Inara tool or the JOSS GitHub
  Action.
- Confirm `make joss-paper-ready`, `make citation-metadata-ready`, and
  `python3 -m pytest python/tests` pass after any paper-related metadata edits.

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
