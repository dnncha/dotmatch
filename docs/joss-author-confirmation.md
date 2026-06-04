# JOSS Author Confirmation

Confirm these details immediately before submitting DotMatch to JOSS. They are
not technical blockers, but they are account-side facts that should come from
the submitting author rather than from repository automation.

## Author Metadata

- Author name: `Donncha O'Toole`
- Corresponding author: yes
- Affiliation: `Independent researcher, Ireland`
- ORCID: `0009-0003-5012-7229`
- Email/contact: use the email associated with the JOSS submitter account.

If the affiliation, ORCID, or author list changes, update `paper/paper.md` and
rerun:

```bash
make joss-paper-ready
```

## Submission Choices

- Target journal: Journal of Open Source Software.
- Repository: <https://github.com/dnncha/dotmatch>.
- Paper path: `paper/paper.md`.
- Software version for the current submission packet: `v0.1.7`.
- Software archive DOI: `10.5281/zenodo.20541629`.
- License: Apache-2.0.
- Public tracker for preparation: <https://github.com/dnncha/dotmatch/issues/43>.

## Author Statements

Review and confirm that these statements are accurate before submission:

- The paper's AI usage disclosure correctly describes Codex assistance.
- The paper's acknowledgement section does not omit funding, institutional, or
  contributor credit that should be declared.
- The software scope is correct: DotMatch is a deterministic known-target
  short-DNA assignment tool, not a genome aligner, adapter trimmer, BCL
  converter, UMI/cell quantification pipeline, variant caller, or downstream
  CRISPR statistics package.
- The current Zenodo DOI archives the release being submitted.

## Submission-Day Checks

Run these checks on the final branch state used for submission:

```bash
make joss-paper-ready
make citation-metadata-ready
make docs-ready
```

Then confirm that the latest GitHub Actions runs for `ci`, `codeql`, `pages`,
and `JOSS draft PDF` are green on `main`.
