# DotMatch Evaluation Protocol

This protocol is for core facilities, CRISPR screen teams, assay-development
groups, CROs, biotech teams, and workflow maintainers assessing DotMatch on a
known-target sequencing workflow.

The goal is to decide whether DotMatch fits an existing assay or pipeline by
checking installation, input mapping, assignment outputs, and methods artifacts
on data the reviewing team is allowed to share or summarize.

## Suitable Workflows

DotMatch is appropriate to evaluate when the workflow has:

- FASTQ or FASTQ.gz reads, or a public/sanitized fixture;
- a known target table such as guides, barcodes, feature tags, primers, or panel
  targets;
- a fixed, scaffolded, or inferable read window;
- an existing output to compare against, such as counts, split FASTQs, or QC
  tables;
- permission to record commands, package versions, and non-sensitive output
  summaries.

DotMatch is not a substitute for genome alignment, basecalling, adapter
trimming, UMI/cell quantification, variant calling, or screen-level statistics.

## Intake Fields

Capture these details before running a comparison:

| Field | Required note |
| --- | --- |
| Project or assay | Public name, anonymized label, or internal reference |
| Assay context | CRISPR, inline barcode, feature barcode, amplicon, adapter prefix, panel, or other known-target workflow |
| Current workflow | Existing command, wrapper, notebook, or manual process |
| Target table shape | Identifier columns, sequence column, expected target length |
| Read window | Start, length, read mate, and whether the window is inferred |
| Correction policy | Exact, Hamming, Levenshtein, radius, quality filter, or no correction |
| Outputs to inspect | QC report, counts, split FASTQs, per-read assignments, top unmatched, methods text |
| Public-use permission | None, anonymized summary, public project name, or approved record text |

Do not place private FASTQ/BAM/BCL files, patient data, customer assay designs,
or restricted screenshots in the public repository.

## Review Steps

1. Install DotMatch from PyPI or Bioconda unless the review is explicitly testing
   an unreleased checkout.
2. Run the closest tutorial or workflow example first.
3. Audit the target table before enabling correction.
4. Run the evaluation on minimized or public data when possible.
5. Inspect `sample_qc.tsv`, `summary.json`, `assignments.tsv`,
   `top_unmatched.tsv`, and the HTML report before comparing headline counts.
6. Check whether ambiguous reads, unmatched reads, invalid windows, unsafe
   correction, or workflow handoff became easier to review.
7. Record blockers as workflow requirements, not as product claims.

## Scorecard

| Dimension | Score | Notes |
| --- | --- | --- |
| Install from released channel | 0-2 | PyPI, Bioconda, source, or container |
| Input mapping | 0-2 | Target table, sample table, read window |
| Assignment review | 0-2 | Ambiguous, unmatched, invalid, unsafe correction |
| Workflow fit | 0-2 | TSV, JSON, FASTQ, HTML, MultiQC, notebook |
| Methods and citation | 0-2 | Version, command, ambiguity policy |
| Public record permission | 0-2 | Approved wording and URL, if any |

Suggested interpretation:

- `0-4`: workflow fit is poor or incomplete.
- `5-8`: useful evaluation; record blockers and repeat after fixes.
- `9-12`: strong candidate for a documented public workflow example or approved
  use record.

## Public Use Records

Public use records belong in `docs/adopters/` only when the reviewing
organization or project has approved the exact name, URL, scope note, and
evidence link. Until then, keep the evaluation as a private review artifact.
