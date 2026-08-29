# DotMatch evaluation request

DotMatch is looking for a small number of independent evaluations on real
known-target short-DNA workflows. This is a request for evidence, not a claim
of adoption or performance.

## Who this is for

- CRISPR guide counting or screen QC;
- barcode demultiplexing or whitelist assignment;
- primer, amplicon, or targeted-panel matching;
- a core facility or production pipeline that already has a deterministic
  target set and a measurable counting or QC step.

The workflow must have a declared target set and an explicit matching contract.
DotMatch is not a genome aligner and does not emit SAM/BAM or CIGAR output.

## Smallest useful evaluation

1. Install a tagged DotMatch release from the documented package or source
   route.
2. Record the target-file format, read window, metric (`hamming` or
   `levenshtein`), threshold, ambiguity policy, and DotMatch version.
3. Run the equivalent existing tool or workflow on the same inputs when one
   exists.
4. Compare at least one outcome and one operational measure: assignment
   agreement, ambiguous/none rates, counts, wall time, peak memory, or failure
   rate.
5. Keep private data private. A public-safe report can describe the workflow
   class and aggregate results without exposing sequences or identifiers.

The project can provide a command template and help interpret a result, but an
independent evaluator should own the input, comparator, and conclusion.

## Evidence record

For a public-safe record, open an issue or pull request containing:

```text
DotMatch version and install route:
Workflow/data class:
Input scope and read window:
Metric, threshold, and ambiguity policy:
Comparator or truth source:
Reads/records evaluated:
Expected result:
Observed result:
Operational metric (if measured):
Hardware/software environment:
What this does not establish:
Publication/permission status:
```

Do not include private sequences, patient information, proprietary sample
identifiers, or an organization's name without permission. A project-authored
example, an open integration submission, a download count, or a repository
view is not an independent impact record.

## How records are counted for the 10x goal

Scientific records require an independently authored public URL, release and
install route, workflow/data scope, expected-versus-observed result, and a
clear boundary on the conclusion. An integration counts only after an upstream
release or public external-use evidence.

Industry records additionally require permission to describe the evaluation
and a measured before/after operational metric. Confidential data may remain
private; the existence and wording of the record still require approval.

To start, use the public repository issue tracker and title the report
`Independent evaluation: <workflow class>`. No data upload is required.
