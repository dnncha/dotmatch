# DotMatch 10x goal: baseline and evidence gates

This document turns the 10x objective into measurable project targets. It is a
planning and evidence boundary, not a claim that any target has already been
met. Accuracy, speed, scientific impact, and industry impact are measured
separately; downloads, stars, open pull requests, and local synthetic tests do
not substitute for the relevant evidence.

## Baseline snapshot

The baseline is the public DotMatch 0.2.2 evidence in this repository. The
canonical working checkout is intentionally not used for this snapshot because
it contains unrelated uncommitted changes.

| Dimension | Current recorded evidence | Baseline used for the goal |
| --- | --- | --- |
| Scientific accuracy | The CRISPR comparison report records six public-data Edlib-oracle rows of 10,000 reads each, all with zero assignment mismatches. Count agreement, ambiguity, and outcome semantics are reported separately. | 60,000 oracle-checked reads, zero observed mismatches. Because zero is an observed floor, not a measurable positive accuracy margin, improvement is defined as 10x broader independent validation with the same zero-mismatch and semantic gates. |
| Speed | The full Sanson/Brunello public Hamming `k=1` row records 634,950.2 reads/s and 388.9288 seconds against guide-counter at 473,609.5 reads/s. The report marks this as a single full-data row, not a repeated claim. | 634,950.2 reads/s and 388.9288 seconds on the recorded workload and hardware. A repeated baseline must be captured before publishing a speed multiplier. |
| Scientific impact | `docs/workflow-adoption.json` is `not_ready` with no accepted external workflow entries. `docs/adopters/` contains only the record template. Open external PRs are not counted as adoption. | Zero verified independent public use records in the project evidence registry. A multiplicative ratio from zero is undefined, so the 10x operational target is ten independently authored, public-safe scientific evaluation or use records. |
| Industry impact | There is no permissioned public industry, core-facility, or production-use record in `docs/adopters/`. Package availability and downloads are distribution evidence, not industry impact. | Zero verified industry-use records in the project evidence registry. The 10x operational target is ten permissioned evaluation or deployment records with workflow scope and at least one before/after operational metric. |

### Baseline boundaries

- The two public CRISPR datasets and existing benchmark lanes demonstrate
  scoped scientific evidence; they do not establish universal accuracy or
  superiority.
- A zero-mismatch result is retained as zero mismatches. It is not converted
  into a claim of perfect accuracy.
- Speed comparisons must preserve input data, target set, read window, metric,
  ambiguity policy, output contract, hardware, software versions, warm-up
  policy, repeat count, and memory measurement.
- Scientific and industry records require an external author or an explicitly
  approved organization. A project-authored integration PR, download counter,
  repository view, or unapproved logo is not an impact record.

## 10x targets

### 1. Scientific accuracy

Pass the accuracy target only when all of the following are true:

1. At least 600,000 reads are checked across at least ten public or
   permission-cleared datasets/lanes that were not used to tune the matcher.
2. Each lane has a declared truth source or independent edit-distance oracle,
   recorded commands, versions, input window, target set, and ambiguity policy.
3. The aggregate assignment mismatch count is zero, and the report retains
   `unique`, `ambiguous`, `none`, and `invalid` outcomes rather than collapsing
   them into a single accuracy number.
4. Count-level agreement and relevant per-read invariants pass for every lane.

With zero mismatches, the simple one-sided Rule-of-3 upper bound moves from
approximately `3/60,000` to `3/600,000`, a tenfold tighter observed-error bound.
This is the defensible interpretation of “10x more accurate” for a system whose
current checked error count is already zero.

### 2. Performance

Pass the speed target only when a five-repeat baseline and five-repeat candidate
are run on the same full Sanson/Brunello Hamming `k=1` workload and the same
declared hardware/software environment:

- candidate throughput: at least 6,349,502 reads/s, or
- candidate wall time: at most 38.8929 seconds,

while preserving zero oracle mismatches, count semantics, output equivalence,
and a separately reported peak-memory budget. If the workload or comparator
changes, it is a new benchmark lane and cannot be called a 10x improvement of
this baseline.

### 3. Scientific impact

Reach ten independently authored public-safe scientific evaluation or use
records. Each record must include a public URL, DotMatch version, install
route, workflow/data class, expected-versus-observed result, and the boundary
of what was actually tested. Accepted nf-core, Galaxy, MultiQC, Snakemake, or
other workflow integrations may count as records only after upstream release
or explicit external use evidence is public.

### 4. Industry impact

Reach ten permissioned records from distinct organizations, core facilities, or
production workflow owners. Each record must document the workflow scope, the
release used, the operational problem, and at least one measured before/after
metric such as runtime, memory, failure rate, QC visibility, or reproducibility.
Confidential data may remain private, but the existence and wording of the
record must be approved before publication. Downloads alone never satisfy this
gate.

## Evidence ledger

The target ledger should keep these fields separate:

- accuracy: dataset, truth source, reads checked, mismatches, ambiguity/outcome
  counts, and reproducibility command;
- performance: workload, comparator, hardware, versions, repeat statistics,
  throughput, wall time, peak memory, and output-equivalence status;
- scientific impact: external author, public URL, workflow, release, result,
  and scope boundary;
- industry impact: permission status, organization or anonymized identifier,
  workflow, release, before/after metric, and approved wording.

Do not combine the four dimensions into a single score until the underlying
records exist. A 10x claim is valid only when the corresponding gate above is
green and the evidence is independently auditable.

## Current next actions

1. Re-run the full Sanson/Brunello Hamming `k=1` lane five times and freeze the
   baseline artifact before optimizing.
2. Expand held-out oracle validation to ten lanes without changing semantics.
3. Convert open workflow submissions into accepted/released integrations before
   counting scientific impact.
4. Recruit permissioned scientific and industry evaluators through the public
   validation request; do not publish private data or infer impact from replies
   alone.
