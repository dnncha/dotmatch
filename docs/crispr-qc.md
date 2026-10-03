# CRISPR count QC

DotMatch CRISPR QC evaluates guide-counting and representation diagnostics for
a pooled CRISPR screen. It does not replace MAGeCK, BAGEL, drugZ, CERES,
CRISPResso2, or other downstream screen/editing analysis methods. Its job is
narrower: help detect likely FASTQ-to-guide counting configuration problems,
surface guide representation issues, and make suspicious samples visible before
downstream statistics.

## Command

```bash
dotmatch crispr qc \
  --counts counts.mageck.tsv \
  --sample-qc sample_qc.tsv \
  --library guides.csv \
  --out crispr_qc.json \
  --summary-tsv crispr_qc.summary.tsv \
  --report crispr_qc.html
```

`dotmatch assay run` writes these CRISPR QC artifacts automatically for
`mode = "count"` and `assay_type = "crispr"`.

The legacy-compatible alias `dotmatch crispr-qc` runs the same QC command.

## CRISPR-First Workflow

```bash
dotmatch crispr infer \
  --library guides.csv \
  --reads sample_R1.fastq.gz \
  --out assay.toml \
  --report inference_report.json

dotmatch crispr plan assay.toml
dotmatch crispr run assay.toml
```

These commands are thin wrappers over AssaySpec so the CRISPR interface and the
general workflow layer produce the same artifacts and validation behavior.

## Metrics

The report computes count-matrix representation metrics directly from the count
matrix instead of trusting precomputed values:

- total guide counts per sample;
- guide coverage and zero-count guide fraction;
- Gini index over guide counts;
- fraction of total assigned guide counts in the top 1% of guides;
- sample assignment, ambiguity, no-match, and invalid rates when
  `sample_qc.tsv` is provided;
- pairwise sample Pearson correlation on `log2(count + 1)`;
- pairwise Spearman correlation on raw guide counts;
- duplicate guide IDs and duplicate guide sequences;
- duplicate guide pairs and guide pairs within one edit;
- non-ACGT guide sequences.

When `sample_qc.tsv` is omitted, assignment, ambiguity, no-match, and invalid
rates are not evaluated and the report records a review warning. When the guide
library is omitted, duplicate/collision/sequence-content checks are not
evaluated and the report records a review warning.

For `--k 2` or higher, CRISPR QC still reports duplicate/one-edit guide
collisions only and records `safe_for_k = null` in the library summary. Use
`dotmatch audit --audit-mode exact` for target-collision auditing at larger
radii. Exact audit reports `safe_at_hamming_k2` and `safe_at_hamming_k3` for
same-length Hamming guide-count libraries; fast audit reports those fields as
not computed.

## Conservative Review Thresholds

DotMatch emits review warnings when any sample crosses these defaults:

```text
assignment_rate < 0.80
ambiguous_rate > 0.05
no_match_rate > 0.15
invalid_rate > 0.02
coverage_fraction < 0.90
zero_count_fraction > 0.10
gini_index > 0.50
top_1pct_fraction > 0.30
pairwise_sample_pearson < 0.80
```

These are diagnostic thresholds, not biological pass/fail laws. Screens differ
by library, cell model, selection pressure, PCR depth, and sampling strategy.
The report is deliberately conservative so suspicious guide-counting or sample
representation problems are reviewed before downstream modeling.
In `dotmatch assay start`, the pairwise correlation threshold only fails the
run for samples that share a `condition` in `assay.toml`; plasmid, T0 and
selected samples are expected to diverge.
`qc_status = "pass"` means no configured DotMatch QC threshold was crossed; it
does not certify screen quality or biological success.

## Scientific Boundary

After reviewing counting QC, the separate experimental CRISPRWorks Fit package
can consume the raw count matrix. See
[From counts to gene effects](tutorials/crispr-fit-first-run.md) for explicit
designs, joint sample maps and control-gene calibration. Keep counting QC and
the inference model's assumptions in the analysis record.

CRISPR QC provides diagnostics for the counting layer. It does not prove that a
count table is biologically correct, call enriched or depleted genes, infer
essentiality, quantify editing outcomes, evaluate off-target biology, or decide
whether a screen is biologically successful. Use MAGeCK-style or other
screen-analysis tools for hit calling and CRISPResso2-style tools for amplicon
editing outcomes.
