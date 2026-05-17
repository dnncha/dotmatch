# Usability Comparison

This table summarizes workflow fit and usability boundaries. It is not a benchmark result.

| Tool | Primary workflow | Direct FASTQ.gz input | Count matrix | One substitution | One insertion/deletion | Explicit ambiguous/no-match | Target audit | Notes |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| DotMatch | known short-target assignment | yes | yes | yes | yes | yes | yes | General engine for guides, barcodes, panels, whitelists |
| guide-counter | CRISPR guide counting | yes | yes | yes | no, per current docs | workflow-specific | no | Serious CRISPR comparator; compare directly for mismatch-only guide counting |
| MAGeCK count | CRISPR guide counting | yes | yes | exact FASTQ mode | no direct mismatch FASTQ route | limited | no | Downstream ecosystem standard |
| Cutadapt | adapter/search/trimming | yes | no | yes | yes | not assignment-centered | no | Workflow comparator, not assignment oracle |
| Bowtie2 | reference alignment | yes | no | yes | yes | mapping-centered | no | Over-general for known short-target assignment |
| Edlib scan | exact pairwise oracle | no workflow shell | no | yes | yes | yes if wrapped | no | Exact semantic comparator; exhaustive over targets |

## Example Workflow

The target user-facing workflow is:

```bash
dotmatch count \
  --targets guides.csv \
  --reads sample.fastq.gz \
  --target-start 23 \
  --target-length 19 \
  --k 1 \
  --metric levenshtein \
  --indel-window 1 \
  --auto-offset 2 \
  --out counts.tsv \
  --summary summary.json \
  --ambiguous-out ambiguous.tsv \
  --unmatched-out unmatched.tsv
```

This should produce:

- count matrix for downstream analysis;
- deterministic assignment policy;
- ambiguity and unmatched diagnostics;
- exact Levenshtein semantics including one-base indels;
- a hamming mode for fair one-mismatch/no-indel guide-counter comparisons;
- selected guide offset in the summary JSON when auto-offset detection is used;
- reproducible validation against native Edlib scan.

## Scope Boundary

DotMatch is not a universal replacement for guide-counter. Its current CRISPR positioning is:

> Compared with mismatch-only guide counters, DotMatch provides a general exact Levenshtein assignment primitive with indel support, ambiguity semantics, target audit, native validation, and multi-domain known-target workflows.

Direct speed comparisons against guide-counter require a pinned guide-counter version, exact commands, and a workflow where the semantics being compared are clearly stated.
