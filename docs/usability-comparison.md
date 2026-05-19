# Why DotMatch For Known-Target Assignment

This table summarizes workflow fit and usability boundaries. It is not a
benchmark result. DotMatch is designed for runs where the target set and read
window are known up front and the assignment policy must remain auditable.

| Tool | Primary workflow | Direct FASTQ.gz input | Count matrix | One substitution | One insertion/deletion | Explicit ambiguous/no-match | Target audit | Offset/length diagnostics | Report output | Notes |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| DotMatch | known short-target assignment | yes | yes | yes | yes | yes | yes | yes | yes | General engine for guides, barcodes, panels, whitelists |
| guide-counter | CRISPR guide counting | yes | yes | yes | no, per current docs | workflow-specific | no | limited | no | Serious CRISPR comparator; compare directly for mismatch-only guide counting |
| MAGeCK count | CRISPR guide counting | yes | yes | exact FASTQ mode | no direct mismatch FASTQ route | limited | no | no | no | Downstream ecosystem standard |
| Cutadapt | adapter/search/trimming | yes | no | yes | yes | not assignment-centered | no | adapter-centered | no | Workflow comparator, not a direct assignment checker |
| Bowtie2 | reference alignment | yes | no | yes | yes | mapping-centered | no | mapping-centered | no | Over-general for known short-target assignment |
| Edlib scan | exact edit-distance checker | no workflow shell | no | yes | yes | yes if wrapped | no | no | no | Independent checker that compares reads against every target |

The practical difference is that DotMatch reports `unique`, `ambiguous`,
`none`, and `invalid` assignment outcomes under the same fixed-window rules used
for counting or demultiplexing. Barcode workflows can also produce collision,
offset, correction-safety, top-unmatched, and run-record outputs for review.

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
- exact Levenshtein matching rules including one-base indels;
- a hamming mode for fair one-mismatch/no-indel guide-counter comparisons;
- selected guide offset in the summary JSON when auto-offset detection is used;
- reproducible validation against an Edlib full scan.

## Scope Boundary

DotMatch is not a universal replacement for guide-counter. Its current CRISPR positioning is:

> Compared with mismatch-only guide counters, DotMatch provides exact Levenshtein matching with one-base insertion/deletion support, explicit ambiguous-read handling, target-library checks, Edlib validation, and the same known-target workflow across guides, barcodes, primers, and panels.

Direct speed comparisons against guide-counter require a pinned guide-counter version, exact commands, and a workflow where the matching rules being compared are clearly stated.
