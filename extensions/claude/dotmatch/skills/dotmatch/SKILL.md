---
name: dotmatch
description: Assign FASTQ reads to a known list of short DNA targets with DotMatch. Use for CRISPR sgRNA guide counting (MAGeCK-compatible counts), inline barcode demultiplexing, feature-barcode or Perturb-seq guide-capture assignment, barcode or guide library collision audits, and short DNA edit distance. Not for genome alignment, UMI or cell quantification, variant calling, or screen statistics.
---

# DotMatch

DotMatch matches one fixed read window against a finite target table (CSV or TSV with id and sequence columns) using Hamming k<=3 or Levenshtein k<=2 correction. It runs locally and never uploads data.

## Pick the route

- Window position known (for example guide at read offset 23, length 20): call `audit_targets`, then `count_targets` or `demux_barcodes`.
- CRISPR screen or Perturb-seq direct capture with an unknown window: call `discover`, `prepare_assay`, `inspect_assay`, `run_assay`, `review_assay`, then `handoff_assay`, following each result's `next_actions`.
- Unsure whether DotMatch fits: call `route_task` with the task in plain words. If nothing matches, use a different tool.

These are tools of the `dotmatch` MCP server. Without MCP, run the same work through the CLI: `dotmatch count --help`, `dotmatch demux --help`, `dotmatch audit --help`, `dotmatch agent tools --json`, `dotmatch capabilities --json`. Install with `python3 -m pip install dotmatch` or `conda install -c bioconda dotmatch`.

## Rules

- Coordinates are zero-based.
- Before Hamming k=2 or k=3, require `safe_at_hamming_k2` or `safe_at_hamming_k3` to be true in the audit summary.
- Report only `unique` assignments as counts. Keep `ambiguous`, `none` and `invalid` separate and state their rates.
- Never edit target sequences or identities to raise assignment rates.
- Outputs must be new paths or empty directories; DotMatch refuses to overwrite.
- A low assignment rate usually means the wrong window. Try `prepare_assay` or `dotmatch inspect-unmatched` before raising k.
- Local tests and demos do not establish biological accuracy; say so when presenting results.
