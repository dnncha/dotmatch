# Barcode Validation Notes

This page describes what has been checked for DotMatch barcode work. The scope
is fixed-window known-target assignment after FASTQ generation. It does not
cover BCL conversion, adapter trimming, UMI/cell aggregation, or downstream
biological effect analysis.

Developer check:

```bash
make barcode-validation-ready
```

For users, the shortest useful path is `dotmatch barcode autopsy`: it produces
one `report.html` for review plus `findings.tsv`, `offset_scan.tsv`,
`correction_safety.tsv`, `top_unmatched.tsv`, and `provenance.json` for
pipeline records and methods review.

The check requires:

- at least five public fixed-window datasets;
- successful DotMatch rows with positive assignments;
- comparator or oracle rows with documented settings;
- zero recorded validation mismatches where validation is part of the row;
- metadata for each public dataset;
- plain notes on what each dataset does and does not support;
- explicit failure-mode fixtures for the barcode diagnostic report vocabulary.

Current public fixed-window examples are listed in
`docs/barcode-science-readiness.json`:

- SRP009896/SRR391079 inline barcode demultiplexing with Cutadapt and exact
  hash-splitter comparator rows;
- 10x TotalSeq-B feature barcode fixed-window assignment;
- 10x GEM-X CRISPR guide-capture fixed-window assignment;
- nf-core viralrecon ARTIC V3 primer-start fixed-window assignment;
- public TruSeq adapter-prefix fixed-window assignment.

These datasets answer different questions and should not be combined into one
broad biological claim. The wet-lab-facing report should explain whether a run
is clean, weakly specified, offset-shifted, collision-prone, low-quality,
ambiguous, invalid, or unmatched, then point the user to the file that supports
that conclusion.
