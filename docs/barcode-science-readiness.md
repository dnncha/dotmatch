# Barcode Science Readiness

DotMatch's barcode story is intentionally narrower than general demultiplexing:
fixed-window known-target assignment after FASTQ generation. The readiness gate
therefore checks evidence for fixed-window biology, not BCL conversion, adapter
trimming, UMI/cell aggregation, or downstream biological effect analysis.

Run:

```bash
make barcode-science-ready
```

For users, the shortest useful path is `dotmatch barcode autopsy`: it produces a
single `report.html` for review plus `findings.tsv`, `offset_scan.tsv`,
`correction_safety.tsv`, `top_unmatched.tsv`, and `provenance.json` for
workflow and methods evidence.

The gate requires:

- at least five public fixed-window evidence datasets;
- successful DotMatch rows with positive assignments;
- comparator or oracle rows with documented semantics;
- zero recorded validation mismatches where validation is part of the row;
- evidence-ready metadata for each public dataset;
- conservative claim boundaries for every dataset;
- explicit failure-mode fixtures for the barcode-autopsy report vocabulary.

Current public fixed-window evidence lanes are listed in
`docs/barcode-science-readiness.json`:

- SRP009896/SRR391079 inline barcode demultiplexing with Cutadapt and exact
  hash-splitter comparator rows;
- 10x TotalSeq-B feature barcode fixed-window assignment;
- 10x GEM-X CRISPR guide-capture fixed-window assignment;
- nf-core viralrecon ARTIC V3 primer-start fixed-window assignment;
- public TruSeq adapter-prefix fixed-window assignment.

These datasets are not interchangeable biological claims. They are deliberately
used as separate fixed-window evidence lanes. The wet-lab-facing autopsy report
should explain whether a run is clean, weakly specified, offset-shifted,
collision-prone, low-quality, ambiguous, invalid, or unmatched, then point the
user to the exact artifact that supports that conclusion.
