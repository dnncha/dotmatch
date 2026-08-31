# CRISPR guide counting

Use intent `crispr-guide-counting` with a finite guide library and a directory of local FASTQ or FASTQ.gz files. The prepared AssaySpec infers a fixed guide window, audits the library at the configured correction radius, preserves ambiguous assignments, writes per-sample QC and a MAGeCK-compatible count matrix, and applies CRISPR representation gates.

With `threads > 1`, rely on aggregate assignment and ambiguity counts in `summary.json` and `sample_qc.tsv`; ordered row-level diagnostic files require one thread.

Required inputs:

- A local TSV or CSV guide library with stable guide identifiers and sequences.
- A local directory containing the intended sample FASTQs.
- An empty output directory.

The evidence supports deterministic guide assignment and counting only. It does not establish screen hits, gene essentiality, differential abundance, biological replication quality beyond the reported gates, or downstream MAGeCK statistics. A target collision is a design or policy block; do not repair it by changing guide sequences.
