# Perturb-seq direct-guide capture

Use intent `perturb-seq-guide-capture` for a finite guide-barcode list and guide-capture FASTQ files. DotMatch evaluates one reviewed fixed window and returns per-read `unique`, `ambiguous`, `none`, or `invalid` outcomes. Only unique assignments contribute to counts.

Required inputs:

- A local TSV or CSV table of known direct-guide sequences.
- A local directory containing the guide-capture FASTQ or FASTQ.gz files.
- An empty output directory.

The GSE146194 reference evidence checks per-read direct-guide assignment for 32 guides at exact and one-mismatch radii. It does not establish cell-barcode correction, UMI handling, guide-per-cell calls, expression processing, perturbation effects, or general performance on indirect capture designs. Keep those downstream steps and claims outside the DotMatch handoff.
