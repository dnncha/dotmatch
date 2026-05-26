# Barcode Troubleshooting Demo

This example shows how DotMatch investigates a fixed-window barcode run. It
uses the public SRP009896/SRR391079 inline-barcode fixture already used by the
barcode comparison check, then writes a QC report explaining the barcode window,
library safety, assignments, unmatched reads, and output files.

Run from the repository root:

```bash
make barcode-autopsy-demo
```

Start with `results/report.html`. The report gives a decision summary, findings,
plain-language interpretation, next actions, and exact commands. Use
`results/findings.tsv` for workflow automation and `results/provenance.json`
for command records.

The target writes:

- `results/report.html`
- `results/report.md`
- `results/findings.tsv`
- `results/offset_scan.tsv`
- `results/collision_graph.tsv`
- `results/correction_safety.tsv`
- `results/assignments.tsv`
- `results/sample_qc.tsv`
- `results/barcode_counts.tsv`
- `results/ambiguous.fastq`
- `results/unmatched.fastq`
- `results/top_unmatched.tsv`
- `results/provenance.json`
- `results/multiqc_dotmatch_barcode_mqc.yaml`

On the bundled SRP009896/SRR391079 fixture, the diagnostic output may mark the
offset inference as `review` because the highest exact-window signal is still
low. That is expected: the demo is meant to show how DotMatch surfaces failure
modes and unsafe rescue, not to turn a weak barcode signal into a silent
assignment.

The demo is intentionally scoped. It does not claim to replace Cutadapt or BCL
Convert. Cutadapt remains the appropriate general-purpose adapter/barcode
processing comparator; BCL Convert remains the appropriate BCL/CBCL-to-FASTQ
tool. DotMatch is for post-FASTQ fixed-window known-barcode assignment QC.
