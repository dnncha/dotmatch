# Barcode Troubleshooting Report

DotMatch reports how fixed-window barcode reads were assigned, rejected, or flagged for review.

Speed is reported only after the comparator settings are documented.

## Decision Summary

- Inference status: `review`
- Assignment rate: `0.72492000`
- Top finding: `low_confidence_offset`
- Primary report: `report.html`

## Comparator check

Use this section with recorded comparator settings for the same barcode window, length, and correction policy.

## Offset check

Highest-scoring sampled barcode window: start=0, length=4.
Exact assignment rate at that window: 0.11283. Inference status: review.
Warning: best scanned window has a low exact assignment rate; review the barcode sheet, read orientation, and window before treating the offset as ready.

## Barcode safety audit

The audit output reports duplicate and nearby barcode pairs before one-edit correction is trusted.

## Assignment summary

Run assignment rate: 0.72492000.
Exact assignments: 11283. Corrected assignments: 61209.
Ambiguous reads: 13911 (0.13911000). Unmatched reads: 13597 (0.13597000). Invalid windows: 0.
Top failure reason: ambiguous_collision. Top unmatched sequences are written to `top_unmatched.tsv`.

### What this means

The sampled fixed-window scan found a best window, but the evidence is weak enough that the assay specification should be reviewed before production use.

### Next action

Check the barcode sheet, read side, offset, barcode length, and expected orientation; rerun inference after correcting the assay specification.

Do not rescue ambiguous reads into either sample without changing the barcode design or assignment policy.

## Findings

- review: low_confidence_offset - The sampled scan has a weak best window; the assay specification may be wrong or incomplete. Next action: Review read side, barcode start, barcode length, target sheet, and orientation before production use.
- review: low_assignment_rate - Fewer than 80% of reads were uniquely assigned under the selected fixed-window rules. Next action: Inspect offset_scan.tsv, top_unmatched.tsv, and the barcode sheet before trusting rescued reads.
- review: high_ambiguity_rate - A material fraction of reads is compatible with more than one barcode. Next action: Do not rescue ambiguous reads into either sample without changing the barcode design or assignment policy.
- review: unsafe_one_edit_correction - At least one barcode is not safe for one-edit correction. Next action: Use k=0 or redesign/fix colliding barcodes before enabling one-edit rescue.

## Workflow handoff

Artifacts are written as stable TSV, JSON, FASTQ, HTML, and MultiQC custom-content inputs.

## QC Checklist

- Exact command records are stored in `provenance.json`.
- Offset evidence is recorded in `offset_scan.tsv`.
- Barcode collision safety is recorded under `audit/` and summarized in `correction_safety.tsv`.
- Ambiguous and unmatched reads are retained when requested instead of being silently assigned.
- Benchmark notes stay tied to documented comparator settings.

## Commands

- `dotmatch audit --targets examples/barcode_demux/data/barcodes.tsv --k 1 --audit-mode auto --out-dir examples/barcode_autopsy/results/audit`
- `dotmatch demux --barcodes examples/barcode_demux/data/barcodes.tsv --reads examples/barcode_demux/data/SRR391079.subsample100000.fastq.gz --barcode-start 0 --barcode-length auto --k 1 --metric hamming --out-dir examples/barcode_autopsy/results/demuxed --summary examples/barcode_autopsy/results/summary.json --assignments examples/barcode_autopsy/results/assignments.tsv --ambiguous-out examples/barcode_autopsy/results/ambiguous.fastq --unmatched-out examples/barcode_autopsy/results/unmatched.fastq`
- `dotmatch inspect-unmatched --targets examples/barcode_demux/data/barcodes.tsv --reads examples/barcode_demux/data/SRR391079.subsample100000.fastq.gz --target-start 0 --target-length 4 --k 1 --top 100 --out examples/barcode_autopsy/results/top_unmatched.tsv`

## Artifacts

- `offset_scan.tsv`
- `audit/`
- `collision_graph.tsv`
- `correction_safety.tsv`
- `demuxed/`
- `summary.json`
- `assignments.tsv`
- `ambiguous.fastq`
- `unmatched.fastq`
- `top_unmatched.tsv`
- `sample_qc.tsv`
- `barcode_counts.tsv`
- `findings.tsv`
- `provenance.json`
- `multiqc_dotmatch_barcode_mqc.yaml`
