# Barcode Demultiplexing Benchmark

This report records the checked barcode-demultiplexing example. It is separate from the CRISPR guide-counting report.

Current status: DotMatch has a checked public SRP009896/SRR391079 exact-prefix
inline-barcode example with five repeats, Cutadapt anchored no-indel demux rows,
and a simple exact-prefix hash splitter for comparison. This is a narrow
barcode example; broader barcode-demultiplexing statements need additional
datasets and comparator settings.

The benchmark script can also emit a simple `hash_splitter_exact` row. This is
an exact-prefix comparison, not an edit-distance demultiplexer.

## Figures

![Throughput](../../../benchmarks/figures/barcode_demux_throughput.svg)

![Peak memory](../../../benchmarks/figures/barcode_demux_peak_memory.svg)

![Assigned reads](../../../benchmarks/figures/barcode_demux_assigned_reads.svg)

![Verified candidates/read](../../../benchmarks/figures/barcode_demux_verified_per_read.svg)

## Raw Rows

| tool | workflow | assignment rules | repeats | reads | barcodes | k | metric | mean seconds | mean reads/sec | peak RSS KB | assigned | ambiguous | unmatched | verified/read | cv | exit |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| cutadapt_demux | real_srp009896_inline_barcode | anchored_cutadapt_demux_no_indels | 5 | 100000 | 48 | 0 | hamming | 0.578074 | 173055.7 | 23408 | 658 |  | 99342 |  | 0.0222 | 0 |
| dotmatch_demux | real_srp009896_inline_barcode | fixed_position_unique_ambiguous_nomatch | 5 | 100000 | 48 | 0 | hamming | 0.063247 | 1581617.3 | 5712 | 658 | 0 | 99342 | 0.0066 | 0.0203 | 0 |
| hash_splitter_exact | real_srp009896_inline_barcode | longest_unique_exact_prefix_no_mismatch | 5 | 100000 | 48 | 0 | exact | 0.308569 | 324920.1 |  | 658 |  | 99342 |  | 0.0570 | 0 |

## Checked Comparison

`make barcode-comparison-gate` passes for the SRP009896/SRR391079 exact-prefix
example shown above. The checked comparison is narrow: Cutadapt is run as an
anchored no-indel demultiplexer after trimming the leading `N`, and
`hash_splitter_exact` is a simple exact-prefix comparison, not an edit-distance
demultiplexer.

Suggested real-data starting point: SRP009896 / SRR391079-SRR391082, a maize GBS dataset described in public Cutadapt demultiplexing examples as 5-prime inline barcode reads with 96 demultiplexed outputs. `scripts/fetch_srp009896_barcode_demo.py --use-public-example-barcodes` extracts the first-member barcode sheet from the public Google Drive example archive with a ranged request instead of downloading the full 7.4 GB ZIP, then filters rows to the requested accession when the run column is present.

Important boundary: the SRP009896 barcode sheet contains variable-length barcodes (`4-8 bp`) and reused barcode sequences across run blocks. SRP009896 reads include a leading `N`, so the public-example benchmark should use `--barcode-start 1`, `--barcode-length auto`, and the exact-prefix `k=0` lane unless a separate fixed-length sheet is supplied.

The broader fixed-window barcode checks use the public fixed-window matrix in
`docs/barcode-science-readiness.json` and the failure-mode fixtures under
`examples/barcode_autopsy/failure_modes/`. Those examples support barcode
troubleshooting, but they do not turn this benchmark into a general BCL,
adapter-trimming, UMI/cell-quantification, or downstream biological-effect
comparison.
