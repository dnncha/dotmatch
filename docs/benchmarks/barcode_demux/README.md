# Barcode Demultiplexing Benchmark

This report is the barcode-demultiplexing evidence track. It is separate from the CRISPR guide-counting report.

Current status: DotMatch has a checked public SRP009896/SRR391079 exact-prefix inline-barcode lane with five repeats, Cutadapt anchored no-indel demux rows, and an exact hash-splitter baseline. Broader barcode-demultiplexing claims require additional datasets and comparator semantics, not only this public lane.

The benchmark script can also emit a simple `hash_splitter_exact` row. This is a transparent exact-prefix baseline, not an edit-distance demultiplexer.

## Figures

![Throughput](../../../benchmarks/figures/barcode_demux_throughput.svg)

![Peak memory](../../../benchmarks/figures/barcode_demux_peak_memory.svg)

![Assigned reads](../../../benchmarks/figures/barcode_demux_assigned_reads.svg)

![Verified candidates/read](../../../benchmarks/figures/barcode_demux_verified_per_read.svg)

## Raw Rows

| tool | workflow | semantics | repeats | reads | barcodes | k | metric | mean seconds | mean reads/sec | peak RSS KB | assigned | ambiguous | unmatched | verified/read | cv | exit |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| cutadapt_demux | real_srp009896_inline_barcode | anchored_cutadapt_demux_no_indels | 5 | 100000 | 48 | 0 | hamming | 1.228280 | 95052.5 | 19296 | 658 |  | 99342 |  | 0.3347 | 0 |
| dotmatch_demux | real_srp009896_inline_barcode | fixed_position_unique_ambiguous_nomatch | 5 | 100000 | 48 | 0 | hamming | 0.073314 | 1373183.2 | 5744 | 658 | 0 | 99342 | 0.0066 | 0.0900 | 0 |
| hash_splitter_exact | real_srp009896_inline_barcode | longest_unique_exact_prefix_no_mismatch | 5 | 100000 | 48 | 0 | exact | 0.382510 | 270395.2 |  | 658 |  | 99342 |  | 0.2028 | 0 |

## Comparison Evidence Gate

This table passes `make barcode-comparison-gate` for the checked SRP009896/SRR391079 exact-prefix lane. Keep wording specific: Cutadapt is run as anchored no-indel demux after trimming the leading `N`, and `hash_splitter_exact` is a transparent exact-prefix baseline, not an edit-distance demultiplexer.

Suggested real-data starting point: SRP009896 / SRR391079-SRR391082, a maize GBS dataset described in public Cutadapt demultiplexing examples as 5-prime inline barcode reads with 96 demultiplexed outputs. `scripts/fetch_srp009896_barcode_demo.py --use-public-example-barcodes` extracts the first-member barcode sheet from the public Google Drive example archive with a ranged request instead of downloading the full 7.4 GB ZIP, then filters rows to the requested accession when the run column is present.

Important boundary: the SRP009896 barcode sheet contains variable-length barcodes (`4-8 bp`) and reused barcode sequences across run blocks. SRP009896 reads include a leading `N`, so the public-example benchmark should use `--barcode-start 1`, `--barcode-length auto`, and the exact-prefix `k=0` lane unless a separate fixed-length sheet is supplied.
