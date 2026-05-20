# Public CRISPR Workflow Comparator

This report tracks the MAGeCK/Yusa public CRISPR benchmark. The single-run table below is a smoke/latest wiring check only; repeated rows and comparison-gated rows are the only rows intended to support user-facing performance statements.

## Smoke/Latest Wiring Table

**Reduced evidence.** These rows are secondary benchmark context. Use the repeated-run statistics below, and use `docs/benchmarks/crispr_comparison/README.md` once `make crispr-comparison-gate` passes for two real CRISPR datasets.

| tool | version | semantics | n_reads | n_targets | seconds | reads_per_sec | peak_rss_kb | assigned_reads | corrected_reads | ambiguous_reads | rejected_reads | overcount_reads | verified_per_read | exit_code |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dotmatch_exact_k0 | local | exact_k0_no_errors | 20000 | 87437 | 0.110329 | 181275.9 | 114448 | 17894 | 0 | 0 | 2106 | 0 | 0.8947 | 0 |
| dotmatch_levenshtein_k1 | local | levenshtein_k1_substitution_insertion_deletion | 20000 | 87437 | 3.402604 | 5877.9 | 117168 | 18552 | 1616 | 996 | 452 | 0 | 2.8277 | 0 |
| dotmatch_hamming_k1 | local | hamming_k1_no_indels | 20000 | 87437 | 0.202272 | 98876.7 | 118560 | 18261 | 522 | 156 | 1583 | 0 | 0.9288 | 0 |
| mageck_count_exact | 0.5.9.5 | exact_fastq_count_trim5_23 | 20000 | 87437 | 1.099061 | 18197.4 | 130576 | 17894 | 0 |  | 2106 | 0 |  | 0 |
| guide_counter_one_mismatch | 0.1.3 | hamming_k1_no_indels_auto_offset | 20000 | 87437 | 1.168958 | 17109.3 | 541360 | 20956 |  |  | 0 | 956 |  | 0 |
| external_competitors_ERR376998.fastq.gz | see_competitor_csv | cutadapt_bowtie2_extracted_workflow | 10000 | 87437 | 10.158057 | 984.4 | 1050160 |  |  |  |  | 0 |  | 0 |
| external_competitors_ERR376999.fastq.gz | see_competitor_csv | cutadapt_bowtie2_extracted_workflow | 10000 | 87437 | 8.913703 | 1121.9 | 1093248 |  |  |  |  | 0 |  | 0 |

![Public CRISPR throughput](../../../benchmarks/figures/public_crispr_throughput.svg)

![Public CRISPR runtime](../../../benchmarks/figures/public_crispr_runtime_seconds.svg)

![Public CRISPR memory](../../../benchmarks/figures/public_crispr_peak_memory.svg)

![Public CRISPR assignment impact](../../../benchmarks/figures/public_crispr_assignment_impact.svg)

![Public CRISPR verified candidates](../../../benchmarks/figures/public_crispr_verified_candidates.svg)

## Repeated-Run Statistics

| tool | semantics | records_per_sample | repeats | mean_reads_per_sec | p50_reads_per_sec | p95_reads_per_sec | mean_seconds | p50_seconds | cv | max_peak_rss_mb | mean_verified_per_read |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dotmatch_exact_k0 | exact_k0_no_errors | 10000 | 5 | 157549.3 | 159496.4 | 170724.6 | 0.1279 | 0.1254 | 0.0928 | 111.8 | 0.895 |
| dotmatch_exact_k0 | exact_k0_no_errors | 100000 | 5 | 891240.9 | 935485.4 | 1005640.6 | 0.2280 | 0.2138 | 0.1332 | 113.3 | 0.894 |
| dotmatch_hamming_k1 | hamming_k1_no_indels | 10000 | 5 | 98031.2 | 97890.4 | 100034.0 | 0.2040 | 0.2043 | 0.0138 | 115.8 | 0.929 |
| dotmatch_hamming_k1 | hamming_k1_no_indels | 100000 | 5 | 201672.3 | 210244.8 | 225625.7 | 1.0144 | 0.9513 | 0.1518 | 121.4 | 0.928 |
| dotmatch_levenshtein_k1 | levenshtein_k1_substitution_insertion_deletion | 10000 | 5 | 7023.6 | 7071.9 | 7137.8 | 2.8483 | 2.8281 | 0.0178 | 114.4 | 2.828 |
| dotmatch_levenshtein_k1 | levenshtein_k1_substitution_insertion_deletion | 100000 | 5 | 6726.7 | 6615.1 | 7237.6 | 29.7944 | 30.2338 | 0.0515 | 117.6 | 2.822 |
| guide_counter_one_mismatch | hamming_k1_no_indels_auto_offset | 10000 | 5 | 29200.3 | 29512.9 | 29966.5 | 0.6853 | 0.6777 | 0.0263 | 528.7 |  |
| guide_counter_one_mismatch | hamming_k1_no_indels_auto_offset | 100000 | 5 | 209503.8 | 222772.8 | 235285.7 | 0.9914 | 0.8978 | 0.1873 | 528.7 |  |
| mageck_count_exact | exact_fastq_count_trim5_23 | 10000 | 5 | 25454.9 | 26579.7 | 27042.9 | 0.7933 | 0.7525 | 0.1018 | 129.5 |  |
| mageck_count_exact | exact_fastq_count_trim5_23 | 100000 | 5 | 109483.7 | 112039.9 | 127978.7 | 1.8864 | 1.7851 | 0.1849 | 153.6 |  |

![Repeated public CRISPR throughput](../../../benchmarks/figures/public_crispr_repeated_throughput.svg)

![Repeated public CRISPR runtime](../../../benchmarks/figures/public_crispr_repeated_runtime_seconds.svg)

![Repeated public CRISPR peak memory](../../../benchmarks/figures/public_crispr_repeated_peak_memory.svg)

![Repeated public CRISPR verified candidates](../../../benchmarks/figures/public_crispr_repeated_verified_candidates.svg)

## DotMatch Hamming Speedup

This table keeps the fair CRISPR speed lane separate: DotMatch Hamming `k=1` versus tools with one-mismatch/no-indel or exact-count semantics.

| baseline | records_per_sample | dotmatch_hamming_reads_per_sec | baseline_reads_per_sec | speedup |
| --- | --- | --- | --- | --- |
| guide_counter_one_mismatch | 10000 | 98031.2 | 29200.3 | 3.36x |
| guide_counter_one_mismatch | 100000 | 201672.3 | 209503.8 | 0.96x |
| mageck_count_exact | 10000 | 98031.2 | 25454.9 | 3.85x |
| mageck_count_exact | 100000 | 201672.3 | 109483.7 | 1.84x |

## Count Agreement

| comparison | status | n_guides | total_left | total_right | total_delta | differing_guides | max_abs_delta | pearson | spearman |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dotmatch_hamming_vs_guide_counter | ok | 87437 | 18261 | 20956 | -2695 | 2409 | 7 | 0.93013799 | 0.93268449 |
| dotmatch_exact_vs_mageck_exact | ok | 87437 | 17894 | 17894 | 0 | 0 | 0 | 1.00000000 | 1.00000000 |

![Public CRISPR count agreement](../../../benchmarks/figures/public_crispr_count_agreement.svg)

## Multi-Sample Scaling

| tool | n_samples | records_per_sample | total_reads | threads | seconds | reads_per_sec | peak_rss_kb | assigned_reads | overcount_reads | exit_code |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dotmatch_hamming_k1_threaded | 2 | 100000 | 200000 | 2 | 0.498202 | 401443.9 | 124416 | 182657 | 0 | 0 |
| guide_counter_one_mismatch | 2 | 100000 | 200000 | 1 | 0.772011 | 259063.6 | 541376 | 208700 | 8700 | 0 |
| dotmatch_hamming_k1_threaded | 4 | 100000 | 400000 | 4 | 0.525465 | 761230.2 | 139696 | 365314 | 0 | 0 |
| guide_counter_one_mismatch | 4 | 100000 | 400000 | 1 | 1.055745 | 378879.3 | 541312 | 417400 | 17400 | 0 |
| dotmatch_hamming_k1_threaded | 8 | 100000 | 800000 | 8 | 0.551423 | 1450791.2 | 178560 | 730628 | 0 | 0 |
| guide_counter_one_mismatch | 8 | 100000 | 800000 | 1 | 1.640900 | 487537.2 | 541312 | 834800 | 34800 | 0 |

![Public CRISPR sample scaling throughput](../../../benchmarks/figures/public_crispr_sample_scaling_throughput.svg)

![Public CRISPR sample scaling memory](../../../benchmarks/figures/public_crispr_sample_scaling_memory.svg)

## Edlib Oracle Validation

| dataset | sample | oracle | checked_reads | mismatches | indel_window | stratum_exact | stratum_corrected | stratum_ambiguous | stratum_unmatched | stratum_contains_n |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mageck_yusa | plasmid | edlib_native | 1000 | 0 | 1 | 839 | 81 | 56 | 24 | 8 |
| mageck_yusa | ESC1 | edlib_native | 1000 | 0 | 1 | 848 | 77 | 50 | 25 | 5 |

## Interpretation

- `dotmatch_hamming_k1` is the fair lane for guide-counter-style one-mismatch/no-indel guide counting.
- `dotmatch_levenshtein_k1` is DotMatch's stronger lane: substitutions plus one-base insertions/deletions with explicit ambiguity reporting.
- `dotmatch_exact_k0` is the fair exact-count lane for MAGeCK's direct FASTQ counting mode.
- MAGeCK is run as exact FASTQ counting with `--trim-5 23`, matching the public Yusa demo workflow.
- guide-counter is fast, but on the 10k Yusa run its own stats report more mapped reads than input reads, consistent with its multi-offset counting loop; DotMatch assigns at most one target per read and reports ambiguity instead.
- In the multi-sample scaling table, DotMatch processes sample batches with threads while staying in the tens of MB. guide-counter uses roughly half a GB and its count total grows beyond input reads.
- Cutadapt and Bowtie2 rows are workflow comparators on extracted guide windows; they are not exact assignment oracles.
- Native Edlib scan remains the exact semantic oracle for assignment correctness.
- Public speed statements should cite only repeated rows with zero validation mismatches and explicit semantics.

## Raw Commands

| tool | command |
| --- | --- |
| dotmatch_exact_k0 | dotmatch count --targets examples/crispr_guides/data/yusa_library.csv --reads examples/crispr_guides/data/ERR376998.fastq.gz --reads examples/crispr_guides/data/ERR376999.fastq.gz --sample-label plasmid,ESC1 --target-start 23 --target-length 19 --k 0 --metric hamming --format mageck --out examples/crispr_guides/output/counts.exact.mageck.tsv --summary examples/crispr_guides/output/summary.exact.json |
| dotmatch_levenshtein_k1 | dotmatch count --targets examples/crispr_guides/data/yusa_library.csv --reads examples/crispr_guides/data/ERR376998.fastq.gz --reads examples/crispr_guides/data/ERR376999.fastq.gz --sample-label plasmid,ESC1 --target-start 23 --target-length 19 --k 1 --metric levenshtein --indel-window 1 --auto-offset 5 --auto-offset-sample 10000 --format mageck --out examples/crispr_guides/output/counts.levenshtein.mageck.tsv --summary examples/crispr_guides/output/summary.levenshtein.json |
| dotmatch_hamming_k1 | dotmatch count --targets examples/crispr_guides/data/yusa_library.csv --reads examples/crispr_guides/data/ERR376998.fastq.gz --reads examples/crispr_guides/data/ERR376999.fastq.gz --sample-label plasmid,ESC1 --target-start 23 --target-length 19 --k 1 --metric hamming --auto-offset 5 --auto-offset-sample 10000 --format mageck --out examples/crispr_guides/output/counts.hamming.mageck.tsv --summary examples/crispr_guides/output/summary.hamming.json |
| mageck_count_exact | build/competitor-env/bin/mageck count -l examples/crispr_guides/data/yusa_library.csv -n mageck_exact_benchmark --sample-label plasmid,ESC1 --trim-5 23 --fastq examples/crispr_guides/data/ERR376998.fastq.gz examples/crispr_guides/data/ERR376999.fastq.gz |
| guide_counter_one_mismatch | build/guide-counter/bin/guide-counter count --input examples/crispr_guides/data/ERR376998.fastq.gz examples/crispr_guides/data/ERR376999.fastq.gz --samples plasmid ESC1 --library examples/crispr_guides/data/yusa_library.csv --output examples/crispr_guides/output/guide_counter --offset-sample-size 10000 |
| external_competitors_ERR376998.fastq.gz | python3 scripts/bench_competitors.py --barcodes examples/crispr_guides/output/targets.tsv --reads examples/crispr_guides/data/ERR376998.fastq.gz --barcode-start 23 --barcode-length 19 --k 1 --dotmatch dotmatch --out examples/crispr_guides/output/competitors_ERR376998.fastq.csv --run-cutadapt --run-bowtie2 |
| external_competitors_ERR376999.fastq.gz | python3 scripts/bench_competitors.py --barcodes examples/crispr_guides/output/targets.tsv --reads examples/crispr_guides/data/ERR376999.fastq.gz --barcode-start 23 --barcode-length 19 --k 1 --dotmatch dotmatch --out examples/crispr_guides/output/competitors_ERR376999.fastq.csv --run-cutadapt --run-bowtie2 |
