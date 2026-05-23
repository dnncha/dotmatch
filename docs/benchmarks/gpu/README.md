# Experimental GPU Acceleration Benchmark

This report is an experimental GPU evidence lane. It is intentionally not a production speed claim: the Metal path brute-forces packed Hamming `k=1` distances and is compared against DotMatch's existing CPU indexed Hamming assignment with identical output checks.

The decision rule is simple: GPU rows must have zero mismatches before any speed result is considered, and CPU-indexed throughput remains the production baseline unless the GPU path wins end-to-end on real workloads.

## Synthetic Figure

![GPU speedup](../../../benchmarks/figures/gpu_metal_speedup.svg)

## Public CRISPR Figure

![Public CRISPR GPU speedup](../../../benchmarks/figures/gpu_crispr_metal_speedup.svg)

## Synthetic CPU vs Metal Rows

| n_reads | n_targets | len | k | gpu_reads_per_sec | cpu_reads_per_sec | gpu_vs_cpu_kernel | gpu_vs_cpu_total | mismatches |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20000 | 737 | 20 | 1 | 4146940.0 | 297594.0 | 13.93x | 13.50x | 0 |
| 50000 | 4096 | 20 | 1 | 3649260.0 | 366657.0 | 9.95x | 10.05x | 0 |


## Public CRISPR CPU vs Metal Rows

| total_reads | packable_reads | n_targets | target_start | target_length | gpu_reads_per_sec | cpu_reads_per_sec | gpu_vs_cpu_kernel | gpu_vs_cpu_total | mismatches | count_delta |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20000 | 19954 | 87437 | 23 | 19 | 433598.0 | 191601.0 | 2.26x | 1.92x | 0 | 0 |


## Synthetic Raw Rows

| tool | backend | status | workload | n_reads | n_targets | len | k | prep_seconds | seconds | total_seconds | reads_per_sec | total_reads_per_sec | pairs_per_sec | checksum | mismatches | device | notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dotmatch_cpu_index | cpu | ok | synthetic_hamming | 20000 | 737 | 20 | 1 | 0.000604375 | 0.0672056 | 0.06781 | 297594 | 294942 | 2.19327e+08 | 124483641 | 0 | Apple M5 | CPU indexed Hamming k=1 baseline; prep is index build |
| dotmatch_gpu_metal | metal | ok | synthetic_hamming | 20000 | 737 | 20 | 1 | 0.0001995 | 0.00482283 | 0.00502233 | 4.14694e+06 | 3.98221e+06 | 3.05629e+09 | 124483641 | 0 | Apple M5 | Metal brute-force packed Hamming k=1; prep is shared-buffer allocation and copy |
| dotmatch_cpu_index | cpu | ok | synthetic_hamming | 50000 | 4096 | 20 | 1 | 0.00291654 | 0.136367 | 0.139284 | 366657 | 358979 | 1.50183e+09 | 1714574030 | 0 | Apple M5 | CPU indexed Hamming k=1 baseline; prep is index build |
| dotmatch_gpu_metal | metal | ok | synthetic_hamming | 50000 | 4096 | 20 | 1 | 0.000155292 | 0.0137014 | 0.0138567 | 3.64926e+06 | 3.60836e+06 | 1.49474e+10 | 1714574030 | 0 | Apple M5 | Metal brute-force packed Hamming k=1; prep is shared-buffer allocation and copy |


## Public CRISPR Raw Rows

| tool | backend | status | workload | total_reads | packable_reads | n_targets | target_start | target_length | k | input_seconds | prep_seconds | seconds | total_seconds | reads_per_sec | total_reads_per_sec | assigned_unique | assigned_exact | assigned_corrected | ambiguous | unmatched | invalid_windows | non_acgt_windows | checksum | mismatches | count_delta | device | notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dotmatch_cpu_index | cpu | ok | public_crispr_yusa_hamming | 20000 | 19954 | 87437 | 23 | 19 | 1 | 0.0747252 | 0.0536802 | 0.104384 | 0.232789 | 191601 | 85914.6 | 18376 | 17894 | 482 | 0 | 1578 | 0 | 46 | 15101300162 | 0 | 0 | Apple M5 | CPU indexed Hamming k=1 public CRISPR baseline |
| dotmatch_gpu_metal | metal | ok | public_crispr_yusa_hamming | 20000 | 19954 | 87437 | 23 | 19 | 1 | 0.0747252 | 0.000191333 | 0.0461257 | 0.121042 | 433598 | 165232 | 18376 | 17894 | 482 | 0 | 1578 | 0 | 46 | 15101300162 | 0 | 0 | Apple M5 | Metal public CRISPR FASTQ extract-pack-dispatch-readback-count lane |


## Scope

This lane tests whether GPU compute is worth productizing. It currently covers fixed-length packed A/C/G/T Hamming `k=1` assignment only. The public CRISPR row includes FASTQ parsing, guide-window extraction, packing, Metal dispatch, readback, and count/QC aggregation. It does not cover Levenshtein indels, BCL conversion, N/IUPAC GPU fallback, CUDA deployment, or production scheduling.

A future production GPU path needs additional real-workload gates for feature-barcode and BCL lanes, plus CPU fallback for non-A/C/G/T windows before promotion out of experimental status.
