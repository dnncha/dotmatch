# Oligo/Adapter Assignment Evidence

This report covers fixed-window assignment of short adapter-like oligos with DotMatch's known-target counting layer.

The synthetic lane checks exact, one-substitution, ambiguous, and unmatched adapter-like oligos. The public lane uses the fast-adapter-trimming TruSeq R1 fixture and validates DotMatch k=0 against a transparent exact-slice hash baseline over the documented fixed window.

## Synthetic Command

```bash
dotmatch count --targets benchmarks/work/oligo_adapter/adapter_oligos.tsv --reads benchmarks/work/oligo_adapter/adapter_reads.fastq --sample-label oligo_adapter_fixture --target-start 0 --target-length 12 --k 1 --metric hamming --format dotmatch --out benchmarks/work/oligo_adapter/adapter_counts.tsv --summary benchmarks/work/oligo_adapter/adapter_summary.json --assignments benchmarks/work/oligo_adapter/adapter_assignments.tsv --ambiguous report --sample-qc benchmarks/work/oligo_adapter/adapter_sample_qc.tsv
```

## Raw Rows

| tool | workflow | status | targets | reads | start | length | k | metric | assigned | exact | corrected | ambiguous | unmatched | validation mismatches |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| dotmatch_count | synthetic_oligo_adapter_fixture | smoke | 3 | 6 | 0 | 12 | 1 | hamming | 4 | 3 | 1 | 1 | 1 | 0 |
| dotmatch_count | public_fast_adapter_truseq_r1 | supported | 9 | 10000 | 229 | 20 | 0 | hamming | 156 | 156 | 0 | 0 | 9844 | 0 |
| dotmatch_count | public_fast_adapter_truseq_r1 | supported | 9 | 10000 | 229 | 20 | 1 | hamming | 159 | 156 | 3 | 0 | 9841 | 0 |
| exact_slice_hash | public_fast_adapter_truseq_r1 | supported | 9 | 10000 | 229 | 20 | 0 | exact | 156 | 156 | 0 | 0 | 9844 | 0 |

## Public Adapter-Prefix Lane

- Dataset: fast-adapter-trimming `788707_20180313_S_R1.small.fastq.gz` with `adapters/truseq.fa` target prefixes.
- Source repository: https://github.com/linsalrob/fast-adapter-trimming
- Source license: MIT, as reported by the upstream GitHub repository metadata.
- Fixed window: DotMatch uses `--target-start 229 --target-length 20` on R1.
- Comparator semantics: the exact-slice hash baseline counts reads whose fixed R1 substring exactly matches a deduplicated TruSeq adapter-prefix target. It validates fixed-window assignment semantics, not trimming correctness.

## Public Commands

```bash
dotmatch count --targets examples/oligo_adapter/data/adapter_oligos.tsv --reads examples/oligo_adapter/data/788707_20180313_S_R1.small.subsample10000.fastq.gz --sample-label fast_adapter_truseq_r1 --target-start 229 --target-length 20 --k 0 --metric hamming --format dotmatch --out benchmarks/work/oligo_adapter/public_fast_adapter_truseq_r1_k0_counts.tsv --summary benchmarks/work/oligo_adapter/public_fast_adapter_truseq_r1_k0_summary.json --assignments benchmarks/work/oligo_adapter/public_fast_adapter_truseq_r1_k0_assignments.tsv --ambiguous report --sample-qc benchmarks/work/oligo_adapter/public_fast_adapter_truseq_r1_k0_sample_qc.tsv
```

```bash
dotmatch count --targets examples/oligo_adapter/data/adapter_oligos.tsv --reads examples/oligo_adapter/data/788707_20180313_S_R1.small.subsample10000.fastq.gz --sample-label fast_adapter_truseq_r1 --target-start 229 --target-length 20 --k 1 --metric hamming --format dotmatch --out benchmarks/work/oligo_adapter/public_fast_adapter_truseq_r1_k1_counts.tsv --summary benchmarks/work/oligo_adapter/public_fast_adapter_truseq_r1_k1_summary.json --assignments benchmarks/work/oligo_adapter/public_fast_adapter_truseq_r1_k1_assignments.tsv --ambiguous report --sample-qc benchmarks/work/oligo_adapter/public_fast_adapter_truseq_r1_k1_sample_qc.tsv
```

```bash
python3 scripts/bench_oligo_adapter.py --include-public --metadata examples/oligo_adapter/data/metadata.json
```


## Evidence Boundary

Use these lanes only to verify fixed-window known-oligo/adapter assignment, one-substitution rescue, and explicit ambiguous/unmatched diagnostics. Run `make oligo-adapter-smoke-gate` for smoke evidence and `make oligo-adapter-public-gate` for the public lane. The public lane supports adapter-prefix assignment wording for the checked R1 window only; it is not adapter trimming evidence. Primer removal, UMI grouping, read merging, or production workflow claims require separate comparator semantics, raw artifacts, validation, and a passing gate.
