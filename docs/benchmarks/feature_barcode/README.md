# Feature Barcode Assignment Evidence

This report covers feature-barcode assignment evidence for DotMatch's known-target counting layer.

The synthetic lane checks exact, ambiguous, and unmatched feature IDs. The public lane uses a 10x Genomics TotalSeq-B antibody Feature Barcode R2 subsample and validates DotMatch k=0 against a transparent exact-slice hash baseline over the documented feature-reference window.

## Synthetic Command

```bash
dotmatch count --targets benchmarks/work/feature_barcode/feature_barcodes.tsv --reads benchmarks/work/feature_barcode/feature_reads.fastq --sample-label feature_barcode_fixture --target-start 0 --target-length 10 --k 1 --metric hamming --format dotmatch --out benchmarks/work/feature_barcode/feature_counts.tsv --summary benchmarks/work/feature_barcode/feature_summary.json --assignments benchmarks/work/feature_barcode/feature_assignments.tsv --ambiguous report --sample-qc benchmarks/work/feature_barcode/feature_sample_qc.tsv
```

## Raw Rows

| tool | workflow | status | features | reads | start | length | k | metric | assigned | exact | corrected | ambiguous | unmatched | validation mismatches |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| dotmatch_count | synthetic_feature_barcode_fixture | smoke | 4 | 6 | 0 | 10 | 1 | hamming | 4 | 4 | 0 | 1 | 1 | 0 |
| dotmatch_count | public_10x_totalseq_b_feature_barcode | supported | 10 | 20000 | 10 | 15 | 0 | hamming | 18575 | 18575 | 0 | 0 | 1425 | 0 |
| dotmatch_count | public_10x_totalseq_b_feature_barcode | supported | 10 | 20000 | 10 | 15 | 1 | hamming | 19101 | 18575 | 526 | 0 | 899 | 0 |
| exact_slice_hash | public_10x_totalseq_b_feature_barcode | supported | 10 | 20000 | 10 | 15 | 0 | exact | 18575 | 18575 | 0 | 0 | 1425 | 0 |

## Public Feature-Barcode Lane

- Dataset: 10x Genomics 1k Human PBMCs with TotalSeq-B Human TBNK Antibody Cocktail, 3' v3.1.
- Source page: https://www.10xgenomics.com/datasets/1-k-human-pbm-cs-with-total-seq-b-human-tbnk-antibody-cocktail-3-v-3-1-3-1-standard-6-0-0
- Feature reference pattern: `^NNNNNNNNNN(BC)NNNNNNNNN`, so DotMatch uses `--target-start 10 --target-length 15` on antibody R2.
- Comparator semantics: the exact-slice baseline counts reads whose fixed R2 substring exactly matches a feature-reference sequence. It validates per-read assignment semantics, not Cell Ranger cell/UMI quantification.

## Public Commands

```bash
dotmatch count --targets examples/feature_barcode/data/feature_barcodes.tsv --reads examples/feature_barcode/data/1k_PBMCs_TotalSeq_B_3p_antibody_S5_L001_R2.subsample20000.fastq.gz --sample-label 10x_totalseq_b_antibody_L001_R2 --target-start 10 --target-length 15 --k 0 --metric hamming --format dotmatch --out benchmarks/work/feature_barcode/public_10x_totalseq_b_k0_counts.tsv --summary benchmarks/work/feature_barcode/public_10x_totalseq_b_k0_summary.json --assignments benchmarks/work/feature_barcode/public_10x_totalseq_b_k0_assignments.tsv --ambiguous report --sample-qc benchmarks/work/feature_barcode/public_10x_totalseq_b_k0_sample_qc.tsv
```

```bash
dotmatch count --targets examples/feature_barcode/data/feature_barcodes.tsv --reads examples/feature_barcode/data/1k_PBMCs_TotalSeq_B_3p_antibody_S5_L001_R2.subsample20000.fastq.gz --sample-label 10x_totalseq_b_antibody_L001_R2 --target-start 10 --target-length 15 --k 1 --metric hamming --format dotmatch --out benchmarks/work/feature_barcode/public_10x_totalseq_b_k1_counts.tsv --summary benchmarks/work/feature_barcode/public_10x_totalseq_b_k1_summary.json --assignments benchmarks/work/feature_barcode/public_10x_totalseq_b_k1_assignments.tsv --ambiguous report --sample-qc benchmarks/work/feature_barcode/public_10x_totalseq_b_k1_sample_qc.tsv
```

```bash
python3 scripts/bench_feature_barcode.py --include-public --metadata examples/feature_barcode/data/metadata.json
```


## Evidence Boundary

Use these lanes to verify fixed-window feature-barcode whitelist counting and explicit ambiguity handling. The public 10x lane supports only per-read assignment against the documented Feature Barcode reference window. Broader CITE-seq, cell-hashing, or cell-level quantification comparisons require public comparator output, UMI/cell aggregation validation, exact commands, and a passing gate.
