# Perturb-Seq And CRISPR Guide-Capture Assignment Evidence

This report covers perturb-seq-adjacent guide assignment evidence for DotMatch's known-target counting layer.

The synthetic lane checks fixed-window guide plus feature-barcode pair assignment through `pair-count`. The public lane uses a 10x Genomics CRISPR Guide Capture R2 subsample and validates DotMatch k=0 against a transparent exact-slice hash baseline over the observed fixed guide window.

Current status: public guide-capture assignment evidence only. This is not single-cell expression quantification or Cell Ranger perturbation-effect validation.

## Synthetic Command

```bash
dotmatch pair-count --left-targets benchmarks/work/perturb_seq/perturb_guides.tsv --right-targets benchmarks/work/perturb_seq/perturb_features.tsv --reads benchmarks/work/perturb_seq/perturb_reads.fastq --left-start 0 --left-length 6 --right-start 6 --right-length 6 --k 1 --metric hamming --out benchmarks/work/perturb_seq/perturb_pair_counts.tsv --summary benchmarks/work/perturb_seq/perturb_summary.json --assignments benchmarks/work/perturb_seq/perturb_assignments.tsv
```

## Synthetic Raw Row

| tool | workflow | guides | features | reads | k | metric | assigned pairs | ambiguous | left unmatched | right unmatched | invalid | validation mismatches |
| --- | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| dotmatch_pair_count | synthetic_perturb_seq_fixture | 3 | 2 | 7 | 1 | hamming | 3 | 1 | 1 | 1 | 1 | 0 |

## Public CRISPR Guide-Capture Rows

| tool | workflow | status | guides | reads | start | length | k | metric | assigned | exact | corrected | unmatched | validation mismatches |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| dotmatch_count | public_10x_crispr_guide_capture | supported | 1 | 20000 | 63 | 19 | 0 | hamming | 15979 | 15979 | 0 | 4021 | 0 |
| dotmatch_count | public_10x_crispr_guide_capture | supported | 1 | 20000 | 63 | 19 | 1 | hamming | 16582 | 15979 | 603 | 3418 | 0 |
| exact_slice_hash | public_10x_crispr_guide_capture | supported | 1 | 20000 | 63 | 19 | 0 | exact | 15979 | 15979 | 0 | 4021 | 0 |

## Public Dataset

- Dataset: 10x Genomics 1k A375 Cells Transduced with Non-Target and Target sgRNA, Chromium GEM-X Single Cell 5'.
- Source page: https://www.10xgenomics.com/datasets/1k-CRISPR-5p-gemx
- Fixture semantics: the fetcher selects the observed fixed-window CRISPR Guide Capture group with the most exact assignments in the copied R2 prefix.
- Comparator semantics: the exact-slice baseline counts reads whose fixed R2 substring exactly matches the selected guide sequence. It validates per-read guide assignment semantics, not Cell Ranger cell/UMI quantification or perturbation effects.

## Public Commands

```bash
dotmatch count --targets examples/perturb_seq/data/crispr_guides.tsv --reads examples/perturb_seq/data/1k_CRISPR_5p_gemx_crispr_S1_L001_R2.subsample20000.fastq.gz --sample-label 10x_crispr_guide_capture_L001_R2 --target-start 63 --target-length 19 --k 0 --metric hamming --format dotmatch --out benchmarks/work/perturb_seq/public_10x_crispr_k0_counts.tsv --summary benchmarks/work/perturb_seq/public_10x_crispr_k0_summary.json --assignments benchmarks/work/perturb_seq/public_10x_crispr_k0_assignments.tsv --ambiguous report --sample-qc benchmarks/work/perturb_seq/public_10x_crispr_k0_sample_qc.tsv
```

```bash
dotmatch count --targets examples/perturb_seq/data/crispr_guides.tsv --reads examples/perturb_seq/data/1k_CRISPR_5p_gemx_crispr_S1_L001_R2.subsample20000.fastq.gz --sample-label 10x_crispr_guide_capture_L001_R2 --target-start 63 --target-length 19 --k 1 --metric hamming --format dotmatch --out benchmarks/work/perturb_seq/public_10x_crispr_k1_counts.tsv --summary benchmarks/work/perturb_seq/public_10x_crispr_k1_summary.json --assignments benchmarks/work/perturb_seq/public_10x_crispr_k1_assignments.tsv --ambiguous report --sample-qc benchmarks/work/perturb_seq/public_10x_crispr_k1_sample_qc.tsv
```

```bash
python3 scripts/bench_perturb_seq.py --include-public --metadata examples/perturb_seq/data/metadata.json
```


## Evidence Boundary

Use these lanes to verify fixed-window guide/feature pair assignment, side-level diagnostics, and narrow public CRISPR guide-capture per-read assignment. Broader Perturb-seq comparisons require public cell barcode and UMI handling, guide-per-cell calls, expression or perturbation-effect comparator output, exact commands, validation artifacts, and a passing gate.
