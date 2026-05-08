# DotMatch Evidence Notes

This page maps user-facing statements to the evidence currently checked into the repository. Keep wording specific unless a raw artifact, validation command, and generated report support broader wording.

The structured assay-lane manifest is `docs/assay-evidence.json`; `make assay-evidence-ready` checks that supported, gated, smoke, and planned assay lanes have explicit artifacts, exact commands, comparator semantics, validation notes, or next-evidence requirements. For raw CSV artifacts named by that manifest, the gate also checks data rows, command and exit-code provenance where present, and zero validation mismatch columns where recorded.

## Current Defensible Statements

| Statement | Status | Evidence | Boundary |
| --- | --- | --- | --- |
| DotMatch provides exact short-DNA global edit distance and threshold matching for known targets. | Supported | `make test`, `make cli-test`, native C tests, Python tests, `docs/benchmarks/native/README.md` | Not a genome aligner; no CIGAR/SAM/BAM support. |
| The core Levenshtein `k=1` threshold predicate avoids heap allocation for exact, one-substitution, one-insertion, and one-deletion checks. | Supported | `make test`; native allocation regression in `tests/test_qdalign_threshold_alloc.c` | This is an implementation-level threshold predicate guarantee, not a cross-platform wall-clock speed claim. |
| Indexed assignment preserves native exhaustive-scan semantics for `unique`, `ambiguous`, `none`, and `invalid` outcomes in the supported `k<=2` lanes. | Supported | `dotmatch validate`, native assignment tests, Edlib validation artifacts under `benchmarks/raw/`, Levenshtein `k=2` CLI regression cases | `k=2` count/demux uses the exhaustive assignment fallback; current `N`/IUPAC behavior is literal-byte matching, not wildcard expansion semantics. |
| Public CRISPR guide-counting rows are validated. | Supported | `make public-crispr-evidence-gate` passes; report at `docs/benchmarks/public_crispr/README.md` | Supports the documented MAGeCK/Yusa public-data workflow, not universal CRISPR superiority. |
| Extended CRISPR comparison rows are validated. | Supported | `make crispr-comparison-gate` passes; report at `docs/benchmarks/crispr_comparison/README.md` | Applies to the recorded CRISPR guide-counting lanes and their documented comparator semantics. |
| FASTQ count and demux workflows can optionally gate one-edit substitution and read-insertion rescue by observed Sanger Phred quality. | Supported | `make cli-test`; `--max-correction-qual` CLI regression cases | This is a deterministic correction filter, not a calibrated sequencing-error probability model. Read-deletion rescue has no observed edited base to score and is not rejected by this gate. |
| FASTQ count and demux workflows support optional Levenshtein `k=2` fixed-window correction. | Supported | `make cli-test`; Levenshtein `k=2` count/demux regression cases | Hamming remains limited to `k=0` and `k=1`; `--indel-window` remains a `k=1` option. |
| DotMatch has a first paired/combinatorial fixed-window counting command. | Supported | `make cli-test`; `pair-count` CLI regression case | Counts only reads where both target windows are uniquely assigned. This is not a perturb-seq expression quantification or guide-pair statistical-analysis workflow. |
| DotMatch has a native fixed-position inline barcode demultiplexing command with fixed-length and auto-length barcode sheet modes, and a checked public SRP009896 exact-prefix comparison lane. | Supported | `make cli-test`, `make bench-barcode-comparison`, `make barcode-comparison-gate`, report at `docs/benchmarks/barcode_demux/README.md` | Applies to the checked SRP009896/SRR391079 `k=0` exact-prefix lane with variable-length public example barcodes, Cutadapt anchored no-indel demux, and exact hash-splitter baseline. Broader barcode demultiplexing claims need additional datasets and comparator semantics. |
| DotMatch has a checked public classic per-cycle tiny-BCL demultiplexing milestone. | Supported | `make cli-test`, `make bench-bcl-10x`, `make bcl-tiny-public-gate`, report at `docs/benchmarks/bcl_demux/README.md` | Applies only to the public 10x tiny-BCL classic per-cycle row and count-total validation where available. CBCL/NovaSeq-style and production Illumina replacement wording is blocked. |
| DotMatch has a checked public CRISPR Guide Capture per-read assignment lane plus a perturb-seq-style guide/feature pair diagnostic lane. | Supported | `make bench-perturb-seq-public`, `make perturb-seq-public-gate`, report at `docs/benchmarks/perturb_seq/README.md` | Applies to the checked 10x GEM-X A375 CRISPR Guide Capture R2 fixed-window assignment lane and exact-slice baseline, plus synthetic pair-count side diagnostics. This is not Cell Ranger UMI/cell-level quantification or perturbation-effect evidence. |
| DotMatch has a checked public feature-barcode per-read assignment lane. | Supported | `make bench-feature-barcode-public`, `make feature-barcode-public-gate`, report at `docs/benchmarks/feature_barcode/README.md` | Applies to the checked 10x TotalSeq-B antibody Feature Barcode R2 fixed-window assignment lane and exact-slice baseline. This is not Cell Ranger UMI/cell-level quantification evidence. |
| DotMatch has a checked public amplicon/panel primer-start assignment lane. | Supported | `make bench-amplicon-panel-public`, `make amplicon-panel-public-gate`, report at `docs/benchmarks/amplicon_panel/README.md` | Applies to the checked nf-core viralrecon ARTIC V3 Illumina R1 fixed-window primer-start assignment lane and exact-prefix baseline, plus synthetic ambiguity diagnostics. This is not amplicon consensus, variant calling, primer trimming, or clinical validation evidence. |
| DotMatch has a checked public fixed-window oligo/adapter prefix assignment lane. | Supported | `make bench-oligo-adapter-public`, `make oligo-adapter-public-gate`, report at `docs/benchmarks/oligo_adapter/README.md` | Applies to the checked fast-adapter-trimming TruSeq R1 fixed-window adapter-prefix assignment lane and exact-slice baseline, plus synthetic ambiguity diagnostics. This is not adapter trimming, primer removal, UMI grouping, or read-merging evidence. |

## Evidence Boundaries

| Area | Check | Boundary |
| --- | --- | --- |
| Barcode demultiplexing comparisons. | `make barcode-comparison-gate` requires real-data rows, assigned reads, Cutadapt rows, and at least one additional relevant comparator. | The checked public lane is exact-prefix SRP009896/SRR391079 evidence; the built-in fixture remains only a workflow smoke test. |
| Raw BCL demultiplexing comparisons. | `make bcl-tiny-public-gate` checks the narrow public 10x tiny-BCL milestone; `make bcl-comparison-gate` requires real run folders, CBCL evidence where relevant, and production comparator validation before broader wording is used. | The public tiny-BCL row is a classic-BCL parser/workflow milestone, not a broad BCL comparison result. |
| Perturb-seq comparisons. | `make perturb-seq-public-gate` requires the synthetic diagnostic lane plus public DotMatch k=0/k=1 guide-capture rows and exact-slice baseline agreement. | The checked public lane supports per-read fixed-window 10x CRISPR Guide Capture assignment only; guide-per-cell calls, expression processing, and perturbation-effect wording remain blocked. |
| Feature-barcode comparisons. | `make feature-barcode-public-gate` requires the synthetic diagnostic lane plus public DotMatch k=0/k=1 rows and exact-slice baseline agreement. | The checked public lane supports per-read fixed-window 10x antibody Feature Barcode assignment only; broader CITE-seq, cell-hashing, and cell/UMI quantification wording remains blocked. |
| Amplicon/panel comparisons. | `make amplicon-panel-public-gate` requires the synthetic diagnostic lane plus public DotMatch k=0/k=1 primer-start rows and exact-prefix baseline agreement. | The checked public lane supports per-read fixed-window ARTIC primer-start assignment only; consensus generation, primer trimming, variant calling, and clinical interpretation wording remain blocked. |
| Oligo/adapter comparisons. | `make oligo-adapter-public-gate` requires the synthetic diagnostic lane plus public DotMatch k=0/k=1 rows and exact-slice baseline agreement for the fast-adapter-trimming TruSeq R1 fixed window. | The checked public lane supports adapter-prefix assignment wording only; adapter trimming, primer removal, UMI grouping, read merging, and production adapter workflow claims require separate comparator evidence. |
| General aligner replacement. | No repository check promotes this. | DotMatch does not currently expose reference-index mapping, traceback/CIGAR, SAM/BAM, paired-end mapping, or genome-scale alignment semantics. |
| Native SeqAn/Parasail comparisons. | `make native-comparator-scope-ready` checks the documented scope in `docs/native-comparator-scope.md`. | SeqAn and Parasail are not completed comparator evidence until equivalent scoring semantics, native dependency/version capture, raw CSV rows, generated reports, and zero assignment mismatches are recorded. |
| Package-channel availability. | `make python-package-test` verifies local Linux/macOS wheel and sdist installability. | PyPI, Bioconda, Docker registry, and Zenodo release artifacts are separate distribution channels. |

## Rules For New Public Statements

New README, website, or release-note statements should include:

- exact command lines;
- raw CSV artifacts under `benchmarks/raw/`;
- a generated report under `docs/benchmarks/`;
- correctness validation against the relevant oracle;
- comparator versions and semantics;
- a gate script when the statement is important enough to appear in the README.
