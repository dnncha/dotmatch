# DotMatch Evidence Notes

Use this page when writing the README, website, or release notes. It says what
we can say plainly today, and where we still need more data before making a
broader claim.

The assay list lives in `docs/assay-evidence.json`. Run
`make assay-evidence-ready` after changing it. The check makes sure each lane
has the basics that matter for a public claim: artifacts, commands, comparator
semantics, validation notes, and a clear next step when the evidence is not
there yet. For raw CSV files, it also checks rows, command provenance, exit
codes where present, and recorded zero-mismatch validation columns.

## System Design Review

DotMatch is a deterministic known-target assignment system, not a learned
classifier. The current design uses exact edit-distance semantics, indexed
candidate generation, explicit ambiguity accounting, public-data benchmarks,
workflow fixtures, and gate scripts because those choices keep assignment
auditable and reproducible.

Applied techniques:

- packed A/C/G/T indexing for exact, Hamming `k=1`, Levenshtein `k=1`, and
  Levenshtein `k=2` candidate pruning where the fixed window can be encoded;
- exhaustive scan fallback for unsupported windows so the optimized path does
  not change assignment semantics;
- independent correctness checks against exhaustive scan and Edlib before
  performance rows are treated as evidence;
- radius ambiguity policy by default, so a read is counted only when exactly
  one target lies inside the configured edit-distance radius;
- Phred-quality gating for one-edit substitution and read-insertion rescue;
- barcode-panel error-sphere enumeration through `k=2`, nearest-neighbor
  checks, reverse-complement warnings, simulation, and machine-checkable
  assignment-collision reports;
- workflow-specific comparator rows for MAGeCK, guide-counter, Cutadapt-style
  barcode checks, exact-slice baselines, native Edlib scan, exact hash lookup,
  and BK-tree baselines where the semantics are comparable;
- experimental Apple Metal GPU benchmarking for packed fixed-window Hamming
  `k=1` workloads, with CPU checksum agreement required before speed ratios
  are reported, including a public CRISPR extract-pack-dispatch-readback-count
  lane;
- HTML/TSV/JSON QC outputs, MultiQC/workflow examples, and autopsy reports that
  expose offset errors, unsafe rescue, ambiguous reads, unmatched reads, and
  invalid extraction windows.

Techniques deliberately not claimed:

- supervised ML assignment is not currently used because the core workflow has
  a known target list and an explicit error model. A learned assignment layer
  would need labeled assay-specific training data, held-out public validation,
  and a clear advantage over the deterministic oracle before it could be
  claimed;
- calibrated statistical screen interpretation is left to downstream tools such
  as MAGeCK, BAGEL, drugZ, or CERES because DotMatch stops at read assignment
  and QC;
- probabilistic basecalling, UMI modeling, cell-level quantification, genome
  alignment, variant calling, adapter trimming, and production BCL conversion
  are outside the current evidence boundary.

## Current Defensible Statements

| Statement | Status | Evidence | Boundary |
| --- | --- | --- | --- |
| DotMatch provides exact short-DNA global edit distance and threshold matching for known targets. | Supported | `make test`, `make cli-test`, native C tests, Python tests, `docs/benchmarks/native/README.md` | Not a genome aligner; no CIGAR/SAM/BAM support. |
| The core Levenshtein `k=1` threshold predicate avoids heap allocation for exact, one-substitution, one-insertion, and one-deletion checks. | Supported | `make test`; native allocation regression in `tests/test_qdalign_threshold_alloc.c` | This describes the implementation of the threshold predicate. It does not make a cross-platform speed claim. |
| Indexed assignment preserves native exhaustive-scan semantics for `unique`, `ambiguous`, `none`, and `invalid` outcomes in the supported `k<=2` lanes. | Supported | `dotmatch validate`, native assignment tests, Edlib validation artifacts under `benchmarks/raw/`, Levenshtein `k=2` CLI regression cases | `k=2` uses packed A/C/G/T hash-neighborhood pruning for windows up to 32 bases, with fallback preserving semantics for unsupported cases; current `N`/IUPAC behavior is literal-byte matching, not wildcard expansion semantics. |
| Public CRISPR guide-counting rows are validated. | Supported | `make public-crispr-evidence-gate` passes; report at `docs/benchmarks/public_crispr/README.md` | Supports the documented MAGeCK/Yusa public-data workflow, not universal CRISPR superiority. |
| Extended CRISPR comparison rows are validated. | Supported | `make crispr-comparison-gate` passes; report at `docs/benchmarks/crispr_comparison/README.md` | Applies to the recorded CRISPR guide-counting lanes and their documented comparator semantics. |
| DotMatch has an experimental GPU acceleration evidence lane. | Experimental | `make bench-gpu`, `make gpu-evidence-gate`, report at `docs/benchmarks/gpu/README.md` | Current evidence is Apple Metal-only for packed Hamming `k=1`, including synthetic rows and a public CRISPR FASTQ extract-pack-dispatch-readback-count row. It is not a production speed claim. Promotion requires additional real-workload gates, non-A/C/G/T fallback, and production scheduling. |
| FASTQ count and demux workflows can optionally gate one-edit substitution and read-insertion rescue by observed Sanger Phred quality. | Supported | `make cli-test`; `--max-correction-qual` CLI regression cases | This is a deterministic correction filter, not a calibrated sequencing-error probability model. Read-deletion rescue has no observed edited base to score and is not rejected by this gate. |
| FASTQ count and demux workflows support optional Levenshtein `k=2` fixed-window correction. | Supported | `make cli-test`; Levenshtein `k=2` count/demux regression cases | Hamming remains limited to `k=0` and `k=1`; `--indel-window` remains a `k=1` option. |
| DotMatch has a first paired/combinatorial fixed-window counting command. | Supported | `make cli-test`; `pair-count` CLI regression case | Counts only reads where both target windows are uniquely assigned. This is not a perturb-seq expression quantification or guide-pair statistical-analysis workflow. |
| DotMatch has a native fixed-position inline barcode demultiplexing command with fixed-length and auto-length barcode sheet modes, and a checked public SRP009896 exact-prefix comparison lane. | Supported | `make cli-test`, `make bench-barcode-comparison`, `make barcode-comparison-gate`, report at `docs/benchmarks/barcode_demux/README.md` | Applies to the checked SRP009896/SRR391079 `k=0` exact-prefix lane with variable-length public example barcodes, Cutadapt anchored no-indel demux, and exact hash-splitter baseline. Broader barcode demultiplexing claims need additional datasets and comparator semantics. |
| DotMatch has a checked public lane-1 classic per-cycle tiny-BCL demultiplexing milestone. | Supported | `make cli-test`, `make bench-bcl-10x`, `make bcl-tiny-public-gate`, report at `docs/benchmarks/bcl_demux/README.md` | Applies only to the public 10x tiny-BCL classic per-cycle lane-1 row and count-total validation where available. Multi-lane operation, CBCL/NovaSeq support, and production Illumina demultiplexing replacement claims require separate evidence. |
| DotMatch has a checked public CRISPR Guide Capture single-guide extraction lane plus a perturb-seq-style guide/feature pair diagnostic lane. | Gated | `make bench-perturb-seq-public`, `make perturb-seq-public-gate`, report at `docs/benchmarks/perturb_seq/README.md` | The public row currently has one guide target, so it supports fixed-window extraction and exact-slice validation only. Useful multi-guide Perturb-seq assignment, Cell Ranger UMI/cell-level quantification, and perturbation-effect claims require additional evidence. |
| DotMatch has a checked public feature-barcode per-read assignment lane. | Supported | `make bench-feature-barcode-public`, `make feature-barcode-public-gate`, report at `docs/benchmarks/feature_barcode/README.md` | Applies to the checked 10x TotalSeq-B antibody Feature Barcode R2 fixed-window assignment lane and exact-slice baseline. This is not Cell Ranger UMI/cell-level quantification evidence. |
| DotMatch has a checked public amplicon/panel primer-start assignment lane. | Supported | `make bench-amplicon-panel-public`, `make amplicon-panel-public-gate`, report at `docs/benchmarks/amplicon_panel/README.md` | Applies to the checked nf-core viralrecon ARTIC V3 Illumina R1 fixed-window primer-start assignment lane and exact-prefix baseline, plus synthetic ambiguity diagnostics. This is not amplicon consensus, variant calling, primer trimming, or clinical validation evidence. |
| DotMatch has a checked public fixed-window oligo/adapter prefix assignment lane. | Supported | `make bench-oligo-adapter-public`, `make oligo-adapter-public-gate`, report at `docs/benchmarks/oligo_adapter/README.md` | Applies to the checked fast-adapter-trimming TruSeq R1 fixed-window adapter-prefix assignment lane and exact-slice baseline, plus synthetic ambiguity diagnostics. This is not adapter trimming, primer removal, UMI grouping, or read-merging evidence. |

## Evidence Boundaries

| Area | Check | Boundary |
| --- | --- | --- |
| Barcode demultiplexing comparisons. | `make barcode-comparison-gate` requires real-data rows, assigned reads, Cutadapt rows, and at least one additional relevant comparator. | The checked public lane is exact-prefix SRP009896/SRR391079 evidence; the built-in fixture remains only a workflow smoke test. |
| Raw BCL demultiplexing comparisons. | `make bcl-tiny-public-gate` checks the narrow public 10x tiny-BCL milestone; `make bcl-comparison-gate` requires real run folders, CBCL evidence where relevant, and production comparator validation before broader comparisons are made. | The public tiny-BCL row shows that the classic-BCL path works on the tiny fixture. It is not enough for broad BCL comparison evidence. |
| Perturb-seq comparisons. | `make perturb-seq-public-gate` requires the synthetic diagnostic lane plus public DotMatch k=0/k=1 guide-capture rows and exact-slice baseline agreement. | The checked public lane currently has one guide target, so it is a public extraction/validation milestone, not useful multi-guide assignment evidence. Guide-per-cell, expression-processing, and perturbation-effect claims require more data and comparator semantics. |
| Feature-barcode comparisons. | `make feature-barcode-public-gate` requires the synthetic diagnostic lane plus public DotMatch k=0/k=1 rows and exact-slice baseline agreement. | The checked public lane supports per-read fixed-window 10x antibody Feature Barcode assignment only. It does not support broader CITE-seq, cell-hashing, or cell/UMI quantification claims. |
| Amplicon/panel comparisons. | `make amplicon-panel-public-gate` requires the synthetic diagnostic lane plus public DotMatch k=0/k=1 primer-start rows and exact-prefix baseline agreement. | The checked public lane supports per-read fixed-window ARTIC primer-start assignment only. It does not support consensus generation, primer trimming, variant calling, or clinical interpretation claims. |
| Oligo/adapter comparisons. | `make oligo-adapter-public-gate` requires the synthetic diagnostic lane plus public DotMatch k=0/k=1 rows and exact-slice baseline agreement for the fast-adapter-trimming TruSeq R1 fixed window. | The checked public lane supports adapter-prefix assignment only; adapter trimming, primer removal, UMI grouping, read merging, and production adapter workflow claims require separate comparator evidence. |
| General aligner replacement. | No repository check promotes this. | DotMatch does not currently expose reference-index mapping, traceback/CIGAR, SAM/BAM, paired-end mapping, or genome-scale alignment semantics. |
| Native SeqAn/Parasail comparisons. | `make native-comparator-scope-ready` checks the documented scope in `docs/native-comparator-scope.md`. | SeqAn and Parasail are not completed comparator evidence until equivalent scoring semantics, native dependency/version capture, raw CSV rows, generated reports, and zero assignment mismatches are recorded. |
| Package-channel availability. | `make python-package-test` verifies local Linux/macOS wheel and sdist installability; Bioconda currently publishes DotMatch 0.1.4 as the latest public package before the 0.1.7 release, with 0.1.6 still under review. | Package-channel availability is distribution evidence only. It does not extend CRISPR, barcode, feature-barcode, BCL, amplicon, adapter, GPU, or aligner-comparison claims; `make distribution-channels` records each release channel only after public metadata and install smoke tests pass. |

## Rules For New Public Statements

Before adding a new README, website, or release-note claim, make sure it has:

- exact command lines;
- raw CSV artifacts under `benchmarks/raw/`;
- a generated report under `docs/benchmarks/`;
- correctness validation against the relevant oracle;
- comparator versions and semantics;
- a gate script when the statement is important enough to appear in the README.
