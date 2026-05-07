# DotMatch Scientific Claim Ledger

This ledger maps public-facing claims to the evidence that currently supports or blocks them. It is deliberately conservative: a claim should stay narrow until a gate, raw artifact, and reproducible command support it.

## Current Defensible Claims

| Claim | Status | Evidence | Boundary |
| --- | --- | --- | --- |
| DotMatch provides exact short-DNA global edit distance and threshold matching for known targets. | Supported | `make test`, `make cli-test`, native C tests, Python tests, `docs/benchmarks/native/README.md` | Not a genome aligner; no CIGAR/SAM/BAM claim. |
| Indexed assignment preserves native exhaustive-scan semantics for `unique`, `ambiguous`, `none`, and `invalid` outcomes in the supported `k<=1` lanes. | Supported | `dotmatch validate`, native assignment tests, Edlib validation artifacts under `benchmarks/raw/` | Current wildcard `N` behavior is literal-byte matching, not IUPAC wildcard semantics. |
| Public CRISPR guide-counting claim gate passes on checked artifacts. | Supported | `make public-crispr-claim-gate` returns `PUBLIC CRISPR CLAIM GATE: PASS`; report at `docs/benchmarks/public_crispr/README.md` | Supports a narrow MAGeCK/Yusa public-data workflow claim, not universal CRISPR superiority. |
| Two-dataset CRISPR SOTA gate passes on checked artifacts. | Supported | `make crispr-sota-gate` returns `CRISPR SOTA GATE: PASS`; report at `docs/benchmarks/crispr_sota/README.md` | Applies to the recorded CRISPR guide-counting lanes and their documented comparator semantics. |
| DotMatch has a native fixed-position inline barcode demultiplexing command. | Supported | `make cli-test`, `make bench-barcode-demux`, report at `docs/benchmarks/barcode_demux/README.md` | Barcode SOTA claims are blocked until real barcode sheet and comparator requirements pass. |
| DotMatch has a first classic per-cycle BCL demultiplexing milestone. | Supported | `make cli-test`, `make bench-bcl-small`, report at `docs/benchmarks/bcl_demux/README.md` | CBCL/NovaSeq-style and production Illumina replacement claims are blocked. |

## Blocked Or Not-Yet Claims

| Claim | Current Gate | Why Blocked |
| --- | --- | --- |
| Barcode demultiplexing state of the art. | `make barcode-sota-gate` currently fails. | The checked artifact lacks a real barcode sheet, real repeated rows, Cutadapt repeated rows, and a second comparator such as Ultraplex, Je, or exact hash splitter. |
| Raw BCL demultiplexing state of the art. | `make bcl-sota-gate` currently fails. | The checked artifacts do not include a successful CBCL row and production comparator validation. |
| General aligner replacement. | No gate should promote this. | DotMatch does not currently expose reference-index mapping, traceback/CIGAR, SAM/BAM, paired-end mapping, or genome-scale alignment semantics. |
| Published package-channel availability. | Packaging verifier present, external channels not released. | `make python-package-test` verifies local Linux/macOS wheel and sdist installability, but PyPI, Bioconda, Docker registry, and Zenodo release artifacts are not published yet. |

## Rules For New Claims

New README, website, manuscript, or release-note claims should include:

- exact command lines;
- raw CSV artifacts under `benchmarks/raw/`;
- a generated report under `docs/benchmarks/`;
- correctness validation against the relevant oracle;
- comparator versions and semantics;
- a gate script when the claim is important enough to appear in the README headline.
