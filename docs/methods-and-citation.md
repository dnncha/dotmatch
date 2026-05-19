# Methods and Citation Template

Use this page as a starting point for methods sections, benchmark notes, and software citations. Keep public statements aligned with `docs/scientific-claims.md`.

DotMatch uses a literal-byte alphabet policy for known-target assignment:
`N` and IUPAC ambiguity symbols are ordinary symbols and are not expanded as
wildcards. This policy is reported by `qdaln_alphabet_policy()` and recorded in
DotMatch count, demux, and pair-count summaries as `alphabet_policy`.

## Software Citation

If you use DotMatch, cite the software release through `CITATION.cff`. Add the Zenodo DOI once a release DOI exists.

Suggested citation before DOI assignment:

> O'Toole D. DotMatch: Streaming Exact One-Edit Barcode and Guide Assignment Without Exhaustive Scanning. Software release v0.1.0. https://github.com/dnncha/dotmatch

## Methods Sentence

For CRISPR guide-counting workflows:

> Reads were assigned to the guide library using DotMatch v0.1.0 with exact known-target assignment and deterministic `unique`, `ambiguous`, and `no-match` outcomes. Count matrices retained only uniquely assigned reads; ambiguous and unmatched reads were excluded from target counts and retained in diagnostic summaries.

For one-edit Levenshtein rescue:

> DotMatch used global Levenshtein distance <=1 over the extracted guide window, including one-base substitutions, insertions, and deletions. Assignments were retained only when the best target was unique under the configured ambiguity policy.

For Hamming-only guide-counter-style comparisons:

> DotMatch used Hamming distance <=1 over fixed-length extracted guide sequences so that the comparison matched one-substitution/no-indel guide-counting rules.

## Reproducibility Commands

Core verification:

```bash
make test
make cli-test
make python-test
make python-package-test
make repository-ready
make citation-metadata-ready
make workflow-examples-ready
make coverage
```

Current CRISPR evidence checks:

```bash
make public-crispr-evidence-gate
make crispr-comparison-gate
```

Current inline-barcode evidence check:

```bash
make barcode-comparison-gate
```

Current feature-barcode assignment evidence check:

```bash
make feature-barcode-public-gate
```

Current public CRISPR guide-capture assignment evidence check:

```bash
make perturb-seq-public-gate
```

Current amplicon/panel primer-start assignment evidence check:

```bash
make amplicon-panel-public-gate
```

Current public tiny-BCL milestone evidence check:

```bash
make bcl-tiny-public-gate
```

Current oligo/adapter fixed-window public evidence check:

```bash
make oligo-adapter-public-gate
```

Blocked broader comparisons:

```bash
make bcl-comparison-gate
```

These checks are deliberately narrow. The barcode check is for the
SRP009896/SRR391079 exact-prefix lane. The feature-barcode check is for the 10x
TotalSeq-B fixed-window per-read assignment lane, not Cell Ranger-style
UMI/cell quantification. The perturb-seq public check is for the 10x CRISPR
Guide Capture fixed-window per-read assignment lane, not guide-per-cell calls,
expression quantification, or perturbation-effect analysis. The amplicon/panel
public check is for the nf-core ARTIC V3 R1 fixed-window primer-start assignment
lane, not consensus generation, primer trimming, variant calling, clinical
panels, or diagnostic interpretation. The tiny-BCL public check is for the
public 10x tiny-BCL classic per-cycle milestone, not production demultiplexing,
CBCL/NovaSeq support, or broad BCL comparison evidence. The oligo/adapter public
check is for the fast-adapter-trimming TruSeq R1 fixed-window adapter-prefix
assignment lane, not adapter trimming, primer removal, UMI grouping, read
merging, or production adapter workflow evidence. Leave broader BCL comparison
statements out until real-data comparator evidence is in the repository.

## Evidence Boundary

Describe DotMatch v0.1.0 as a known-target short-DNA assignment engine. It is
not a genome aligner, general Edlib replacement, production Illumina
demultiplexer, full Perturb-seq analysis pipeline, adapter trimmer, UMI grouper,
read merger, or amplicon consensus/variant-calling workflow. Current public
evidence supports CRISPR guide counting, exact-prefix SRP009896/SRR391079
inline-barcode demultiplexing, per-read 10x TotalSeq-B feature-barcode
assignment, per-read 10x CRISPR Guide Capture assignment, nf-core ARTIC
amplicon primer-start assignment, the public 10x tiny-BCL classic per-cycle
milestone, fast-adapter-trimming TruSeq adapter-prefix assignment, and exact
short-DNA assignment statements for the documented workflows.
