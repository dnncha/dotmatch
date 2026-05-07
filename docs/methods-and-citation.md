# Methods and Citation Template

Use this page as a starting point for manuscripts, methods sections, benchmark notes, and software citations. Keep claims aligned with `docs/scientific-claims.md`.

## Software Citation

If you use DotMatch, cite the software release through `CITATION.cff`. Add the Zenodo DOI once a release DOI exists.

Suggested citation before DOI assignment:

> O'Toole D. DotMatch: Streaming Exact One-Edit Barcode and Guide Assignment Without Exhaustive Scanning. Software release v0.1.0. https://github.com/Dnncha/dotmatch

## Methods Sentence

For CRISPR guide-counting workflows:

> Reads were assigned to the guide library using DotMatch v0.1.0 with exact known-target assignment and deterministic `unique`, `ambiguous`, and `no-match` semantics. Count matrices retained only uniquely assigned reads; ambiguous and unmatched reads were excluded from target counts and retained in diagnostic summaries.

For one-edit Levenshtein rescue:

> DotMatch used global Levenshtein distance <=1 over the extracted guide window, including one-base substitutions, insertions, and deletions. Assignments were retained only when the best target was unique under the configured ambiguity policy.

For Hamming-only guide-counter-style comparisons:

> DotMatch used Hamming distance <=1 over fixed-length extracted guide sequences so that the comparison matched one-substitution/no-indel guide-counting semantics.

## Reproducibility Commands

Core verification:

```bash
make test
make cli-test
make python-test
make python-package-test
make publication-ready
make coverage
```

Current CRISPR claim gates:

```bash
make public-crispr-claim-gate
make crispr-sota-gate
```

Blocked broader claims:

```bash
make barcode-sota-gate
make bcl-sota-gate
```

The blocked gates are expected to fail until the required real-data and comparator evidence is present.

## Claim Boundary

DotMatch v0.1.0 should be described as a known-target short-DNA assignment engine. Do not describe it as a genome aligner, general Edlib replacement, or production Illumina demultiplexing replacement. The current public evidence supports narrow CRISPR guide-counting and exact short-DNA assignment claims; barcode and raw BCL state-of-the-art claims remain gated.
