# Methods and citation

Use this page as a starting point for methods sections, benchmark notes, and software citations. Keep public statements aligned with `docs/scientific-claims.md`.

DotMatch uses a literal-byte alphabet policy for known-target assignment:
`N` and IUPAC ambiguity symbols are ordinary symbols and are not expanded as
wildcards. This policy is reported by `qdaln_alphabet_policy()` and recorded in
DotMatch count, demux, and pair-count summaries as `alphabet_policy`.
User-facing assignment workflows default to the `radius` ambiguity policy:
a read is `unique` only when exactly one target lies inside the configured
edit-distance radius. The `best` policy is available as an explicit
compatibility mode.

## Software Citation

If you use DotMatch, cite the software release through `CITATION.cff`.
Installed packages also provide `dotmatch citation` for a copyable citation.
Use the Zenodo concept DOI `10.5281/zenodo.20541628` for general software
citation. For an exact release, use its version-specific DOI: 0.4.0
(`10.5281/zenodo.22214073`), 0.5.0 (`10.5281/zenodo.22431034`), 0.6.0
(`10.5281/zenodo.22969176`), 0.6.1 (`10.5281/zenodo.22985881`), 0.6.3
(`10.5281/zenodo.22998373`), or 0.6.4 (`10.5281/zenodo.23043070`).

Published 0.6.4 packages and containers embed the concept DOI because its
version DOI was minted only after those immutable artifacts were published.
The 0.7.0 candidate has no minted version DOI; its source metadata defers
that identifier until publication. Both resolve to
the same Zenodo release family; use the version DOI when exact-release
provenance matters.

Suggested citation:

> O'Toole D. DotMatch: deterministic known-target short-DNA assignment for sequencing workflows. Software release v0.6.4. https://github.com/dnncha/dotmatch

Current release DOI (v0.6.4): <https://doi.org/10.5281/zenodo.23043070>

Concept DOI (all versions): <https://doi.org/10.5281/zenodo.20541628>

## Methods Sentence

For CRISPR guide-counting workflows:

> Reads were assigned to the guide library using DotMatch v0.6.4 with known-target assignment, literal-byte sequence semantics, and the radius ambiguity policy. Count matrices retained only reads for which exactly one guide lay inside the configured edit-distance radius; ambiguous and unmatched reads were excluded from target counts and retained in diagnostic summaries.

For one-edit Levenshtein rescue:

> DotMatch used global Levenshtein distance <=1 over the extracted guide window, including one-base substitutions, insertions, and deletions. Assignments were retained only when exactly one target was inside the configured radius, unless an explicit best-distance compatibility policy was recorded in the run summary.

For Hamming-only guide-counter-style comparisons:

> DotMatch used Hamming distance <=1 over fixed-length extracted guide sequences so that the comparison matched one-substitution/no-indel guide-counting semantics.

## Gene-level inference with CRISPRWorks Fit

Record counting and inference separately. Fit is an experimental source
package; retain `crisprworks-fit --version`, the checkout commit or source
hashes, numerical dependency versions, design or sample map, and control-gene
selection. The [Count-to-Fit tutorial](tutorials/crispr-fit-first-run.md)
describes the artifacts to keep.

For MLE, cite [MAGeCK MLE](https://doi.org/10.1186/s13059-015-0843-6) and record
MAGeCK2 0.3.0, normalization, guide-efficiency settings, permutation mode,
round count and seed. For the joint engine, cite
[JACKS](https://doi.org/10.1101/gr.238923.118) and record its shared guide-efficacy
model, matched baselines, pseudocount 32 and log2 effect scale.

If using gene-control calibration, record the control declarations, effect or
z score, exact guide-count strata, minimum control count, resolution, adjustment
and testing family. Cite the selected multiple-testing method:
[BH](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) or
[BY](https://doi.org/10.1214/aos/1013699998). Report directional or two-sided
testing as predeclared. The Count release DOI identifies its versioned counting
software; the new experimental Fit source needs its own precise provenance.

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

Current CRISPR evidence gates:

```bash
make public-crispr-evidence-gate
make crispr-comparison-gate
```

Current inline-barcode evidence gate:

```bash
make barcode-comparison-gate
```

Current feature-barcode assignment evidence gate:

```bash
make feature-barcode-public-gate
```

Current public CRISPR guide-capture assignment evidence gate:

```bash
make perturb-seq-public-gate
```

Current amplicon/panel primer-start assignment evidence gate:

```bash
make amplicon-panel-public-gate
```

Current public tiny-BCL milestone evidence gate:

```bash
make bcl-tiny-public-gate
```

Current oligo/adapter fixed-window public evidence gate:

```bash
make oligo-adapter-public-gate
```

Blocked broader comparisons:

```bash
make bcl-comparison-gate
```

These gates are deliberately narrow. The barcode gate is for the
SRP009896/SRR391079 exact-prefix lane. The feature-barcode gate is for the 10x
TotalSeq-B fixed-window per-read assignment lane, not Cell Ranger-style
UMI/cell quantification. The perturb-seq public gate is for the 10x CRISPR
Guide Capture fixed-window per-read assignment lane, not guide-per-cell calls,
expression quantification, or perturbation-effect analysis. The amplicon/panel
public gate is for the nf-core ARTIC V3 R1 fixed-window primer-start assignment
lane, not consensus generation, primer trimming, variant calling, clinical
panels, or diagnostic interpretation. The tiny-BCL public gate is for the
public 10x tiny-BCL classic per-cycle milestone, not production demultiplexing,
CBCL/NovaSeq support, or broad BCL comparison evidence. The oligo/adapter public
gate is for the fast-adapter-trimming TruSeq R1 fixed-window adapter-prefix
assignment lane, not adapter trimming, primer removal, UMI grouping, read
merging, or production adapter workflow evidence. Leave broader BCL comparison
statements out until real-data comparator evidence is in the repository.

## Evidence Boundary

Describe DotMatch v0.6.4 as a known-target short-DNA assignment engine. It is
not a genome aligner, general Edlib replacement, production Illumina
demultiplexer, full Perturb-seq analysis pipeline, adapter trimmer, UMI grouper,
read merger, or amplicon consensus/variant-calling workflow. Current public
evidence supports CRISPR guide counting, exact-prefix SRP009896/SRR391079
inline-barcode demultiplexing, per-read 10x TotalSeq-B feature-barcode
assignment, per-read 10x CRISPR Guide Capture assignment, nf-core ARTIC
amplicon primer-start assignment, the public 10x tiny-BCL classic per-cycle
milestone, fast-adapter-trimming TruSeq adapter-prefix assignment, and exact
short-DNA assignment statements for the documented workflows.
