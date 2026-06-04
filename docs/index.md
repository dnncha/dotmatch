# DotMatch Documentation

DotMatch is a deterministic command-line and Python toolkit for known-target
short-DNA assignment. It is designed for computational biologists and
bioinformaticians who already have a table of expected sequences and need to
count, demultiplex, audit, or diagnose reads without hiding ambiguous cases.

Use DotMatch when the biological question is:

> Which known guide, barcode, primer, feature tag, adapter, or panel target did
> this read contain?

DotMatch is intentionally narrower than a genome aligner, basecaller, UMI
pipeline, or screen-level statistics package. It works on extracted short
windows and known target lists. That narrow scope is what makes its assignment
contract easy to inspect: each read is reported as `unique`, `ambiguous`,
`none`, or `invalid`.

## Start Here

- New users should begin with [Getting Started](getting-started.md).
- CRISPR users can follow the [first-run CRISPR guide-counting tutorial](tutorials/crispr-count-first-run.md).
- Labs evaluating scientific claims should read [Trust, Scope, and Evidence](trust-and-scope.md).
- Workflow and pipeline authors should use the [public output schemas](schemas.md).

## Core Ideas

DotMatch compares a fixed read window with a known target table under explicit
edit-distance rules. By default, a read is counted only when exactly one target
falls inside the configured radius. If several targets are compatible, the read
is reported as ambiguous rather than assigned by accident.

This behavior matters in real assays. Unsafe one-mismatch correction, shifted
barcode positions, duplicate targets, low-quality rescued bases, and ambiguous
near-neighbors can all create plausible but wrong counts. DotMatch makes those
states visible in TSV, JSON, and HTML reports so results can be reviewed by
people and consumed by workflow systems.

```{toctree}
:maxdepth: 2
:caption: User Guide

getting-started
tutorials/crispr-count-first-run
assayspec
crispr-qc
barcode-panel-design
workbench
```

```{toctree}
:maxdepth: 2
:caption: Reference

schemas
methods-and-citation
joss-reviewer-guide
joss-author-confirmation
joss-submission
packaging
release-process
```

```{toctree}
:maxdepth: 2
:caption: Evidence and Boundaries

trust-and-scope
scientific-claims
barcode-science-readiness
usability-comparison
native-comparator-scope
benchmarks/README
evidence-gallery/README
```

```{toctree}
:maxdepth: 1
:caption: Detailed Reports
:hidden:

benchmarks/amplicon_panel/README
benchmarks/barcode_demux/README
benchmarks/barcode_panel_design/README
benchmarks/bcl_demux/README
benchmarks/crispr_comparison/README
benchmarks/feature_barcode/README
benchmarks/gpu/README
benchmarks/native/README
benchmarks/oligo_adapter/README
benchmarks/perturb_seq/README
benchmarks/public_crispr/README
benchmarks/real/README
evidence-gallery/report-zoo/README
evidence-gallery/scenarios/amplicon_artic_primer_start
evidence-gallery/scenarios/barcode_autopsy_review
evidence-gallery/scenarios/barcode_srp009896_comparator
evidence-gallery/scenarios/barcode_unsafe_correction
evidence-gallery/scenarios/barcode_wrong_offset_fixture
evidence-gallery/scenarios/bcl_tiny_classic
evidence-gallery/scenarios/feature_barcode_10x
evidence-gallery/scenarios/oligo_adapter_truseq_prefix
evidence-gallery/scenarios/perturb_seq_10x_guide_capture
evidence-gallery/scenarios/public_crispr_yusa
evidence-gallery/snapshots/barcode_autopsy/report
```
