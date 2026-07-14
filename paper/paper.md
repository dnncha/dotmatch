---
title: 'DotMatch: ambiguity-aware known-target DNA assignment and auditable assay workflows'
tags:
  - bioinformatics
  - CRISPR
  - FASTQ
  - barcode demultiplexing
  - edit distance
  - assay reproducibility
authors:
  - name: Donncha O'Toole
    orcid: 0009-0003-5012-7229
    corresponding: true
    affiliation: 1
affiliations:
  - name: Independent researcher, Ireland
    index: 1
date: 14 July 2026
bibliography: paper.bib
repository: https://github.com/dnncha/dotmatch
archive_doi: 10.5281/zenodo.20541628
---

# Summary

DotMatch is open-source software for assigning short sequencing-read regions to
a known set of DNA targets. Such targets include CRISPR guides, inline sample
barcodes, feature barcodes, primer or adapter prefixes, and identifiers in
designed amplicon panels. These problems differ from general genome alignment:
the target set is known before analysis, the informative sequence is usually
short, and an incorrect rescue can directly contaminate a count matrix or
sample split.

DotMatch defines four assignment outcomes: unique, ambiguous, unmatched, and
invalid. Under its default radius policy a read is counted only when exactly one
target lies within the selected edit-distance radius. Multiple compatible
targets are reported as ambiguous rather than resolved by target order. Invalid
extraction windows and unmatched sequences remain visible in quality-control
artifacts.

The implementation combines native C distance and indexed-assignment kernels
with command-line and Python interfaces. It provides FASTQ counting,
demultiplexing, target-library audit, CRISPR-compatible count matrices,
barcode-panel design, unmatched-read inspection, assay inference, and
self-contained reliability reports. Reproducible public-data comparisons and
independent exhaustive or Edlib validation constrain performance and
correctness claims.

The same distribution now introduces AssayCode as an additive assay-level
identity. DotMatch remains the package, engine, scientific citation, and
compatibility contract. AssayCode supplies concise assay workflow commands and
an experimental AssayScript v2 compiler for multi-read segment descriptions.
The compiler validates target libraries and allowed combinations, fingerprints
inputs, records safety findings, and selects a deterministic execution strategy
for each segment. Experimental calibration and sequential-monitoring modules
are separated from default deterministic assignment until dedicated evidence
gates are satisfied.

# Statement of need

Known-target sequence assignment is frequently implemented as an undocumented
local script or embedded inside a larger analysis workflow. A simple nearest
match is not generally sufficient. Correction safety depends on the target
codebook, metric, allowed error radius, sequencing quality, extraction position,
and ambiguity policy. A barcode can be one substitution from several samples;
a guide window can be shifted by library construction; a read can be too short
to contain the requested region; independently plausible barcode components can
form a biologically impossible tuple.

# State of the field

These states affect downstream scientific conclusions. CRISPR screen analysis
software such as MAGeCK [@li2014mageck] expects reliable guide counts.
Single-cell CRISPR guide-assignment strategies can change the number of assigned
cells and downstream discoveries [@braunger2024crispat]. Demultiplexing systems
such as Pheniqs demonstrate the value of quality-aware probabilistic confidence
for complex barcode designs [@galanti2021pheniqs]. Cutadapt
[@martin2011cutadapt] and Flexiplex [@cheng2024flexiplex] address adjacent
adapter, barcode, or flexible sequence-search tasks. Edlib
[@sosic2017edlib] provides an exact edit-distance alignment library.

DotMatch does not seek to replace these tools. It provides a narrow,
inspectable layer between raw sequencing reads and downstream biological
analysis: extract declared regions, compare them with declared known targets,
preserve ambiguity, audit correction safety, and record enough provenance for
the decision to be reproduced.

# Software design

## Assignment contract and algorithms

The deterministic assignment contract consists of a target table, read window,
metric, maximum distance, and ambiguity policy. Exact, Hamming, and Levenshtein
modes are supported within documented command-specific bounds. The default
radius policy returns a unique assignment only when one target is compatible
within the complete radius. A separate best-distance policy is available for
explicit compatibility with workflows that choose the closest target.

For supported fixed A/C/G/T windows, indexed candidate generation reduces the
number of targets requiring distance verification. Exact hash lookup,
packed-neighbourhood indexes, seeded candidate generation, specialized Hamming
paths, and bit-parallel Myers distance kernels are selected according to metric,
length, and radius. Unsupported inputs fall back to semantics-preserving paths
rather than silently changing the assignment rule.

The implementation reports candidate-verification counts and can compare
indexed results against exhaustive native assignment or Edlib. Public evidence
gates require zero oracle disagreement for the validated lanes. Performance
claims are tied to raw benchmark artifacts, comparator semantics, software
versions, and generated reports instead of being generalized to all alignment
or demultiplexing workloads.

# Auditable assay workflows

AssaySpec v1 wraps assignment in a reviewable project. It validates referenced
inputs, compiles native commands, audits target libraries, executes counting or
demultiplexing, and writes normalized specifications, manifests, methods text,
software versions, reliability findings, suggested fixes, and HTML reports.
Production profiles can block unsafe correction or an unreviewed inferred
specification before assignment begins.

The autopsy workflow examines low assignment, ambiguity, invalid windows, and
frequent unmatched sequences. It can identify evidence for a wrong extraction
offset and propose a specification change. Proposed changes are written as
review artifacts; they are not silently applied.

Barcode-panel design uses the same assignment outcomes. It checks duplicates
and near neighbours, enumerates configured error spheres within documented
bounds, reports collision pairs, simulates reads, and exports lab-facing panel
and plate artifacts. This connects pre-sequencing codebook design with the
post-sequencing rules used to decode it.

# AssayCode and AssayScript

AssayCode is an additive product identity installed by the DotMatch Python
distribution. Existing commands, imports, native artifacts, schemas, DOI, and
citations remain supported. AssayCode provides assay-oriented shortcuts while
retaining an explicit escape hatch to the DotMatch engine.

AssayScript v2 is an experimental multi-segment description. A segment declares
its source read (R1, R2, I1, or I2), target library, fixed position or anchor,
length, positional jitter, orientation, metric, radius, ambiguity policy, and
whether it is required. A constraint table can enumerate allowed combinations
across segments.

The compiler currently produces a portable JSON plan rather than claiming a
complete universal assay runtime. The plan contains source and input SHA-256
fingerprints, target counts and lengths, safety status, selected matching
strategy, execution order, and review findings. This creates a testable
foundation for joint combinatorial execution without overstating current
capability.

# Experimental uncertainty and run monitoring

DotMatch's deterministic behavior remains the default. An experimental
quality-aware module fits per-cycle error rates and substitution patterns from
independently trusted observed/expected training pairs. It combines empirical
evidence with Phred probabilities, supports abundance priors with smoothing,
reports posterior mass and likelihood ratios, and abstains when configured
selective-decoding thresholds are not met.

A joint decoder combines calibrated component probability tables only over
declared allowed tuples. This can resolve evidence using assay constraints while
still returning ambiguous when posterior separation is insufficient.
Calibration metrics include Brier score and expected calibration error, and a
held-out threshold selector can maximize accepted calls subject to an empirical
false-discovery ceiling. These APIs remain experimental until public datasets
show calibration and yield improvements at a fixed measured error rate.

A deterministic design-time simulator perturbs fixed-length panels under a
declared scalar or per-cycle substitution model and applies the same
ambiguity-preserving radius rule. It reports usable yield, ambiguity, no-call
rate, false-discovery rate, and truth-by-call confusion. This digital twin is a
reproducible stress test for panel geometry, not a replacement for held-out
platform data.

AssayCode also includes a bounded-memory sequential monitor for assignment
events. It reports assignment, ambiguity, unmatched, and invalid rates with a
Wilson confidence interval and emits machine-readable on-track, review, or
insufficient-data decisions. It is a workflow monitoring primitive, not yet a
sequencer-control or production adaptive-sampling claim.

# Validation and evidence boundaries

Tests cover native kernels, deterministic fuzzing, Python APIs, CLI workflows,
AssaySpec projects, AssayScript compilation, calibration mathematics,
panel simulation, sequential monitoring, packaging, and public workflow fixtures. The repository
separates supported, experimental, and unsupported statements in a
machine-checked evidence inventory.

Public examples include CRISPR guide counting, inline barcode demultiplexing,
feature-barcode and guide-capture assignment, primer-prefix workloads, and
barcode-panel design. Comparative reports keep exact, Hamming, and Levenshtein
semantics separate. Accelerated or probabilistic paths are not promoted by
association with deterministic evidence; each requires its own correctness,
calibration, performance, and public-data gate.

# Research impact statement

DotMatch is distributed as source, Python wheels, and a Bioconda package
following common life-science distribution practices [@gruning2018bioconda].
The package includes the DotMatch command, Python APIs, native library
artifacts, and the additive AssayCode command. Release metadata and a Zenodo
archive support precise software citation [@dotmatch_zenodo_017].

The intended impact is a reusable reliability boundary for laboratories, core
facilities, workflow authors, and assay developers. A documented assignment
contract reduces the risk that plausible but ambiguous reads become
untraceable counts. The longer-term AssayCode direction connects assay
description, pre-run safety, deterministic execution, calibrated uncertainty,
diagnosis, and monitoring while preserving the narrower evidence boundary of
each implemented component.

# AI usage disclosure

OpenAI Codex was used to assist with implementation, tests, documentation,
benchmark infrastructure, and preparation of this manuscript. The author
remains responsible for the software, scientific claims, release decisions, and
manuscript. AI-assisted changes are retained only after repository tests and
the relevant evidence or release gates pass; experimental features are labeled
as such where those gates are incomplete.

# Acknowledgements

DotMatch uses public sequencing datasets, public CRISPR libraries, and
open-source tools as comparators and validation references. No external
financial support is claimed.

# References
