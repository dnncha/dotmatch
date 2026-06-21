---
title: 'DotMatch: deterministic known-target short-DNA assignment for sequencing workflows'
tags:
  - bioinformatics
  - CRISPR
  - FASTQ
  - barcode demultiplexing
  - edit distance
  - sequence analysis
authors:
  - name: Donncha O'Toole
    orcid: 0009-0003-5012-7229
    corresponding: true
    affiliation: 1
affiliations:
  - name: Independent researcher, Ireland
    index: 1
date: 4 June 2026
bibliography: paper.bib
repository: https://github.com/dnncha/dotmatch
archive_doi: 10.5281/zenodo.20541628
---

# Summary

DotMatch is an open-source command-line and Python package for assigning short
DNA read windows to a known target table. In many sequencing assays the
expected sequences are already defined: for example a CRISPR guide library,
sample barcodes, feature barcodes, primer-prefix tags, or a designed panel of
assay identifiers. In these settings the immediate task is not genome
alignment, but deciding whether a fixed window in each read matches one of the
expected targets. DotMatch extracts that window, compares it with the supplied
target table under exact, Hamming, or Levenshtein edit-distance rules, and
reports each read as uniquely assigned, ambiguous, unmatched, or invalid.

The software provides native C routines, a command-line interface, Python
wrappers, workflow examples, and checks that connect documented claims to
reproducible examples. It is distributed from GitHub with release citation
metadata and a Zenodo archive [@dotmatch_zenodo_017]. The current examples cover
CRISPR guide counting, fixed-position barcode demultiplexing,
feature-barcode assignment, primer-start assignment, and barcode panel design.
DotMatch is deliberately narrower than a read aligner or a complete assay
analysis pipeline; it is intended to provide the assignment step that such
workflows often need before downstream analysis.

# Statement of need

Short-DNA assignment is a common step in sequencing analysis, but it is often
handled by a local script, a spreadsheet-derived lookup, or an option buried in
a larger pipeline. The practical question is simple: which known guide,
barcode, primer, or whitelist sequence is present at this fixed position in the
read? The details matter. A one-base correction may be acceptable for one target
table and unsafe for another; an apparently close read may be compatible with
two targets; an extraction window may fall outside a short read. These cases
should be visible in the output because they affect count matrices, sample
splits, and quality-control interpretation.

Existing tools address adjacent needs. MAGeCK provides statistical analysis for
CRISPR screens and includes guide-counting utilities [@li2014mageck]. Cutadapt
identifies and trims adapters and primers from reads [@martin2011cutadapt].
Edlib provides a fast exact edit-distance alignment library [@sosic2017edlib].
These tools are valuable, but they are aimed at different layers of the
workflow. DotMatch is focused on the smaller assignment problem: a known target
set, a read window, an edit-distance policy, and explicit handling of ambiguous
or invalid reads. This makes it useful for laboratories and workflow authors
who want a deterministic assignment layer before statistical or biological
interpretation.

# State of the field

The bioinformatics ecosystem already contains mature tools for read alignment,
adapter trimming, CRISPR screen analysis, and workflow distribution. DotMatch is
intended to sit beside these tools rather than replace them. It uses
edit-distance rules for short known targets, but it does not produce SAM/BAM
records, CIGAR strings, or reference-index mapping output. It can produce guide
counts, but downstream phenotype modeling remains the responsibility of
MAGeCK, BAGEL, drugZ, CERES, or related screen-analysis software. It can
classify a fixed adapter-prefix window, but it does not trim reads or replace
Cutadapt.

This limited scope is also a validation choice. The repository includes public
example data, generated reports, benchmark rows, and checks that constrain what
is claimed in documentation and release notes. Packaging work follows common
life-science software conventions, including Bioconda and container-oriented
release checks [@gruning2018bioconda].

# Software design

DotMatch treats assignment as a contract over four inputs: a target table, a
read window, an edit metric, and an ambiguity policy. Under the default radius
policy, a read is counted only when exactly one target lies within the
configured edit-distance radius. If more than one target is compatible, the read
is reported as ambiguous rather than forced to an arbitrary best match.
Unmatched reads and invalid extraction windows are reported separately so they
can be inspected during quality control.

The implementation combines native C kernels with Python and CLI interfaces.
For supported A/C/G/T windows, indexed candidate generation avoids scanning the
entire target table for every read while preserving the same assignment
semantics as exhaustive validation. Cases outside that indexed path fall back to
semantics-preserving code paths. Outputs include count tables, split FASTQ
files, target audits, top-unmatched tables, summary JSON, and HTML reports. A
barcode-panel design mode enumerates edit-distance neighborhoods, reports
nearest-neighbor risks, flags reverse-complement hazards, and writes collision
summaries that can be checked automatically.

Tests cover native edit-distance routines, CLI workflows, Python wrappers,
packaging metadata, workflow fixtures, and the public examples used in the
documentation. The project documentation separates supported statements from
experimental or unsupported claims, so performance and use-case language is
tied to reproducible artifacts rather than broad assertions.

# Research impact statement

DotMatch is intended for researchers building or auditing sequencing workflows
where known short sequences determine sample, guide, feature, primer, or panel
identity. Its immediate contribution is practical: reproducible assignment
decisions, explicit ambiguity accounting, and a documented assignment step that
can be cited or reused in methods descriptions. Public examples demonstrate
CRISPR guide-counting and barcode-style workloads using checked input data and
comparator semantics.

The project is early-stage, but it is already packaged with tests, release
metadata, and an archived DOI for citation. Future impact will depend less on
additional assignment features than on adoption in standard workflow contexts:
Bioconda/BioContainers propagation, nf-core-style modules, Galaxy wrappers, and
documented downstream uses.

# AI usage disclosure

OpenAI Codex was used to help prepare repository documentation, release
metadata, and this JOSS paper draft. The author remains responsible for the
software and manuscript. AI-generated edits were checked against repository
tests, evidence gates, release metadata, and public DOI records before being
retained.

# Acknowledgements

DotMatch uses public sequencing datasets, public CRISPR library resources, and
open-source bioinformatics tools as comparators and validation references. No
external financial support is claimed in this draft.

# References
