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
archive_doi: 10.5281/zenodo.20541629
---

# Summary

DotMatch is an open-source command-line and Python package for assigning short
DNA read windows to a known target table. Many sequencing assays already know
the exact short sequences they expect: a CRISPR guide library, an inline sample
barcode set, a feature-barcode whitelist, an amplicon primer prefix, or a panel
of designed assay tags. DotMatch focuses on this narrower problem rather than
general genome alignment. It extracts a configured read window, compares it with
the supplied target list under exact, Hamming, or Levenshtein edit-distance
semantics, and reports whether each read is uniquely assigned, ambiguous,
unmatched, or invalid.

The software provides native C routines, a command-line interface, Python
wrappers, workflow examples, and machine-checkable evidence gates. It is
distributed from GitHub with release citation metadata and a Zenodo archive
[@dotmatch_zenodo_017]. The current evidence supports use cases including
CRISPR guide counting, fixed-position barcode demultiplexing, feature-barcode
assignment, CRISPR guide-capture extraction, amplicon primer-start assignment,
adapter-prefix checks, and barcode panel design. DotMatch deliberately does not
claim to be a genome aligner, adapter trimmer, UMI/cell quantification pipeline,
production BCL converter, variant caller, or downstream CRISPR screen
statistics package.

# Statement of need

Short-DNA assignment is a recurring step in sequencing workflows, but it is
often implemented as a small script or as a side effect of a larger pipeline.
Researchers may need to answer a practical question before downstream analysis:
which known guide, barcode, primer, or whitelist sequence is present at this
fixed read position? The answer must be reproducible because ambiguous or unsafe
corrections can alter count matrices, sample splits, and quality-control
interpretation.

Existing tools address adjacent needs. MAGeCK provides statistical analysis for
CRISPR screens and includes guide-counting utilities [@li2014mageck]. Cutadapt
identifies and trims adapters and primers from reads [@martin2011cutadapt].
Edlib provides a fast exact edit-distance alignment library [@sosic2017edlib].
These tools are valuable, but none is specifically framed as an auditable
known-target assignment engine with explicit ambiguity states, fixed-window
FASTQ outputs, panel-design collision checks, and repository-level public-data
evidence gates. DotMatch fills that gap for laboratories and workflow authors
who need a small, deterministic assignment layer before statistical or
biological interpretation.

# State of the field

The bioinformatics ecosystem already contains mature tools for alignment,
adapter trimming, CRISPR screen statistics, and workflow distribution. DotMatch
is intentionally complementary to these tools. It uses edit-distance semantics
for short known targets, but it does not produce SAM/BAM, CIGAR strings, or
reference-index mapping output. It can count CRISPR guides, but downstream
phenotype modeling remains the responsibility of MAGeCK, BAGEL, drugZ, CERES,
or other screen-analysis software. It can classify fixed adapter-prefix
windows, but it does not trim reads or replace Cutadapt.

This narrow scope makes the software easier to validate. The repository records
public-data lanes, raw benchmark rows, generated reports, and gates for claims
that are allowed in documentation and release notes. Distribution work follows
life-science packaging norms, including Bioconda and container-oriented release
checks [@gruning2018bioconda].

# Software design

DotMatch represents assignment as a contract over a target table, a read window,
an edit metric, and an ambiguity policy. The default radius policy counts a read
only when exactly one target lies within the configured edit-distance radius.
If multiple targets are compatible, the read is reported as ambiguous instead
of being forced to an arbitrary best target. Invalid extraction windows and
unmatched reads are surfaced in QC outputs rather than silently discarded.

The implementation combines native C kernels with Python and CLI interfaces.
For supported A/C/G/T windows, indexed candidate generation avoids exhaustive
target scans while preserving the same assignment semantics as exhaustive
validation. Unsupported or out-of-scope cases fall back to semantics-preserving
paths. Output formats include counts, split FASTQ files, target audits,
top-unmatched tables, summary JSON, and HTML reports. The barcode-panel design
mode enumerates error spheres, reports nearest-neighbor risks, warns about
reverse-complement hazards, and writes machine-checkable collision summaries.

Correctness and claim boundaries are part of the project structure. Tests cover
native edit-distance routines, CLI workflows, Python wrappers, packaging
metadata, workflow fixtures, and public evidence gates. Documentation separates
supported statements from experimental or unsupported claims, so benchmark
language is tied to reproducible artifacts rather than broad performance
assertions.

# Research impact statement

DotMatch is intended for researchers building or auditing sequencing workflows
where known short sequences determine sample, guide, feature, primer, or panel
identity. The immediate impact is practical: reproducible assignment decisions,
explicit ambiguity accounting, and reusable methods text for the assignment
stage of an analysis. Public examples demonstrate CRISPR guide-counting and
barcode-style workloads using checked input data and comparator semantics.

The project is early-stage but citable and auditable. The Zenodo DOI allows
users to cite the exact software release, while the repository records the
commands and gates that justify current claims. Future impact should come from
standard workflow integrations, especially Bioconda/BioContainers propagation,
nf-core-style modules, Galaxy wrappers, and documented downstream uses.

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
