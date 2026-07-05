# Evidence Packet v1

An evidence packet is a reviewer- and customer-ready bundle that explains what
was run, which deterministic assignment rules were used, what evidence was
produced, and where ambiguity or failure modes remain visible. This document
defines the public v1 shape so the open-source engine and DotMatch Pro can use
the same reliability language.

Evidence packets should include:

- run identity: DotMatch version, command, input manifest, timestamps, host
  platform, and workflow context;
- assay specification: target table, barcode or guide metadata, fixed-window
  coordinates, distance metric, correction radius, and quality filters;
- deterministic outcomes: counts for `unique`, `ambiguous`, `none`, and
  `invalid` reads, with no silently forced assignments;
- provenance: hashes or immutable references for public or synthetic inputs,
  generated TSV/JSON/HTML outputs, and report build commands;
- safety checks: duplicate targets, ambiguous near-neighbors, invalid windows,
  unsafe correction radii, and panel-collision findings when applicable;
- validated scope: what the packet supports, what it does not support, and which
  public evidence gate or benchmark report is relevant;
- review artifacts: summary tables, top-unmatched tables, panel or target
  audits, QC summaries, and a human-readable report.

Evidence packets must not include private customer FASTQ/BAM/BCL files or other
raw biological data in the public repository. Public examples should use
synthetic, minimized, or public reference data with source and scope clearly
documented.

DotMatch Pro may add signing, team review state, private assay registry links,
access controls, and audit-pack export around this public shape. Those additions
should preserve the open engine's assignment semantics and make provenance
clearer, not less inspectable.
