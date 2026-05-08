# Galaxy CRISPR Counting Wrapper Example

This directory contains an example wrapper for running `dotmatch crispr-count`
from Galaxy. It is an example wrapper for local evaluation and wrapper
development, not a ToolShed release.

Validate the XML with Planemo from the repository root:

```bash
planemo lint examples/workflows/galaxy/dotmatch_crispr_count.xml
planemo test examples/workflows/galaxy/dotmatch_crispr_count.xml
```

The wrapper exposes a two-sample CRISPR guide-counting surface: guide library,
two FASTQ inputs, sample labels, guide offset, guide length, edit-distance
threshold, metric, ambiguity policy, and optional one-base Levenshtein indel
window. It writes a MAGeCK-compatible count table, DotMatch summary JSON, and a
`sample_qc.tsv` table suitable for MultiQC custom content. The embedded Planemo
test uses `test-data/` fixtures copied from `examples/workflows/fixtures/`.

Before publishing to a Galaxy ToolShed, pin a released DotMatch package from
Bioconda or a container, run Planemo lint/test against the target Galaxy
environment, and keep the help text aligned with `docs/scientific-claims.md`.
