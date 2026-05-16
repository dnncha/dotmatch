# Galaxy CRISPR Counting Wrapper Example

This directory contains local example wrappers for running DotMatch from Galaxy.
`dotmatch_crispr_count.xml` keeps the native command interface, and
`dotmatch_assay_run.xml` runs an AssaySpec and exposes the assay report and
manifest summary. These examples have not been published to a Galaxy ToolShed.

Validate the XML with Planemo from the repository root:

```bash
planemo lint examples/workflows/galaxy/dotmatch_crispr_count.xml
planemo test examples/workflows/galaxy/dotmatch_crispr_count.xml
planemo lint examples/workflows/galaxy/dotmatch_assay_run.xml
planemo test examples/workflows/galaxy/dotmatch_assay_run.xml
```

The wrapper exposes a two-sample CRISPR guide-counting surface: guide library,
two FASTQ inputs, sample labels, guide offset, guide length, edit-distance
threshold, metric, ambiguity policy, and optional one-base Levenshtein indel
window. It writes a MAGeCK-compatible count table, DotMatch summary JSON, and a
`sample_qc.tsv` table suitable for MultiQC custom content. The embedded Planemo
test uses `test-data/` fixtures copied from `examples/workflows/fixtures/`.

The AssaySpec wrapper builds a reviewed `status = "ready"` TOML spec from
Galaxy-staged library and FASTQ inputs, then writes `assay_report.html`,
`assay_manifest.json`, `assay_manifest.summary.tsv`, `sample_qc.tsv`, counts,
and native summary JSON. The report is the primary human-readable artifact;
`sample_qc.tsv` and `assay_manifest.summary.tsv` remain plain workflow-friendly
tables.

Before adapting this for a Galaxy ToolShed, pin an available DotMatch package or
container release, run Planemo lint/test against the target Galaxy environment,
and keep the help text aligned with the documented evidence boundaries.
