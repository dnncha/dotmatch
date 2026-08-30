# Galaxy CRISPR Counting Wrapper Example

This directory contains local example wrappers for running DotMatch from Galaxy.
The scoped IUC candidate is `dotmatch_crispr_count.xml`, which accepts one or
more CRISPR FASTQ datasets and is pinned to the publicly available Bioconda
package `dotmatch=0.2.2`. `dotmatch_assay_run.xml` remains a local AssaySpec
example with separate review scope.

Validate the XML with Planemo from the repository root:

```bash
planemo lint examples/workflows/galaxy/dotmatch_crispr_count.xml
planemo test --install_galaxy examples/workflows/galaxy/dotmatch_crispr_count.xml
planemo lint examples/workflows/galaxy/dotmatch_assay_run.xml
planemo test examples/workflows/galaxy/dotmatch_assay_run.xml
```

The wrapper targets the published DotMatch 0.2.2 Bioconda package and exposes a guide library, one or more FASTQ inputs, guide
offset, guide length, edit-distance threshold, metric, ambiguity policy, and an
optional one-base Levenshtein indel window. It writes a MAGeCK-compatible count
table, DotMatch summary JSON, and a `sample_qc.tsv` table suitable for MultiQC
custom content. The embedded Planemo test covers two samples plus unique,
ambiguous, unmatched, and invalid fixture reads.

The AssaySpec wrapper builds a reviewed `status = "ready"` TOML spec from
Galaxy-staged library and FASTQ inputs, then writes `assay_report.html`,
`assay_manifest.json`, `assay_manifest.summary.tsv`, `sample_qc.tsv`,
`crispr_qc.html`, `crispr_qc.json`, `crispr_qc.summary.tsv`, counts, and native
summary JSON. The assay report is the primary human-readable artifact;
`sample_qc.tsv`, `crispr_qc.summary.tsv`, and `assay_manifest.summary.tsv`
remain plain workflow-friendly tables.

Before adapting this for another Galaxy environment, pin an available DotMatch
package or container release and run Planemo lint/test against that environment.
