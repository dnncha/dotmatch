# nf-core-style Module Candidate

This directory contains local nf-core-style module candidates for running
DotMatch from a Nextflow DSL2 workflow. `dotmatch_crispr_count` keeps the native
command path, while `dotmatch_assay_run` runs an AssaySpec and emits the assay
report and manifest summary. These examples have not been submitted to or
accepted by nf-core.

These examples demonstrate workflow integration. They are not benchmark results,
package-channel releases, or publication evidence.

The native module emits a MAGeCK-compatible count matrix, JSON summary,
`sample_qc` table for MultiQC custom content, and `versions.yml`. The AssaySpec
module emits `assay_report.html`, `assay_manifest.json`,
`assay_manifest.summary.tsv`, `sample_qc.tsv`, counts, summary, and
`versions.yml`. Its input tuple includes the AssaySpec plus the target/FASTQ
files referenced by that spec so workflow engines stage the required files into
the task work directory. Local nf-test candidates use the shared fixtures in
`examples/workflows/fixtures/`.

Before adapting this for nf-core or another workflow repository:

- pin a released DotMatch package or container;
- run and adapt the nf-test candidate against the target repository conventions;
- add CI linting against the target workflow repository conventions;
- keep module help text aligned with the documented evidence boundaries.
