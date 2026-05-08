# nf-core-style Module Candidate

This directory contains an nf-core-style module candidate for running
`dotmatch crispr-count` from a Nextflow DSL2 workflow. It is a local integration
starting point, not an upstream nf-core module and not external adoption.

The candidate follows the same evidence boundary as the other workflow examples:
it demonstrates how to wrap DotMatch in a workflow manager, but it is not a
benchmark result, package-channel release, or publication claim.

The module emits a MAGeCK-compatible count matrix, JSON summary, `sample_qc`
table for MultiQC custom content, and `versions.yml`. A local nf-test candidate
in `modules/local/dotmatch/crispr_count/tests/main.nf.test` uses the shared
fixtures in `examples/workflows/fixtures/`.

Before proposing this to nf-core or another external workflow repository:

- pin a released DotMatch package or container;
- run and adapt the nf-test candidate against the target repository conventions;
- add CI linting against the target workflow repository conventions;
- keep module help text aligned with `docs/scientific-claims.md`.
