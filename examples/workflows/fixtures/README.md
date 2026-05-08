# Workflow Test Fixtures

These tiny CRISPR guide-counting fixtures are shared by local workflow examples
and upstream-submission test stubs. They are intentionally small enough for
Planemo, nf-test, and CI smoke checks.

`sample_a.fastq` covers the core DotMatch outcomes:

- `unique`: `ACGT` exactly matches `guide_a`.
- `ambiguous`: `ACGG` is one Hamming edit from both `guide_a` and `guide_b`.
- `unmatched`: `CCCC` is outside the configured edit threshold.
- `invalid`: `AC` is shorter than the configured guide window.

`sample_b.fastq` adds a second sample with one exact `guide_c` assignment and
one unmatched read so MAGeCK-style multi-sample output is exercised.
