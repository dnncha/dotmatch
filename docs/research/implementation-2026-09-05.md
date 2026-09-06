# Assignment-sensitivity implementation and review

This iteration adds a one-pass policy comparison around the existing native
Hamming index, fixes Python target-table parsing and hardens benchmark count
interpretation. The matching kernel and its default semantics are unchanged.

## Engineering measurement

`python scripts/bench_policy_sensitivity.py` compares the fused path with three
separate native queries, using a prebuilt warmed index. Seed 63104; 4,000 random
20-base targets; 100,000 synthetic reads; batches of 4,096; five alternating
repetitions. All statuses and unique target IDs agree in every repetition.

In the recorded local run, median time was 1.2029 seconds for three separate
queries and 0.6345 seconds for the fused comparison, a ratio of 1.90. The raw
measurements are in `benchmarks/raw/policy_sensitivity_local.json`.

This is a same-process synthetic engineering measurement of the policy-comparison
step, including Python result projection. It excludes index construction and
FASTQ file I/O. It is not a competitor comparison, biological accuracy result,
or a claim of 1.90× end-to-end speedup. Re-run before extending the claim to
other hardware, workloads or release versions.

## Interface and reporting choices

The website replaces the decorative workflow picture with a native-generated
interactive read example. A worked comparison shows an equal unique-read total
with different guide counts. The example is intentionally synthetic and its
source inputs, native results and independent candidate counts are tested.

The README is shorter and starts with the scientist's task, installation and
first run. The six published agent tools remain documented. The new sensitivity
command is explicitly labelled source-only until the next package release.

Output schemas retain statuses separately from biological claims. Malformed count
tables are rejected, not rounded or silently zeroed. Aggregate guide totals are
labelled separately from named-sample matrix identity. Historical benchmark files
are not retroactively regenerated or made to appear equivalent.

## Release gates

Run the complete native/CLI/Python suites and existing evidence gates. Check the
four static routes, page-specific canonical URLs and sitemap. Browser tests cover
390, 768 and 1440px widths, actual hydrated example controls, mobile navigation,
library validation/export, stale-result clearing and the absence of network
requests during library interactions. Inspect the resulting screenshots before
merging. Verify deployment separately from source/PR validation.
