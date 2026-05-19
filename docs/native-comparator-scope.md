# Alignment Library Comparison Scope

DotMatch currently has one checked alignment-library comparison: Edlib full
global edit-distance assignment. The generated Edlib report records
`EDLIB_MODE_NW`, `EDLIB_TASK_DISTANCE`, fixed threshold `k`, and zero assignment
mismatches before speedups are reported.

SeqAn and Parasail are not part of the checked README, website, or release-note
performance comparison set yet. Before either name belongs in those comparisons,
the repository needs all of the following:

- equivalent global edit-distance or documented semi-global scoring rules for the exact workload being claimed;
- fixed threshold `k` and identical assignment policy for unique, ambiguous, no-match, and invalid reads;
- dependency name, version, build flags, and platform in the raw file or generated report;
- raw CSV rows under `benchmarks/raw/` plus a generated report under `docs/benchmarks/`;
- zero assignment mismatches against DotMatch and the selected comparison tool;
- a check script that fails when only scaffolding, small fixtures, or unmatched-scoring rows are present.

Until that evidence exists, comparison evidence is limited to Edlib full
global edit-distance assignment scans plus the exact-hash and BK-tree
comparisons recorded in `docs/benchmarks/native/README.md`.
