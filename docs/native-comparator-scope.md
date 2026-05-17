# Native Comparator Scope

DotMatch currently has one native alignment-library comparator: Edlib exhaustive global edit-distance assignment. The generated native report records Edlib through `EDLIB_MODE_NW`, `EDLIB_TASK_DISTANCE`, fixed threshold `k`, and zero assignment mismatches before speedups are reported.

SeqAn and Parasail are not part of the checked README, website, or release-note
performance comparison set yet. Before either name belongs in those comparisons,
the repository needs all of the following:

- equivalent global edit-distance or documented semi-global scoring semantics for the exact workload being claimed;
- fixed threshold `k` and identical assignment policy for unique, ambiguous, no-match, and invalid reads;
- native dependency name, version, build flags, and platform in the raw artifact or generated report;
- raw CSV rows under `benchmarks/raw/` plus a generated report under `docs/benchmarks/`;
- zero assignment mismatches against DotMatch and the selected comparator;
- a gate script that fails when only scaffold, smoke, or unmatched-scoring rows are present.

Until that evidence exists, the supported native comparison wording is limited to Edlib exhaustive global edit-distance assignment scans plus the exact-hash and BK-tree baselines recorded in `docs/benchmarks/native/README.md`.
