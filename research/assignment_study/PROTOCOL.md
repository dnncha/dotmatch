# Frozen exploratory study: multi-offset accounting and CRISPR abundance contrasts

Frozen before new full-run outcomes. Not a prospective preregistration: older small-prefix count discrepancies were already known. Question: does per-offset counting rather than per-read counting alter guide and gene abundance contrasts? Fixed-window mismatch correction is a separate contrast. Computational agreement is not molecular-origin truth.

## Selected data
Use complete files ERR376998, ERR376999, ERR377000 (technical set A) and ERR377001, ERR377002, ERR377003 (set B). Each triplet comprises plasmid SAMEA2188756, mESC library 1 SAMEA2188757, and mESC library 2 SAMEA2188758. B repeats the SAME three sample identities; it is not additional biological replication. Primary ENA sample metadata determined selection. All other PRJEB4038 runs are outside this bounded study. No toxin-selection or drug-resistance contrast is analyzed. Separate library names do not establish independently randomized biological replicates.

The unmodified MAGeCK-distributed Yusa library contains 87,437 unique 19nt targets; SHA256 d41a4122a46fdedb47deca61081dde9037d461a1fce1e6f9744522e21005ceba. No target exclusion or deduplication is allowed. Full FASTQ bytes and ENA MD5 and observed record counts must agree.

## Arms
Fixed start 23, length 19: exact unique, Hamming radius-one unique, unique nearest within one mismatch. Run actual pinned DotMatch native CPU and independently enumerate every ACGT radius-one neighborhood in C++. A single non-ACGT read symbol is a literal mismatch (not a wildcard); the independent oracle checks its four substitutions against exact targets.

In a separate experiment, recreate pinned guide-counter offsets from first 100,000 reads, retaining each offset with >=0.0025 of total matched WINDOWS (not reads). Exact and best-distance modes learn separate offset sets. Require equality with actual upstream per-offset counts. For each mode use the same offsets and successful window hits but count each distinct guide only once per read. This controlled arm may still count different guides from one read; quantify conflicts. Joint nearest unique across offsets is a separately labelled exploratory arm, not the same intervention.

Retain full and calibration-excluded counts, every multi-hit record ordinal and read ID, offsets/guide IDs, canonical/alternate associations, complete state denominators. Canonical window 23 is an analysis reference, not origin truth. Separate same-guide, same-gene cross-guide and cross-gene multiplicity.

## Gates
Synthetic seeded fixture: all eight count policies independently recomputed using pure-Python exhaustive Hamming distances. For real data compare indexed calls against full-library brute force on systematic record ordinals (modulo 10007, capped at 1000 total checks) and first 20 held-out examples per outcome stratum. Require full count-vector equality for three DotMatch and two guide-counter policies on all six files: 30 comparisons. Preserve failures and block interpretation on mismatches.

## Downstream analysis
Primary contrasts: distinct-guide minus per-offset effects for exact and best-distance; ESC1/plasmid and ESC2/plasmid. Secondary: fixed-best minus fixed-exact, fixed-radius minus fixed-best. Analyze technical sets A and B separately, then sum counts by the three sample identities.

Common median-ratio size factors from fixed-exact guides positive in all samples; pseudocount 1 in normalized-count space. Policy-specific median-ratio normalization is a sensitivity analysis. Eligibility: >=30 fixed-exact plasmid counts per technical set, >=60 when summed. Gene median log2 effects require >=3 eligible guides. Descriptive flag abs(delta log2 effect)>=0.5; stronger candidates need >=2 guide changes in the same direction in each ESC library and consistent sign across technical sets. These are not p-values or proof of dependency.

Report nulls, all denominators, medians, quantiles and largest effects. No gene false-discovery claim or new dependency is justified by these library abundance contrasts alone. No external cohort, biological ground truth, or clinical validation is included. Whole-study hit-calling and causal validation are outside this bounded measurement study.

## Provenance
DotMatch base 5c3934d648ff55928db5a9219d7703845941d28e; guide-counter f24de175282064773328e35baf03373e7e3b895d. Record source/binary/input hashes, commands, failed comparisons, and complete outputs. Do not redistribute full raw FASTQs. Conflict: study owner authored DotMatch. No publication acceptance or novelty is assumed; no journal/preprint submission is authorized by execution alone.
