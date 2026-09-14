# From read assignment to gene-level robustness

Protocol AR-001, version 1.0, 6 September 2026.

## Status and disclosure

This is an exploratory public-data research protocol, not an externally registered report. The historical Yusa/Brunello discrepancy summaries were seen before this protocol. No newly executed pilot results were inspected before writing this version. Any changes must be dated and explained in AMENDMENTS.md. DotMatch is developed by the investigator; it is not the authority for biological truth.

Software baseline: dnncha/dotmatch commit 11d159fa1648365f2a4e96917b483c33aa5d9fe7 (0.5.0). Keep production code and historical benchmark artifacts unchanged. Research lives on its own branch.

## Question and estimands

Does fixed-window read assignment alter counts and downstream effect estimates in bulk pooled CRISPR screens, and do library sequence relationships explain those alterations? More assigned reads is not necessarily more correct reads. Policy disagreement is not an estimated error rate.

Stage A primary estimand: the fraction of evaluated read records whose status or uniquely assigned target differs between exact, unique-within-Hamming-one (radius_k1), and unique-nearest-within-Hamming-one (best_k1). Report every pair, the denominator, all four read states, per-guide and per-sample count deltas, and conservation checks. Secondary descriptive estimands: assignment fractions, duplicated reference sequences, distance-one guide neighbours, same-gene versus cross-gene ambiguity, and effects of extraction offsets. No statistical significance claims from reads treated as independent biological replicates.

## Inputs and sampling

Pilot inputs deliberately reproduce the existing examples: Yusa ERR376998 and ERR376999 (19 nt, zero-based start 23); Brunello PRJNA508200 plasmid SRR8297997 and RepA/B/C source run pairs SRR8297837+SRR8297836, SRR8297839+SRR8297838, SRR8297841+SRR8297840 (20 nt; extraction must be audited). Libraries are retrieved from the source URLs recorded in the baseline fetch scripts, never generated as a replacement for real references.

Stage A uses the first 102,000 complete records per sample: first 2,000 for extraction diagnostics and the next 100,000 for evaluation. These are archival prefixes, NOT random samples or whole-screen estimates. In a paired-run biological sample the prefix may come entirely from its first listed run. Report which source contributed. This stage is too shallow for credible genome-wide hit calling. Unavailable sources are recorded as unavailable, never silently replaced. A separate unchanged historical 100,000-record-prefix lane may be used to investigate prior benchmark discrepancies; label that lane explicitly and do not confuse it with the discovery-excluded lane.

## Extraction controls

Use the full supplied library, fixed strand and fixed windows for all policies. Yusa's documented start 23 is the prespecified primary window. For Brunello, tabulate exact matches at every possible offset from 0 through 40 on the 2,000 discovery reads only. Choose the modal exact-hit offset (lowest offset breaks a tie) for a controlled fixed-window pilot, record the full distribution, and flag mixed-offset assays. This controlled lane is NOT a recommended complete assay pipeline. Evaluate all prespecified policies on identical windows. Analyse a multi-offset or vector-flank extraction lane separately; never attribute extraction differences to assignment alone. Do not select a configuration by gene-level results.

## Independent validation and failure rules

Implement a separate Hamming-one reference matcher by enumerating each query's single-position A/C/G/T substitutions against a sequence-to-all-guide-IDs dictionary; do not import DotMatch assignment functions. For reference libraries containing non-ACGT symbols, fail this pilot rather than silently change alphabet semantics. Query N is literal, not a wildcard. Preserve duplicate reference sequences and distinguish duplicate IDs (invalid) from duplicate sequences (ambiguous). Validate the independent matcher against brute-force Hamming distance on exhaustive small alphabets and adversarial fixtures, including exact-plus-neighbour, equidistant ties, duplicate sequences, N and short windows. Brute-force-check deterministic stratified real windows against every reference target as an additional validation layer.

Compare DotMatch with the independent implementation at every guide count and read-state total for all policies, plus every changed-record call available in the read-change output. Require exact agreement under matching semantics. Count all input records exactly once; counts sum to unique states; state totals sum to evaluated records; per-guide deltas reconcile with unique-total deltas. No summary labelled complete when a required gate fails. Record failures and unresolved discrepancy classes, including possible DotMatch defects, without favouring either program.

## Provenance and resource boundaries

Record source URLs, ENA metadata, input sizes, record counts, local SHA-256, library content hashes, protocol and harness hashes, actual git commit, interpreter/platform, DotMatch version and executable/implementation hashes, comparator version/commit and exact command lines. For prefixes, record full-archive MD5 as NOT locally verified; derived prefix SHA-256 is a different check. Use bounded network timeouts, checked subprocesses and immutable result directories. Do not publish private data or raw public FASTQs in git. Deposit aggregates and code; retrieve source data at run time. No new paid infrastructure, journal submission or production release is authorized by this protocol.

## Historical discrepancy investigation

The existing summaries aggregate guide counts across samples; their 'ok' status means execution succeeded, not count identity. Recreate named sample axes before interpreting them. Check whether reported assigned counts exceed records, whether a read is counted at multiple offsets, reference duplicate handling, exact-first versus radius semantics, and stale/mismatched inputs. Reproduction failure is an outcome, not grounds for adjusting inputs to obtain the historical numbers. Do not treat historical or current disagreement as biological misassignment without independent evidence.

## Downstream stage B (not enabled by a shallow pilot)

Before calculating gene-level discoveries: obtain full checksummed runs, verify biological replicate/condition metadata against the source study, lock a cohort manifest and a new analysis-plan version, and pin MAGeCK or another established downstream method. Include at least three biological endpoint replicates where available; technical lanes are combined, never treated as biological replicates. A shared plasmid baseline is not three independent controls. Keep gene mappings, filtering and downstream settings identical between assignment policies. Do not synthesize p-values from pooled prefixes.

Primary downstream comparisons: changes in gene effect estimates and direction, ranks and hit sets at FDR 0.05; separately report threshold crossings near 0.05 and substantial effect changes (prespecified absolute log2-effect difference >=0.5). Estimate uncertainty across independent studies/biological replicates rather than across repeated pipeline runs. Use matched-depth controls, no invented gold standard, and independent validation for biological claims. Keep pilot-selected libraries separate from held-out studies. Define held-out inclusion criteria before inspecting their policy-dependent gene results.

## Publication gate

No claim of improved biological accuracy, a new mechanism, genome-wide prevalence, full historical replication, or a publication-ready manuscript from Stage A alone. A publishable report requires an explained mechanism, validation and adequate independent datasets, and must retain null/robustness findings. Uncertainty bounds over compatible assignments are a separate future method requiring a formal model, conservation constraints and calibration; they are not confidence intervals by default.
