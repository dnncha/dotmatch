# AR004: replicated Brunello dropout screen

Prospective analysis amendment, 6 September 2026. The complete AR003 results, prior Brunello plasmid counts, primary Sanson et al. paper and study run metadata were inspected before this amendment. No new cellular FASTQs or gene outcomes have been inspected. This is a versioned research protocol, not an external preregistration.

## Dataset and nonindependence

Primary source: Sanson et al. 2018, Nature Communications, doi:10.1038/s41467-018-07901-8, PRJNA508200. The paper describes A375 biological triplicates for Brunello with modified tracrRNA, 3-week dropout and a matched plasmid reference.

Use ALL single-end run archives belonging to the following original sample aliases:
- Brunello_mod_tracr_pDNA: SRR8297997.
- Brunello_mod_tracr_RepA_Dropout_A375: SRR8297836 and SRR8297837.
- Brunello_mod_tracr_RepB_Dropout_A375: SRR8297838 and SRR8297839.
- Brunello_mod_tracr_RepC_Dropout_A375: SRR8297840 and SRR8297841.

Verify original sample accessions agree within each alias. Combine technical runs into A/B/C before statistics; never treat six files as six biological replicates. The plasmid baseline is shared and not biologically replicated. This is an independent replicated cellular cohort relative to AR003, but not an independent study or library: the plasmid and reference were already used. Include all runs rather than selecting the smaller lane or a favourable result. No gene-discovery claim follows from a policy-specific p-value.

Reference: the original corrected 77,441-target Brunello table with SHA-256 0d2906187829ea9f736de94a47369bd94d42cde5f348fea9d12a385625cc2ca1. Original annotations retained. Verify archive MD5, byte counts and full record counts; lock SHA-256. No prefix substitution, synthetic fallback or reconstructed FASTQs. Use raw URLs from the previously committed complete ENA metadata. Expected aggregate read count is computed from these frozen rows, not guessed. Upper transport budget 8 GB.

## Assignment arms and principal contrast

Use zero-based starts 21 through 30 inclusive and 20-base windows, frozen from the previously inspected Brunello plasmid discovery, identically for every run. Complete all declared windows or mark a record invalid. No indels, reverse-complement search, Phred weighting, or implicit choice of the best-looking offset.

Enumerate complete Hamming-radius-one candidates at every position, retaining duplicate target identities. The optimized research implementation must reproduce the AR003 Python joint decoder and pinned DotMatch 0.5.0 on validation inputs before public execution.

Five guide-level arms, with separate original-gene resolution for the three joint arms:
1. event_exact: increment one unique exact target at each accepted position (can count one read multiple times).
2. event_best: increment one uniquely nearest target within one mismatch at each accepted position (can count one read multiple times).
3. joint_exact: retain all exact position/target explanations, count a guide once only if unique globally.
4. joint_best: retain globally minimum-distance explanations within one mismatch, count a guide once only if unique globally.
5. joint_radius: retain all radius-one explanations, count a guide once only if unique globally.

Primary matched-domain contrast: event_best versus joint_best. Secondary: event_exact versus joint_exact, joint_exact versus joint_best, joint_best versus joint_radius. Unlike the prior fixed-window comparison, these arms share precisely the same permitted positions. Record event and distinct matched-read totals, repeated same-target versus additional distinct-target events, same-gene versus cross-gene ambiguity, and complete affected guide pairs. event_best is an explicit rule implemented here, not silently labelled actual guide-counter execution.

Native-workflow comparator: attempt the published ReCo workflow using its required Cutadapt 2.8 and Bowtie2 2.3.0 and record exact source/dependencies. Its automatic context/extraction differs and must remain a separate arm; do not call mapping yield accuracy. Preserve failures and never replace its results with a hand-written imitation. A successful current version is not automatically equivalent to the paper implementation.

## Downstream statistics

MAGeCK 0.5.9.5 RRA, pinned source/environment. Per-library target identities and annotations are fixed. Baseline-only common inclusion: joint_exact plasmid >=30 reads. Use exactly the same retained target IDs for all five arms and the same A/B/C versus plasmid contrast. Primary normalization median; remove-zero none; gene-lfc-method median; negative-selection FDR threshold 0.05. Do not use non-targeting guides as the sole normalization/null population for this growth screen, because the known cutting effect makes them unlike targeted nonessential genes. Record omitted controls/NA annotations explicitly.

Primary endpoints: all gene LFC differences; median and maximum absolute difference; genes with |delta LFC|>=0.5; Spearman rank agreement; FDR<=0.05 call overlap (both, only one, neither); thresholds 0.01 and 0.10 as sensitivity; distribution of FDR proximity for discordant calls. Retain all genes/outliers. Report at least two analyses of an identical baseline table to expose any stochastic RRA variation; where the software supports a seed, set and record it. Call-set changes are computational sensitivity, not newly validated biology or FDR calibration evidence.

Descriptive replication endpoints: per-replicate input-read-normalized gene-count log2 ratios with 0.5 pseudocount; direction consistency across A/B/C; same-guide/gene changes observed in all three replicates. Do not count guides or reads as biological replicates. Gene-unique tables and conditional bounds are diagnostics, not fake single-guide input to a guide-based statistical model.

Secondary equal-depth check: deterministically hash-thin each technical run to the minimum total biological-replicate exposure (use one sample-specific rate per A/B/C; plasmid unchanged), preserving the same retained records across assignment arms. Hash selection based only on run accession and record ordinal, not guide identity or outcome. Full-file counts remain primary. If this secondary is not executed, say so explicitly.

## Validation and publication gates

Before new public results, validate optimized all-position candidates and assignments against AR003 and independent exhaustive enumeration on constructed panels, duplicates, literal N, missing windows, order/cache invariance and error-model violations. Recount available AR003 Brunello plasmid and require complete guide/joint gene counts to match. Validate randomly selected public windows against native DotMatch and a full-reference exhaustive oracle. Exact count-budget identities, all input joins, technical-run aggregation, integer counts, file hashes and failure-visible completion markers are mandatory.

Preserve a negative finding. This stage can establish downstream sensitivity in one replicated screen. It does not establish general superiority, false-discovery inflation, independent-study replication, true read origin or a new gene mechanism. Do not change production defaults, submit a manuscript, or contact third parties automatically. Research branch/PR, executed report, methods draft and complete derived evidence are authorized outputs.
