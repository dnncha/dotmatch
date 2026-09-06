# AR004: shared-read evidence and calibration

Version 1, 6 September 2026. Publicly versioned before executing new AR004 numerical analyses. Not a registered biological study. Prior AR001–AR003 results and the Barcas prior-art overlap have been inspected. No claim that positional ambiguity or grouped assignment is new.

## Questions

A. Does record reuse imply measurable technical dependence that ordinary guide count tables cannot encode? Can a read-class sufficient statistic restore a correctly calibrated technical null without asserting biological truth?

B. Do fixed-window versus multi-offset counting results differ under the actual MAGeCK executable when the additional original ESC2 archive is included? This is a within-study descriptive sensitivity extension, not a held-out external study or an estimate of which hit list is true.

## Locked old inputs

Use the complete verified replay artifacts from run 34030217143 for ERR376998, ERR376999, SRR8297997 and original raw references. ZIP and internal file hashes must match provenance. Reconstruct singleton record classes from upstream counts minus full multiple-target class memberships only after verifying repeated-same-target events=0. Account separately for unmatched records. All guide counts and event/record budgets must reconcile. Candidate classes from AR003 are not interchangeable with actually counted upstream events.

## A: empirical provenance and conditional technical null

For each gene and sample compute C=sum_r m_rg, U=sum_r 1[m_rg>0], Q=sum_r m_rg^2, where m_rg is the number of accepted guide rows of that gene receiving record r. Under independent fair coin assignment of each original sequencing record to two artificial groups, the event-count difference D_g has exact conditional mean zero and variance Q, whereas an event-independent calculation uses C. Report Q/C, C^2/Q, U, complete guide-pair overlap and source-class budgets, including unaffected genes. C^2/Q is an effective record-weight quantity under this model, not biological sample size.

Generate 2,000 null label splits per archive with seed 20260906, in reproducible bounded batches. The same record receives one random label shared by all its events. Primary eligibility: C>=100, defined solely from this single archive before label splits. Compare two-sided normal-tail tests using sqrt(C) versus sqrt(Q), threshold 0.05, recording all per-gene rates and BH<=0.05 discovery counts per split. These tests diagnose the defined technical null and are explicitly NOT MAGeCK, biological replicates, sequencer error estimates, or clinical p-values. Monte Carlo repetitions are not new biological evidence. Inspect discrete/low-count effects and compare exact weighted-binomial tails for a prespecified structural extreme (largest Q/C among eligible genes, ties lexical), plus an unaffected coverage-matched control. The selection uses provenance, not favourable null outcomes.

Primary summaries: eligible populations, fraction with Q/C>=1.25, >=1.5, >=2, global and structurally affected false-positive fractions at nominal 0.05, their uncertainty across random splits, complete per-gene values. Show both the worst structural example and the entire population. Compare the same count totals under a synthetic event-independent origin assignment to establish count-matrix non-identifiability of dependence; label this as a mathematical counterexample, not a second biological sample. Include balanced synthetic controls and enumeration tests. No new theorem claim.

## B: actual replicated screen sensitivity

Acquire ERR377000 completely from ENA, verify metadata accession, archive MD5, bytes, record count and SHA-256; no prefix or synthetic fallback. The NIH MAGeCK-VISPR example identifies ERR376998 as plasmid, ERR376999 as ESC1 and ERR377000 as ESC2; check original metadata and document any ambiguity. Keep original 87,437-entry reference and 19-base start 23 fixed. The 87,897 designed-guide figure in the source publication is a separate unresolved library-history reconciliation, not evidence of a wrong reference.

Pin actual MAGeCK 0.5.9.5 to the source archive and SHA-256 in the retrieved Bioconda recipe, compile it and retain source/binary/environment hashes and help. Pin DotMatch 0.5.0 original engine, unchanged guide-counter v0.1.3 original comparator where run. Compare plasmid versus both ESC replicates using the same MAGeCK options, library and sample roles: median normalization, FDR=0.05, gene LFC median. Keep fixed exact/best/radius and native multi-offset-event arms distinct. Report rank/effect/hit agreement and all changed results, not just new significant genes. Code-level runtime adaptations must be recorded, not hidden. If the actual executable is unavailable, do not describe an emulation as MAGeCK.

The additional archive is from the same previously inspected study, not a new independent validation cohort. Cross-policy hit changes do not determine biological truth. Biological accuracy and generalization still require external replicated screens, context-aware comparators (including Barcas/ReCo/CRISPR-Correct as applicable), and orthogonal validation.

## Reporting and reproducibility

Preserve prior results unchanged. Report failed executions and protocol deviations. Record which parts ran on raw reads, which reused archived complete sufficient statistics, and which are simulated or algebraic. Independent arithmetic/enumeration checks are not biological validation. Ship code, complete tables, checksums and a research-only PR; do not merge production defaults, contact authors or submit a manuscript automatically.
