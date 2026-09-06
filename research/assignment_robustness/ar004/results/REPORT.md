# Replicated-screen assignment sensitivity: AR004

**Executed technical study, 6 September 2026. Not peer reviewed.**

Seven complete public archives supplied **246,950,411 sequencing records**: one matched modified-tracr plasmid baseline and three A375 dropout biological replicates, each sequenced in two runs. Technical runs were summed within biological samples before statistics.

## Primary matched-position result

At MAGeCK negative-selection FDR <=0.05, per-position best-distance event counting called **1,124** gene annotations and joint best-distance counting called **1,119**. **1,118** were common, **6** were event-only and **1** were joint-only. These are method-dependent calls, not verified true or false discoveries.

Among 19,112 common reported genes, **20** had absolute MAGeCK log2-effect differences >=0.5. Negative-selection rank correlation was **0.999383**. Both counting arms used the same ten allowed positions, one-mismatch model, guide reference, common baseline-filtered guide subset and downstream settings.

## All prespecified contrasts

| Scope | Left | Right | Common genes | Effect changes >=0.5 | Only left at 0.05 | Only right at 0.05 |
|---|---|---|---:|---:|---:|---:|
| full | event_best | joint_best | 19,112 | 20 | 6 | 1 |
| full | event_exact | joint_exact | 19,112 | 3 | 0 | 0 |
| full | joint_exact | joint_best | 19,112 | 355 | 65 | 41 |
| full | joint_best | joint_radius | 19,112 | 99 | 16 | 25 |
| thin | event_best | joint_best | 19,112 | 19 | 4 | 3 |
| thin | event_exact | joint_exact | 19,112 | 2 | 0 | 2 |
| thin | joint_exact | joint_best | 19,112 | 352 | 65 | 38 |
| thin | joint_best | joint_radius | 19,112 | 99 | 17 | 21 |
| identical_input_repeat | event_best | event_best | 19,112 | 0 | 0 | 0 |
| identical_input_repeat | joint_best | joint_best | 19,112 | 0 | 0 | 0 |

The thinned analysis uses deterministic outcome-independent Bernoulli selection targeting the minimum expected cellular-replicate exposure. Realized counts are recorded and need not be exactly equal. Identical selected records are used across policies. Repeated-input comparisons measure variation from running the same statistical program again; they are not independent experiments.

## Reproducibility and validation

The baseline-only common filter retained **74,293 guides** spanning **19,112 original gene labels**. Eligibility requires >=30 joint-exact plasmid counts; non-targeting and absent-annotation groups are explicitly excluded from gene inference. All reference rows remain in original and aggregated count tables, with every exclusion in `eligibility.tsv`.

Input archive byte lengths, MD5, SHA-256 and complete record totals were verified. The optimized C++ counter first reproduced every AR003 Brunello plasmid joint guide count and gene bound. Every new archive was checked with 200 seeded records against all reference targets at every allowed position and against pinned DotMatch 0.5.0. Complete count-state budgets and technical-run/annotation aggregation were verified, including independent SQL sums. Sampled all-target validation is not represented as exhaustive validation of all 247 million records.

MAGeCK 0.5.9.5 was built from the original source. Median normalization, no zero-count removal, median gene log-fold changes and the same A/B/C versus shared plasmid contrast were used for all arms. Full commands, executable/source hashes, environment, original gene outputs, normalized counts and intermediate files are retained. The primary FDR threshold is 0.05; 0.01/0.10 are sensitivity thresholds, not separate experiments.

## Scientific interpretation

This stage answers a narrower and stronger question than a mapping-rate comparison: do read-accounting choices change statistical conclusions in a replicated screen when the search domain is held fixed? The answers above quantify that sensitivity. They do **not** identify which changed calls are biologically correct, establish false-discovery calibration, or demonstrate general superiority.

The cellular outcomes were unseen before the AR004 protocol, but this remains the same Sanson study and Brunello library as the prior plasmid audit, not independent-study replication. The three biological replicates share one plasmid baseline. Gene annotations are retained as supplied. No copy-number correction, independent gene-function validation or experimentally known molecular origin is added. Event counting is the explicitly implemented per-position rule, not an undeclared emulation presented as actual guide-counter execution. ReCo is a separate native-workflow comparison with different extraction and must be reported separately.

All significant/non-significant and discrepant rows are retained. Genes near an FDR boundary should not be called discoveries merely because one configuration moves them across it. Stable effects alongside changed p-values require inspection of guide evidence, normalization and statistical variability rather than an automatic biological narrative.

## Sources and status

Primary study: Sanson et al., *Nature Communications* (2018), doi:10.1038/s41467-018-07901-8. Original run metadata: ENA PRJNA508200, archived before outcomes. MAGeCK source/documentation: https://sourceforge.net/projects/mageck/ and https://sourceforge.net/p/mageck/wiki/usage/. Frozen AR004 protocol commit: `dc1971a9a9b3938dea17e6961ff16a515a52cd96`. Complete counting commit: `862bdc1e0c026b2ce7475870af1a4a2171c98cfb`, workflow `34034542345`.

No production defaults, manuscript submission, author outreach or new biological mechanism is claimed.
