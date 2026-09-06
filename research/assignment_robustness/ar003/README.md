# AR003: retain gene information without inventing guide certainty

**Executed research results, 6 September 2026. Not peer reviewed. No production behaviour change.**

The joint guide/position decoder has now run on all **30,215,791 original sequencing records** across Yusa ERR376998, Yusa ERR376999 and Brunello SRR8297997. All three full-archive jobs passed. The analysis independently reconstructed the published counts and conditional bounds from the complete candidate classes, and all ten canonical report files reproduced byte-for-byte from a clean source checkout.

[Read the numerical report](results/REPORT.md) · [Validation and limitations](results/VALIDATION.md) · [Frozen protocol](PROTOCOL.md) · [Research PR #108](https://github.com/dnncha/dotmatch/pull/108)

## What we can say

Under global best-distance matching, **334,089 Yusa records** are identifiable at gene-annotation resolution but not at individual-guide resolution: 169,983 plasmid and 164,106 cellular records. Under radius-one uniqueness the corresponding total is 1,392,988. These are separate policies applied to the same records; their counts must not be added. Gene and guide tables likewise overlap and must not be added.

The gain is resolution-specific information, not newly created reads, deduplicated molecules, independent biological evidence, or proof of improved accuracy. Multiple compatible guides belonging to one gene can leave that gene identifiable even when guide identity is unresolved.

Large gene-representation changes remained uncommon. Compared with the fixed-window best-distance baseline, the joint gene-unique result changed by at least 0.5 log2 units for 12 of 18,721 eligible annotation groups under best-distance matching, 70 under exact matching, and 48 under radius-one uniqueness. No strong direction reversal occurred under the declared descriptive definition. This is a comparison of measurement configurations on an unreplicated plasmid/cellular pair, not a gene-discovery experiment.

The constructed error-free Yusa controls retained correct gene identity for all 87,437 records under exact and global-best policies, while only 85,824 records identified a unique guide. The controls deliberately include the correct source in the reference and use known assay contexts. They do not establish real-data accuracy, calibrate uncertainty, or cover errors outside the declared model.

## Why the eligible denominator differs from the preceding report

The earlier standalone continuation at commit `e80619f2c963cdcbf12b998e3307a9a83388974f` used exact-policy plasmid gene counts **>=100 and at least two reference guides**. Its primary eligible population was **18,564** annotation groups.

AR003's separate prospective protocol, frozen at `de0002a37259cf27000015f43dd335343cd01234` before its new full-data analysis, uses exact-policy plasmid gene counts **>=50 and at least two guides with >=10 plasmid counts each**. That produces **18,721** annotation groups. Neither eligibility definition uses the cellular outcome. AR003 also reports baseline thresholds 20 and 100 while retaining its two-supported-guide requirement; these are not identical to the preceding phase's rule.

The configurations being compared also changed: the preceding primary comparisons used fixed-window policies and multi-offset event counts; AR003 compares joint gene-unique resolution against a fixed-window best-distance baseline. Therefore the two reports' headline numbers are **not a matched-cohort performance comparison**, and smaller differences in AR003 must not be promoted as proof of improved accuracy. The original results and protocols remain unchanged. This explanatory note was written after the new results were inspected; it does not modify endpoints or recompute a preferred cohort.

## Reproduce or review

- Full raw-read audit: [run 34032669590](https://github.com/dnncha/dotmatch/actions/runs/34032669590), exact code commit `68a406b41805b808e6428927c2b8ccd2a2e55f9f`, pinned native DotMatch commit `11d159fa1648365f2a4e96917b483c33aa5d9fe7`.
- Completed independent report/replay and packet: [run 34033408569](https://github.com/dnncha/dotmatch/actions/runs/34033408569), report source commit `41d7036b483d9cc163ae3643db1e9658d192bac9`.
- [Complete review packet](https://github.com/dnncha/dotmatch/actions/runs/34033408569/artifacts/9989366169): 69,324,993 bytes; ZIP SHA-256 `338fc2ece0c6e9847d60ebf2bda172b968546d858dde0d71fe4fe585ededb293`. This immutable packet precedes the present explanatory README; its numerical evidence and source are unchanged. GitHub artifact retention expires 5 December 2026, so retain a local copy.

The packet contains the full source, all derived candidate/count evidence, prior replay evidence, complete numerical tables, validation logs and a standard-library-only report reproducer. Original public FASTQs are not redistributed. Their official archive URLs, SHA-256, MD5, byte sizes and record counts are retained. Clean-report replay uses the same verified audit evidence; it is not a second raw-FASTQ execution.

The first report-build attempt failed on a formatting syntax error before producing results. That failure is retained in run 34033081891. The fix changed report rendering, not the completed raw decoder or frozen numerical endpoints. Unfinished local runs from the failed local runtime are excluded from evidence.

## What would make the next paper strong

The result worth testing next is not that 'published screens are broken.' It is whether apparently separate guide support recycles the same sequencing records, and whether that changes conclusions in independently replicated screens. A stable gene point estimate does not settle that question, and this study has not demonstrated false-discovery inflation.

A stronger methods paper needs outcome-blind independent replicated screens, appropriate assay-aware comparators, and known-origin or orthogonal accuracy/calibration evidence. The standalone Brunello plasmid is not a biological replicate; matched scaffold/library/baseline identities must be confirmed. Gene-grouping ideas and elementary candidate bounds are established concepts, not claimed as new mathematics. See [PUBLISHABILITY.md](../PUBLISHABILITY.md).

No journal submission, external-author outreach, merge to main, or production release has been performed. This branch and PR expose completed research for review without changing users' installed counting behaviour.
