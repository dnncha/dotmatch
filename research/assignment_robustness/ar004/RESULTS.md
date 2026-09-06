# AR004: shared reads can create false confidence in CRISPR guide evidence

Executed 6 September 2026. Not peer reviewed. Research-only; no production counting changes or journal submission.

## Main finding

Using complete, verified empirical Yusa read-event classes, the **actual unchanged RRA executable distributed with MAGeCK 0.5.9.5** produced at least one reported FDR<=0.05 discovery in **81/100 all-null shared-record experiments**, versus **9/100 independent-event experiments** with identical per-guide count totals and exact uniform marginal guide p-values. The corresponding Brunello comparison was **3/100 in both arms**.

These are controlled component experiments, **not the biological false-discovery rate of the full MAGeCK pipeline, and not a claim that 81% of genes or papers are wrong**. They show an integration failure mode: uniform marginal guide p-values need not give valid gene-level significance when apparently separate guide evidence shares records.

## Actual RRA component: all 400 prespecified runs completed

| Design | Arm | Null trials with any reported FDR<=0.05 discovery | Exact binomial 95% interval | Mean discoveries/trial |
|---|---|---:|---:|---:|
| Yusa plasmid ERR376998 | Shared records | 81/100 | 71.9–88.2% | 2.56 |
| Yusa plasmid ERR376998 | Independent events | 9/100 | 4.2–16.4% | 0.12 |
| Brunello plasmid SRR8297997 | Shared records | 3/100 | 0.6–8.5% | 0.04 |
| Brunello plasmid SRR8297997 | Independent events | 3/100 | 0.6–8.5% | 0.03 |

Every retained guide has the exact randomized binomial null p-value `P(K<k)+U*P(K=k)`. The record-aware arm gives one fair label to every record, shared by all its counting events; the counterfactual independently labels events while retaining the same per-guide totals. Positive-count guide populations are fixed before randomization. The mean fraction of guide p-values<=0.05 was 5.0029% versus 4.9839% in the Yusa arms. The independent arm's 9% result has substantial Monte Carlo uncertainty and does not establish exact 5% component calibration.

The actual command was `RRA -i input.tsv -o output.tsv -p 0.1 --permutation 100`. All 400 full input tables, output tables and logs are preserved in the downloadable research evidence. The protocol amendment was frozen before these results. The Brunello outcome is retained as important negative evidence against a universal-failure claim.

## Exact dependence and the simpler technical null

For gene g, let m_rg count its guide rows that receive record r. Then C=sum(m) is count events, U=sum(1[m>0]) is distinct supporting records, and Q=sum(m^2) is the exact conditional variance of the difference between two fair random record-label groups. A naive event-independent variance uses C instead. The ratio Q/C describes this technical-null model, not biological variance.

Two identical guide count tables, with the same total input-record budget, can have different Q. The tests explicitly demonstrate this. A count table and an overall input-read total therefore cannot reconstruct the missing dependence. This is an elementary counterexample, not claimed as a new theorem.

| Archive | Eligible original gene annotations, C>=100 | Q/C>=1.25 | Q/C>=1.5 | Q/C>=2 |
|---|---:|---:|---:|---:|
| Yusa plasmid | 18,864 | 3,050 | 2,074 | 184 |
| Yusa ESC1 | 18,287 | 2,798 | 1,938 | 167 |
| Brunello plasmid | 19,113 | 18 | 13 | 0 |

In 2,000 label experiments per archive, the fraction with any BH<=0.05 result was:

| Archive | Naive normal test on shared records | Record-covariance normal test, same labels | Independent-event counterfactual |
|---|---:|---:|---:|
| Yusa plasmid | 73.00% | 4.85% | 5.35% |
| Yusa ESC1 | 73.10% | 4.85% | 4.15% |
| Brunello plasmid | 4.45% | 4.05% | 4.05% |

These simple tests are **not MAGeCK or RRA**. The covariance correction is not an implemented RRA repair. Under the artificial all-null experiment, probability of any rejection equals its all-null FDR; it is not the fraction of all genes rejected or a biological FDR estimate.

The structurally selected extreme ZFP59 in Yusa plasmid has 3,376 count events, 1,118 supporting records and Q=11,440 across five populated guide rows. Exact weighted-binomial convolution gives 29.08% rejection for a naive nominal-5% test, versus 5.07% for the covariance-normal approximation and 4.85% for an exact discrete test. This is an error-model example, not a ZFP59 function discovery. All genes and coverage-matched controls remain in the evidence.

## Fresh full ESC2 and actual full MAGeCK

The new ERR377000 archive contains **10,820,594 records**. It was acquired completely from ENA and checked by MD5, SHA-256, bytes and record count. Actual native DotMatch 0.5.0, actual unchanged guide-counter v0.1.3 and an independent matcher reconciled **349,748 per-guide count cells**; 200 seeded records were exhaustively checked against all reference guides, with zero disagreements. New measured reuse: 758,109 multi-guide records, 750,146 within-gene cases, 797,521 extra events.

The full MAGeCK 0.5.9.5 pipeline was run on the original plasmid plus ESC1 and ESC2, with identical median normalization, FDR adjustment and median gene-LFC options. Only ESC2 was freshly counted in this phase; older full plasmid/ESC1 counts were reused after verification. There are two cellular samples and one plasmid reference, not a replicated control arm or a new external study.

| Counting configuration | Negative-selection classifications, reported FDR<=0.05 | Positive-selection classifications, reported FDR<=0.05 |
|---|---:|---:|
| Fixed-window exact | 369 | 2 |
| Fixed-window radius-one | 357 | 2 |
| Fixed-window best-distance | 357 | 2 |
| Native multi-offset guide-counter events | 397 | 11 |

Fixed best-distance versus multi-offset events: 348 common negative classifications, nine fixed-only and 49 event-only, hence **58 changed negative classifications**, plus **nine changed positive classifications**. Negative-rank Spearman correlation remains 0.97269; median absolute gene LFC difference is 0.03329 log2. The full 19,149-annotation population is retained in all paired outputs. No changed classification is declared biologically true or false. Extraction, multiplicity and normalization inputs all differ, so this is not an isolated causal estimate of read reuse.

## Post hoc worked example: apparently independent support survives replication

Two ABCB1B guide sequences differ by 14 positions in the same frame but are shifted overlaps. Their actual multi-offset rows contain:

| Sample | Guide A events | Guide B events | Same records in both |
|---|---:|---:|---:|
| Plasmid | 110 | 111 | 110 |
| ESC1 | 848 | 844 | 837 |
| ESC2 | 140 | 144 | 140 |

The gene's full-MAGeCK median LFC is +0.24354 with fixed-window best counting and +0.23745 with multi-offset events, yet reported positive-selection FDR changes from **0.650558 to 0.006188**. The two guide rows demonstrably reuse almost all supporting records in both cellular samples. This example was selected after inspecting results. It is not a new ABCB1B mechanism, proof the call is false, or a complete causal decomposition of the FDR difference.

## Verification completed

- Ten new focused unit tests passed, including count-table non-identifiability, cross-gene covariance, exact laws, checksum and missing-provenance failures.
- All 400 actual RRA inputs and outputs were independently reread and every reported discovery list, population and summary reconstructed.
- Twelve first/middle/last RRA runs across both designs and arms reproduced full output bytes.
- The complete 6,000 technical-null label experiments reproduced all 19 result files from a clean source copy, using the same verified archived evidence.
- All four actual full MAGeCK analyses reran from a clean copy of the pinned original source/binaries on the same count tables. All eight gene/guide summary files reproduced byte-for-byte.
- Exact randomized marginal-uniformity identities were checked at 42 count/threshold settings.

These replays are computational reproducibility, not additional biological replication. An initial short-timeout technical attempt was interrupted and excluded; only completed manifests are evidence. Original software source was not edited. Build and local execution environments are recorded separately.

## Prior art, claim boundary and provenance

[Barcas](https://doi.org/10.1186/s12859-016-1326-9) already studied Yusa library similarity, shifting and imperfect mapping. [MAGeCK's primary paper](https://doi.org/10.1186/s13059-014-0554-4) describes alpha-RRA and guide-to-gene permutation. [NIH's example](https://hpc.nih.gov/apps/mageck-vispr.html) maps the three original samples; its subsampled files were not used here. Grouped read assignment and elementary covariance are established concepts. No global novelty claim, invalidation of the original biological paper, therapeutic target discovery or general indictment of MAGeCK follows.

The candidate contribution is the empirical chain from source-record reuse to dependence, component calibration and actual replicated-screen sensitivity. Stable counts or effects and replicated guide agreement are insufficient checks when the same underlying evidence has been copied across guide rows. Biological truth, error/overdispersion generalization, context-aware comparators and independent replicated studies remain publication gates.

Primary protocol: `c41be48cbea1a330ed93c3dae79922f5824e5f67`; RRA amendment: `9b090cf74123fa881440157b27a5f3e5248e5edd`.

[Verified ESC2 and pinned-source acquisition/build](https://github.com/dnncha/dotmatch/actions/runs/34035284814). MAGeCK source archive SHA-256 `b06a18036da63959cd7751911a46727aefe2fb1d8dd79d95043c3e3bdaf1d93a`; executed RRA binary SHA-256 `719fa93279388f875202bf1f2811109f72f3b04751c3026a7a36c5baa48d73b8`. New ESC2 SHA-256 `cc3bced7a6cc524a4f81d5bcf337fe0cba6c2e5eb4858d1743636f22d380cf81`. Complete source, outputs, full commands, all selection frequencies and checksum manifests accompany the review packet. No production merge, journal submission or external author outreach was performed.
