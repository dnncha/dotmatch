# Counting reads or counting matches?

## AR-001 exploratory findings — 6 September 2026

**Status: completed counting-methods investigation; not a completed gene-discovery study or a submission-ready manuscript.** The protocols were committed before their respective new experiments. The initial historical summaries were already known; the exact-overlap descriptor was developed after observing the prefix results. This is not an externally preregistered study. Production DotMatch code was not modified.

## Main finding

Across two complete Yusa sequencing archives containing **20,394,663 read records**, guide-counter 0.1.3 produced **20,270,172 matching events from 18,770,788 distinct contributing read records** under the explicitly tested multi-offset configuration. Independent reconstruction reproduced every comparator guide count, and a separate independent decoder reproduced all three fixed-window DotMatch policy count matrices.

The difference is **1,499,384 additional matching events above one event per contributing read**, or 7.99% of the contributing-read denominator. This establishes a difference in counting units, not a biological error rate. Importantly, 20,270,172 is below the total number of input reads: a simple check that total counts do not exceed input records would not detect the multiplicity in the combined full-archive result.

Most affected Yusa read records link different guide IDs with the same gene annotation. In the selected Brunello prefixes, far fewer records have multiple matching offsets, but most of those records link different gene annotations. Thus the magnitude of a guide-count discrepancy alone is not a measure of its gene-level consequence. Whether any of these differences changes a supported biological conclusion remains untested.

## Definitions and scope

A **read record** is one complete record in the input FASTQ; it is not necessarily a distinct original molecule. A **matching event** is a unique target match at one selected offset within that read. One read can contribute multiple events. A **multi-offset record** contributes at least two such events. **Cross-gene** means its matched guide IDs have more than one gene annotation in the supplied reference; annotation is not proof of molecular origin.

There are separate, overlapping analysis lanes. The six-sample fixed-window pilot evaluates records 2,001–102,000 after using records 1–2,000 for extraction diagnostics. Historical reconstructions use records 1–100,000, including calibration. The full Yusa lane uses complete archives, including the prefix records. Do not add these analysis denominators together as independent reads. No complete Brunello archive was processed in this investigation.

## 1. The full Yusa result survives complete-archive verification

| Full archive / sample | Input read records | Records contributing an event | Matching events | Additional events | Multi-offset records | Cross-gene multi-offset records |
|---|---:|---:|---:|---:|---:|---:|
| ERR376998 / plasmid | 10,093,905 | 9,404,898 | 10,161,871 | 756,973 | 720,995 | 7,445 |
| ERR376999 / ESC1 | 10,300,758 | 9,365,890 | 10,108,301 | 742,411 | 707,092 | 8,040 |
| Combined | 20,394,663 | 18,770,788 | 20,270,172 | 1,499,384 | 1,428,087 | 15,485 |

All 1,428,087 multi-offset records link different guide IDs; 1,412,602 (98.92%) link guides sharing one gene annotation, and 15,485 (1.08%) span annotations. These are complete-archive descriptive counts, not estimates from the bounded diagnostic traces. The saved traces contain only the first at most 1,000 multi-offset records per full archive and are explicitly not random samples.

Both compressed archives matched ENA's byte count and MD5, and the number of parsed FASTQ records matched ENA metadata. Their additional SHA-256 hashes are:

- ERR376998: `7f79b76cec12b70319744417282f963c00818a5f0ae61497bd7b64790ac55f2f`.
- ERR376999: `cf2bc10938e178d16dfb81ca2f9fda805cae892290fac3cfc243bb637f8cda17`.

The comparator selected zero-based offsets 21, 22, 23 and 24 from the first 100,000 records of each file, with `--offset-min-fraction 0.0025` explicitly supplied. This is the tested setting, not a claim about every possible configuration or another version. Calibration records remain in the full-file event counts to reproduce the comparator's algorithm.

The aggregate identity is `20,394,663 input records - 1,623,875 records with no selected-offset event + 1,499,384 additional events = 20,270,172 events`. Unmatched records can therefore mask repeated contributions when only grand totals are examined.

## 2. Fixed-window policy differences are independently reproducible

At fixed zero-based start 23 and length 19, all policies use the same window and reference. Exact requires a unique exact match; radius-one requires exactly one reference target within one substitution; best-distance requires a unique nearest target within one substitution. An exact hit with a distance-one alternative is unique under best-distance but ambiguous under radius-one.

| Full Yusa sample | Exact unique reads | Radius-one unique reads | Best-distance unique reads | Records whose status or unique target changes across policies |
|---|---:|---:|---:|---:|
| Plasmid | 8,615,587 | 8,910,700 | 8,980,661 | 435,399 |
| ESC1 | 8,475,790 | 8,874,959 | 8,948,925 | 547,420 |
| Combined | 17,091,377 | 17,785,659 | 17,929,586 | 982,819 |

The changed-record denominator is all 20,394,663 input records. Best-distance adds 838,209 unique counts relative to exact; radius-one withholds 143,927 exact-hit records that have a distance-one alternative. Another 683 nonexact records become ambiguous rather than unique. These mechanisms reconcile to the 982,819 changed records. None of these quantities identifies how many biological assignments are correct.

The comparator's multi-offset event matrix and the fixed-window count matrices intentionally answer different questions. Their entire difference must not be called double counting: it also includes extraction coverage and assignment semantics. Full-run membership tables explicitly retain records counted by both approaches, only by the multi-offset event approach, only by fixed-window best-distance, or by neither.

## 3. Both historical discrepancies reproduced exactly

| Historical lane | Prefix records | DotMatch counts | guide-counter events | Guides with differing totals after summing samples |
|---|---:|---:|---:|---:|
| Yusa, two prefixes | 200,000 | 184,167 | 208,700 | 13,537 / 87,437 |
| Brunello, four prefixes | 400,000 | 349,184 | 350,374 | 255 / 77,441 |

Every guide-counter count in every named sample equals an independent reconstruction; the historical aggregated numbers are not merely approximately similar. Sample-specific matrices are retained because aggregation can conceal opposing changes.

The event/read reconciliation is:

| Historical lane | Events | Distinct contributing records | Additional events | Multi-offset records | Same-gene multi-offset records | Cross-gene multi-offset records |
|---|---:|---:|---:|---:|---:|---:|
| Yusa | 208,700 | 193,205 | 15,495 | 14,790 | 14,645 | 145 |
| Brunello | 350,374 | 349,415 | 959 | 881 | 109 | 772 |

Only 0.98% of Yusa multi-offset records span gene annotations, versus 87.63% of Brunello multi-offset records. These proportions concern selected archival prefixes and a conditional multi-offset denominator; they are not whole-Brunello rates, error rates, or statistical estimates across independent studies.

Algebraically, the historical Yusa count difference of 24,533 comprises 15,495 additional events plus a 9,038 difference between distinct contributing-record counts. The Brunello difference of 1,190 similarly comprises 959 additional events plus 231. This is a count reconciliation, not complete per-read causal attribution of the remaining assignment differences.

## 4. Positional overlap helps explain the discrepancy

A conventional same-offset Hamming-one audit found 265 unordered guide pairs in Yusa (264 cross-gene) and 459 in Brunello (426 cross-gene). That does not capture matches created by extracting neighbouring windows at different read positions.

We therefore implemented a post-hoc exact-overlap descriptor: for shifts `s` in 1, 2 or 3, retain the directed guide pair A, B when `A[s:] == B[:-s]`. This expresses exact sequence compatibility between an earlier and a later window. It is not a complete one-mismatch ambiguity graph and does not model vector flanks or establish true origin.

| Full reference library | Directed exact-overlap edges | Guide IDs participating | Cross-gene edges |
|---|---:|---:|---:|
| Yusa, 87,437 unique 19-base targets | 5,282 | 9,522 | 19 |
| Brunello, 77,441 unique 20-base targets | 225 | 399 | 197 |

In the historical prefixes, at least one exact short-shift overlap pair is present in **14,678 of 14,790 Yusa multi-offset records (99.24%)** and **749 of 881 Brunello multi-offset records (85.02%)**. This is structural support for an explanation of the observed matching pattern, not held-out predictive validation. The conventional Hamming graph is unordered, whereas the overlap graph is directed; their edge totals should not be interpreted as interchangeable density measures.

The overlap descriptor was chosen after examining the prefix findings. We have not evaluated a prospectively locked predictor on held-out libraries, and have not computed full-archive overlap-support prevalence from the bounded full-archive traces.

## 5. The six-sample pilot is an extraction control, not a hit-calling dataset

| Evaluation prefix | Selected zero-based start | Read records | Exact unique | Radius-one unique | Best-distance unique | Changed records |
|---|---:|---:|---:|---:|---:|---:|
| Yusa plasmid | 23 | 100,000 | 89,834 | 91,587 | 92,353 | 3,289 |
| Yusa ESC1 | 23 | 100,000 | 88,878 | 91,056 | 91,816 | 3,698 |
| Brunello plasmid | 28 | 100,000 | 11,965 | 13,068 | 13,088 | 1,143 |
| Brunello RepA | 26 | 100,000 | 10,360 | 11,559 | 11,578 | 1,237 |
| Brunello RepB | 25 | 100,000 | 10,242 | 11,335 | 11,350 | 1,123 |
| Brunello RepC | 23 | 100,000 | 11,602 | 12,857 | 12,895 | 1,331 |

Brunello guide positions are mixed. The low assigned fraction at one selected window is therefore not a complete-assay performance result. All three policy matrices, all four state totals, all 48 transition rows, and every changed-record call/ordinal matched the independent oracle in each of the six prefixes. Discovery records were excluded from evaluation. These shallow, nonrandom prefixes were not used for gene-level p-values or hit calling.

## Methods, validation and evidence

The production baseline is DotMatch commit `11d159fa1648365f2a4e96917b483c33aa5d9fe7`, version 0.5.0. The external comparator is guide-counter 0.1.3, installed with Cargo's locked dependency configuration. Its downloaded crate records source commit `5e4973017f8a022226c7b003d082c8c456c91ef3`. The comparator's CLI lacks a version flag; installed-package provenance, binary hashes and source hashes are retained instead.

The pilot oracle enumerates query substitutions against a sequence-to-all-reference-IDs dictionary without importing DotMatch assignment code. It is checked against exhaustive Hamming distance on adversarial fixtures and all 625 ACGTN four-base queries across five small libraries. Selected real windows are additionally compared with every reference target. The full-run decoder uses a separate enumerated reference codeword table, checked against the query oracle and brute-force distance. Correct handling of exact-plus-neighbour cases, ties, reference duplicates, literal N, invalid windows, named sample axes and corrupted evidence is covered by regression tests.

All generated matrices preserve zero-count guides and original gene annotations. Read and event conservation, named sample axes, archive metadata, file hashes and equality of every guide count are required before a run is marked complete. The offline summary verified **190 manifest-listed files** across five successful evidence bundles, as well as reference consistency between analysis lanes. This is independent implementation validation within the project, not review by an independent research group.

Executed evidence runs:

- Six-sample controlled pilot: [34026563667](https://github.com/dnncha/dotmatch/actions/runs/34026563667), commit `923a3290f1fe733983aff8b1ba3fc8978184e851`.
- Historical Yusa reconstruction: [34026723503](https://github.com/dnncha/dotmatch/actions/runs/34026723503), commit `cc00e5b0c18202ab6e19ddbe790cb55d8e31c87f`.
- Full Yusa archives: [34027101859](https://github.com/dnncha/dotmatch/actions/runs/34027101859), commit `2780baab967f4206d196c9cb3a7ce30588966afb`.
- Historical Brunello reconstruction: [34027183669](https://github.com/dnncha/dotmatch/actions/runs/34027183669), commit `25df10480a220431b7faff5339c039f07d77b2ed`.
- Complete-reference overlap audit: [34027309060](https://github.com/dnncha/dotmatch/actions/runs/34027309060), commit `8aaeca7d212b8eb0e78eede6b93e7860185b7d8a`.

Two initial research-harness failures remain disclosed in [AMENDMENTS.md](AMENDMENTS.md): the Brunello `Seq` header was initially misread as data, and a comparator version probe used an unsupported flag. Correcting them did not alter reference sequences or counting rules. No source rows were dropped to obtain agreement.

## What remains before biological or publication claims

The full Yusa inputs provide one plasmid and one cellular library, not replicated treatment groups. The Brunello inputs here are prefixes, not complete biological samples. No gene-level significance, changed hit list, new pathway, improved biological assignment accuracy, or validated uncertainty bound has been demonstrated.

The next study stage must lock a complete assay-aware extraction method, full-run cohort manifest and downstream analysis plan before examining policy-dependent gene results. Biological replicates must remain distinct; technical lanes must not become extra biological replicates. Matched-depth analyses, stable downstream settings, held-out studies and independent biological evidence are needed to judge consequences rather than only counting differences.

The scientific contribution to pursue is an audit of **joint offset/guide ambiguity and counting-unit conservation**, with a demonstrated connection to downstream robustness. It is not the generic claim that assignment decisions matter. Related prior work includes [crispat, Bioinformatics 2024, btae535](https://doi.org/10.1093/bioinformatics/btae535), which compares guide-to-cell strategies in single-cell screens, and [ReCo, Bioinformatics 2023, btad448](https://doi.org/10.1093/bioinformatics/btad448), which addresses read counting from staggered CRISPR sequencing. The Brunello source study is [Sanson et al., Nature Communications 2018, 5416](https://doi.org/10.1038/s41467-018-07901-8). These references are a starting novelty check, not a completed systematic review.

Code, protocols, aggregates and checksums should be durably deposited before submission. Current Actions artifacts expire after 30 days; the downloadable evidence bundle preserves copies but is not itself a DOI deposit. Raw sequencing reads are not redistributed. Final scientific review and authorship approval remain outstanding; AI-assisted implementation and analysis should be disclosed under the eventual venue's policy. No journal submission, upstream accusation, production merge or release was made.
