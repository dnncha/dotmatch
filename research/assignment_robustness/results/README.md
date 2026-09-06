# One read, several guides

## Technical summary

**Executed evidence, 6 September 2026:** independent full-read audits of 30,215,791 reads across three public archives and two library/construct designs. Every fixed-window guide count reconciled with pinned DotMatch 0.5.0; every multi-offset count reconciled with the unchanged guide-counter v0.1.3 executable. No assignment disagreements were found.

In the two Yusa archives, 1,428,087 reads entered more than one guide row. Of these, 1,412,602 (98.92%) involved guides labelled with the same gene. There were 1,499,384 extra count events. These are computational read-reuse measurements, not estimates of incorrect biological origins.

**Interpretation:** apparent support from several guide rows can reuse the same reads. Fixed-position Hamming separation alone does not diagnose aliases created by shifted extraction windows. This is an assay-dependent effect, not evidence that every screen or every guide-counting program is affected.

## Scope and definitions

ERR376998 and ERR376999 are the Yusa plasmid/ESC1 tutorial pair, not independent replicated treatment arms. SRR8297997 is a Brunello plasmid library run, not a treatment contrast. The analysis counts sequencing records, not deduplicated original molecules; no UMI deduplication or gene-level hypothesis testing is performed.

Let N be input reads, M reads with at least one accepted window, C retained guide-window count events, U=N-M unmatched reads, and E=C-M extra events. Then C-N=E-U. A count-table total below the input-read total does not rule out repeated counting. A multiple-guide read is a record with at least two distinct accepted target IDs; it is not necessarily a cross-gene ambiguity.

## Complete-archive results

| Archive | Reads N | Matched reads M | Count events C | Extra events E | Multiple-guide reads | Same-gene cases | Cross-gene cases |
|---|---:|---:|---:|---:|---:|---:|---:|
| ERR376998 | 10,093,905 | 9,404,898 | 10,161,871 | 756,973 | 720,995 | 713,550 | 7,445 |
| ERR376999 | 10,300,758 | 9,365,890 | 10,108,301 | 742,411 | 707,092 | 699,052 | 8,040 |
| SRR8297997 | 9,821,128 | 9,289,214 | 9,318,181 | 28,967 | 26,561 | 3,436 | 23,125 |

For Yusa, the pooled count-event total is 20,270,172, below 20,394,663 input reads, despite 1,499,384 excess events. Unmatched reads conceal the duplication in this aggregate check.

The table is a census of the retrieved archives under the specified rules. The sample-to-sample comparison is descriptive; treating millions of reads as millions of biological replicates would be invalid.

## Discovery-only construct model and held-out evaluation

Each sample uses its first 100,000 records solely for offset discovery and modal-flank estimation. Tied flank bases become N. Each reference guide is inserted into that inferred context, and its accepted guide-window events are recorded. The entire model is written before later records are evaluated. The prediction target is whether a read produces multiple count events—not its true biological source.

| Archive | Templates predicting multiple events | Evaluated exact-window reads | Excluded evaluation reads | TP | FP | FN | TN |
|---|---:|---:|---:|---:|---:|---:|---:|
| ERR376998 | 7,057 | 8,525,759 | 1,468,146 | 679,336 | 164 | 66 | 7,846,193 |
| ERR376999 | 7,057 | 8,386,903 | 1,813,855 | 663,212 | 164 | 63 | 7,723,464 |
| SRR8297997 | 131 | 56,204 | 9,664,924 | 476 | 0 | 8 | 55,720 |

The model is evaluated only where the chosen fixed window has a unique exact target. Excluded records remain in the full-read audit. This restriction and the near-constant construct sequence make the mechanistic prediction task much easier than biological-origin inference; no general assignment-accuracy claim follows.

## Error-free interventions isolate the mechanism

The following controls contain constructed reads with known source guide labels. They are not newly measured biological samples. Balanced controls contain exactly one error-free template per reference guide. The Yusa PHF23 source-only witness was selected after inspecting the real-data findings and is explicitly a post-hoc mechanistic test.

| Reference/sample | Synthetic experiment | Upstream rule | Known source reads | Count events | Zero-origin guides receiving counts |
|---|---|---|---:|---:|---:|
| ERR376998 | balanced | one_mismatch | 87,437 | 94,865 | 0 |
| ERR376998 | balanced | exact | 87,437 | 89,074 | 0 |
| ERR376998 | witness_control | one_mismatch | 128 | 256 | 1 |
| ERR376998 | witness_control | exact | 128 | 256 | 1 |
| ERR376998 | witness_depleted | one_mismatch | 16 | 32 | 1 |
| ERR376998 | witness_depleted | exact | 16 | 32 | 1 |
| ERR376998 | witness_flank_changed | one_mismatch | 128 | 256 | 1 |
| ERR376998 | witness_flank_changed | exact | 128 | 128 | 0 |
| ERR376999 | balanced | one_mismatch | 87,437 | 94,865 | 0 |
| ERR376999 | balanced | exact | 87,437 | 89,074 | 0 |
| SRR8297997 | balanced | one_mismatch | 77,441 | 77,441 | 0 |
| SRR8297997 | balanced | exact | 77,441 | 77,441 | 0 |

In the selected witness, the two 19-base targets differ at 17 positions when compared in the same frame. Nevertheless the construct flank plus an 18-base target overlap yields an exact match to the second guide at a shifted offset. Changing the last left-flank base removes that exact alias while a one-mismatch match can remain. Fixed-window exact DotMatch counts agree with the known source counts in every synthetic control.

An eightfold change in source-only synthetic reads produces the same apparent change in the alias guide although that guide contributes no source reads. This is a demonstration of duplicated guide-level evidence, not a new gene-function finding or a calibrated false-discovery-rate result.

## Methods and independent checks

The full input archives were verified against ENA byte counts and MD5s, then locked with SHA-256. The corrected references retain target and gene identities. Fixed-window policy comparisons use zero-based start/length 23/19 for Yusa and 21/20 for Brunello. These are controlled policy-comparison windows, not claims of optimal extraction for staggered reads.

The independent index uses two disjoint exact seeds: any observation within one Hamming substitution must share at least one complete seed. Every candidate is checked by full Hamming distance, with duplicate target identities retained. It does not use DotMatch candidate lists. A 12,500-window constructed grid checks all candidates and all three policies against exhaustive all-target enumeration, including N and duplicate cases. The real-read audit also checks 200 prespecified pseudorandom ordinals per archive against every reference target.

The replay reconciles 756,945 native policy/guide count values and 252,315 upstream guide count values, with 600 exhaustive real-window checks. Each value is derived from the complete archive, not just the oracle sample. Count-table identity, gene labels, all four read states, changed-read totals, and event-conservation identities are checked.

Upstream reproduction intentionally includes its exact-match precedence, last-inserted duplicate-sequence behaviour, ACGT-only mismatch lookup, event-based offset threshold and increment-at-every-accepted-offset loop. These semantics are not silently equated with fixed-window Hamming matching or with true molecular origin.

## Prior art and what could be publishable

Mapping uncertainty is not a new idea: [bcSeq](https://doi.org/10.1093/bioinformatics/bty402) models sequencing-error-aware barcode assignment, and [crispat](https://doi.org/10.1093/bioinformatics/btae535) studies downstream sensitivity to guide-to-cell assignment. [Buschmann and Bystrykh](https://doi.org/10.1186/1471-2105-14-272) already describe why DNA sequence context matters for barcode decoding. The general boundary/uncertainty concepts must not be claimed as novel.

**The strongest candidate contribution is a bulk-screen, read-to-count evidence audit showing when apparent multi-guide support reuses the same records, an assay-context diagnostic, and independent validation that separates this effect from sequencing errors and mismatch rescue.** The current study supports a methods/technical-results draft. A stronger downstream-statistics paper still requires replicated biological contrasts, a locked normalization/filtering plan, held-out studies, and direct tests of hit/rank/effect changes.

## Limitations and unresolved questions

Two library/construct designs are not an atlas of all CRISPR screens. The Yusa results motivated part of the subsequent mechanism work; that adaptation is documented rather than presented as preregistered discovery. A single-read witness cannot establish empirical prevalence, and algorithm agreement cannot certify biological truth. Guide-window aliases can be same-gene without being harmless to guide-level support, but gene-level statistical distortion has not been established here.

The independently generated original pilot bounds remain conditional on candidate-set completeness. They are not confidence intervals and are not used to assert calibrated gene-level coverage in this replay. Likewise, no confidence intervals from millions of read records are substituted for biological replication.

The next decisive question is whether the documented reuse changes replicated guide-support and gene-level conclusions under otherwise identical analysis. Full differential-screen curation, outcome-blind holdouts and appropriate comparator coverage remain publication gates, not completed claims.

## Reproducibility and evidence access

[Completed workflow and complete per-sample artifacts](https://github.com/dnncha/dotmatch/actions/runs/34030217143). Source and protocols live beside this report. Job artifacts contain complete count tables, template models, read-class tables, synthetic controls, commands, source and binary hashes, and failure-visible logs. Public raw FASTQs are not recommitted. Archive transport artifacts expire; original ENA accessions and acquisition scripts remain the source of input data.

The formal portable HTML renderer was unavailable after a local execution-service failure. This Markdown report is the explicitly labelled alternative, not a claim of completed HTML rendering. No production release, manuscript submission, accepted paper or new biological discovery is claimed.
