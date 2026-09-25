# When a read becomes several guide counts

Completed computational measurement study, 25 September 2026. This is a research result, not a new release or a biological-accuracy claim. The complete manuscript, figures, preserved outputs and all postprocessing source were delivered as a reproduction capsule; the core results can also be recomputed with `verify_core_results.py` here.

## Data and design

All **62,831,612** records in six complete public FASTQs were processed, using **87,437** unique 19nt targets from the unmodified MAGeCK-distributed Yusa library. The original experiment is Koike-Yusa et al., Nature Biotechnology (2014), DOI [10.1038/nbt.2800](https://doi.org/10.1038/nbt.2800), ENA PRJEB4038 / ERP003292. No toxin-selected or drug-resistance samples were analyzed.

| Material | Technical set A | Technical set B | Shared biosample |
|---|---|---|---|
| Plasmid | ERR376998: 10,093,905 reads | ERR377001: 10,246,121 reads | SAMEA2188756 |
| ESC1 | ERR376999: 10,300,758 reads | ERR377002: 10,393,396 reads | SAMEA2188757 |
| ESC2 | ERR377000: 10,820,594 reads | ERR377003: 10,976,838 reads | SAMEA2188758 |

A and B resequence the same three biosample identities. They are **not additional biological replicates**. The distributed reference list is not identical in size to the nominal 87,897-guide library described in the source paper. The analyzed CSV SHA-256 is `d41a4122a46fdedb47deca61081dde9037d461a1fce1e6f9744522e21005ceba`.

The protocol and core code were frozen before the new full-run outcomes, at [5872c2ba1b6b](https://github.com/dnncha/dotmatch/tree/5872c2ba1b6b92a5525d7d72373f5a8f1724c859/research/assignment_study). This is an exploratory study, not a prospective preregistration: earlier small-prefix count discrepancies were already known. Fixed window: zero-based offset23, length19. Multi-offset calibration: first100,000 records, retain an offset with >=0.0025 of matched windows. Exact/radius-one/nearest-one and per-offset/distinct-guide/joint policies were kept separate.

## The primary hypothesis was falsified

There were **zero repeated contributions to the same guide from a single read**, in either exact or one-mismatch multi-offset counting, across all six files. Counting each specific guide only once per read therefore changed no count and no gene estimate. This null result is retained, not re-labelled as a positive discovery.

The multiplicity instead involved different guides:

| One-mismatch multi-offset accounting | Observed total |
|---|---:|
| Input records | 62,831,612 |
| Records with at least one guide match | 57,278,925 |
| Guide-count contributions | 61,871,230 |
| Records supporting multiple distinct guides | 4,371,978 |
| Multi-guide records involving one gene label | 4,324,884 |
| Multi-guide records involving more than one gene label | 47,094 |
| Repeated same-guide contributions | 0 |

Thus 6.958% of all records supported multiple guide identities, and 98.923% of these multi-hit records involved guides annotated to the same gene. These are **assignment relationships, not error rates**. Contributions can exceed matched records while remaining below total input records. In ERR376998 alone, 10,161,871 contributions exceeded its 10,093,905 input records.

An exploratory library audit found 3,463 directed pairs with exact 18-base overlap after a one-base shift; 3,461 shared a gene label. In pooled one-mismatch ESC1 counts, 3,603 same-gene guide pairs shared at least half of the smaller read set, requiring >=60 reads for both guides; they involved 3,017 genes. ESC2 had 3,622 such pairs involving 3,033 genes. Read-set intersections were reconstructed from every multi-hit ledger row, not inferred from count correlations. This does not establish interchangeability of biological perturbations.

## Separate fixed-window mismatch sensitivity

Fixed-exact assignment retained 52,069,456 reads; unique nearest within one mismatch retained 54,711,633, adding 2,642,177. Radius-one unique assignment retained 54,265,906. Radius-one ambiguity was 448,044 records; nearest-one ambiguity was 2,317. Neither policy was assumed to identify biological origin.

The primary downstream endpoint was a descriptive median-guide log2 abundance ratio, ESC/plasmid. Common median-ratio size factors were estimated from fixed-exact guides positive in all three pooled samples. Pseudocount1 in normalized count space; eligible guides had >=60 fixed-exact plasmid counts; eligible genes required >=3 guides. There were 73,567 eligible guides and 16,654 eligible genes. A flag means |change in median-guide log2 ratio|>=0.5, **not statistical significance or a newly discovered dependency**.

| Comparison | Status | ESC1 flags | ESC2 flags |
|---|---|---:|---:|
| Same-guide deduplication, one mismatch | Primary | 0 | 0 |
| Same-guide deduplication, exact | Primary | 0 | 0 |
| Fixed nearest-one minus exact | Secondary | 178 | 199 |
| Fixed radius-one minus nearest-one | Secondary | 43 | 41 |
| Multi-offset minus fixed nearest-one | Exploratory | 324 | 379 |
| Joint nearest-one minus per-offset nearest-one | Exploratory | 489 | 547 |

Most fixed-window mismatch effects were small: median absolute gene changes were approximately 0.0238 and 0.0237 log2 units. Policy-specific normalization gave 163/191 mismatch flags; excluding calibration records gave 174/197 with common normalization and 163/192 with policy-specific normalization. Held-out analyses had 16,622 eligible genes. The primary null remained exact in all conditions.

RPL4 and EMG1 met the stronger descriptive secondary criterion in both normalization modes: same-sign gene changes >=0.5 in both technical sets and both ESC libraries, and >=2 concordantly changed guides per ESC contrast in pooled data. These are not novel biological findings.

| Gene | ESC1 exact | ESC1 nearest-one | ESC2 exact | ESC2 nearest-one |
|---|---:|---:|---:|---:|
| RPL4 | -6.052494 | -4.981827 | -6.581786 | -6.018328 |
| EMG1 | -4.764940 | -1.829096 | -5.932256 | -1.880356 |

The values are median-guide log2 abundance ratios, not fitness estimates. Pseudocount sensitivity is material near zero; the capsule preserves values for 0.5,1,2,5 without selecting the most dramatic result.

## Post-result sequence-quality inspection

An explicitly [post-result follow-up](https://github.com/dnncha/dotmatch/tree/4f00844296b988a1f82599f04d8260862475b80e/research/assignment_followup) re-read the SAME six files, not a new cohort. It examined all ten RPL4/EMG1 guides and thirty controls selected by plasmid exact abundance, GC content and a fixed random seed, without selecting on control ESC results.

For RPL4 guide `chr9:64176974-64176997`, ESC1 contained **zero exact reads and 76 reads of `CAAAATGAGAAACCGAGGG`**, one mismatch at guide position17. Technical A/B counts were 40/36. All76 mismatch bases had reported Phred>=30; mean35.895. Plasmid had4 copies of that alternate sequence and323 exact-reference reads. These are amplified-read observations, not independent molecules.

For EMG1 guide `chr6:124705610-124705633`, ESC2 had zero exact reads and26 one-mismatch reads;25 shared `CAAGCCAGCTCAGTTCAGG`, differing at position12, and24 of those25 had reported quality>=Q30.

This pattern was not unique to the selected genes. Four of the thirty matched controls also had a dominant alternate sequence observed >=20 times with >=90% of its mismatch bases at Q30 or above in an ESC library. That descriptive rule was applied after inspection and is not a formal case-control test.

**High reported quality does not prove biological guide mutation, loss of guide activity, selection for escape, or sequencing-error correction.** Library/synthesis variants, PCR errors and amplification, contamination, reference limitations and systematic base-calling errors remain alternatives. The scientific output is an inspectable sequence-level distinction, not a claim that the more permissive policy is biologically better.

## Checks actually completed

- Full compressed-file byte counts, ENA MD5s, SHA-256s and all FASTQ record counts verified.
- Actual native DotMatch and actual pinned guide-counter executed; 30 full count-vector comparisons with an independent Hamming reference passed: 2,623,110 entries. This is not exhaustive per-read biological validation.
- 6,000 real-read full-library brute-force checks passed, systematically/stratified sampled after calibration. A separate synthetic fixture checked eight policies on194 records.
- All six full multi-hit ledgers independently re-read and primary intervention count changes reconstructed.
- An independent standard-library statistical implementation checked 482,966 numeric outputs; maximum export-rounding discrepancy <5e-10.
- Six statistical helper tests passed. The public core checker here separately recomputed all principal pooled results and passed.
- Sequence follow-up: all480 selected-guide exact/corrected comparisons with the first study passed. An independent synthetic quality fixture checked130 output cells on94 records, including positions, sums, thresholds and literal-N handling.

## Reproduce and inspect

Successful primary run: [36191069868](https://github.com/dnncha/dotmatch/actions/runs/36191069868), source `a5bde2d171f8e71d89851a5255e130d0ee345b48`. Primary artifact SHA-256: `295a66d127f1fb74315ea40990bec9a106c05b1f08194cb4eb6a8f2036d3dff1`.

Successful sequence follow-up: [36192950910](https://github.com/dnncha/dotmatch/actions/runs/36192950910), source `4f00844296b988a1f82599f04d8260862475b80e`. Follow-up artifact SHA-256: `9d938328367242199aaf28c696daf3008fc25249c294a9810fcb2c5bc828546f`.

Native DotMatch base: `5c3934d648ff55928db5a9219d7703845941d28e`. guide-counter: `f24de175282064773328e35baf03373e7e3b895d`. Full original outputs and manifests are in the delivered capsule; CI artifacts expire and are not a permanent scientific archive.

```bash
python3 verify_core_results.py --primary /path/to/unmodified-primary-output
```

The checker uses Python's standard library, verifies the input output-manifest and recalculates the principal outcomes. To regenerate raw results, execute `research/assignment_study/execute.py` at the pinned primary source; then the pinned follow-up driver. The first study attempt failed before analysis with HTTP403; checksum-verified mirrors fixed retrieval without changing the scientific rules.

## Claim boundaries and disclosure

This completed computational study is one dataset and three source biosamples with technical re-sequencing. It did not run a replicated, validated gene-hit-calling experiment, establish FDR, validate molecular origin, or correct the source paper's conclusions. Matching and shifted extraction are longstanding methods (e.g. [Dai et al. 2014](https://doi.org/10.12688/f1000research.3928.2)); novelty priority for the general concept is not claimed. The owner authored DotMatch and has an interest in adoption. AI assistance was used for code, analysis and drafting. No external peer review, wet-laboratory validation, journal/preprint submission or publication acceptance is claimed. Main and production code were not modified by this research branch.
