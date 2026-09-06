# Joint guide/position resolution in pooled CRISPR counting

**AR003 technical results; 6 September 2026. Not peer reviewed.**

The new joint decoder was executed on all **30,215,791 original sequencing records** across ERR376998, ERR376999 and SRR8297997. It considered every permitted position jointly, rather than adding a count for each accepted window. Complete fixed-window baseline counts were independently reconciled with the prior pinned DotMatch 0.5.0 results.

## What is identifiable at each resolution?

| Archive | Policy | Guide-unique reads | Gene-unique reads | Additional gene-identifiable reads |
|---|---|---|---|---|
| ERR376998 | exact | 8,870,635 | 9,031,919 | 161,284 |
| ERR376998 | radius_k1 | 8,623,040 | 9,327,033 | 703,993 |
| ERR376998 | best_k1 | 9,238,068 | 9,408,051 | 169,983 |
| ERR376999 | exact | 8,727,284 | 8,879,967 | 152,683 |
| ERR376999 | radius_k1 | 8,594,091 | 9,283,086 | 688,995 |
| ERR376999 | best_k1 | 9,204,749 | 9,368,855 | 164,106 |
| SRR8297997 | exact | 8,745,892 | 8,746,140 | 248 |
| SRR8297997 | radius_k1 | 9,146,555 | 9,151,127 | 4,572 |
| SRR8297997 | best_k1 | 9,282,354 | 9,282,675 | 321 |

The last column is gene-unique minus guide-unique within the same policy. These are overlapping views of the same records, not totals to add together. A gene-identifiable read may remain ambiguous between guides of the same gene. More identifiable reads do not, by themselves, establish more accurate assignments. Complete ambiguous, unmatched, invalid and position-level budgets are in `read-resolution.tsv`.

## Yusa representation sensitivity

| Joint policy versus fixed best | Eligible annotation groups | Absolute change >=0.5 log2 | Share | Strong sign reversals |
|---|---|---|---|---|
| exact | 18,721 | 70 | 0.374% | 0 |
| radius_k1 | 18,721 | 48 | 0.256% | 0 |
| best_k1 | 18,721 | 12 | 0.064% | 0 |

Eligibility is fixed-window exact plasmid sum >=50 with at least two guides each having >=10 baseline counts. The same population is used for every policy. Ratios use original input-read exposure and pseudocount 0.5. Strong sign reversals are a descriptive secondary check requiring opposite signs with both absolute ratios >=0.5. Thresholds 20/100 and pseudocounts 0.1/1 are reported separately. This changes extraction and resolution as well as assignment; it is not a pure mismatch-policy comparison.

All 19,149 original Yusa annotation labels are retained in the complete tables, including zero counts and ineligible groups. They are original annotation groups, not automatically distinct validated genes. `all-primary-outliers.tsv` contains every qualifying outlier, not only illustrative examples. No gene-level p-values or phenotype discoveries are inferred from this unreplicated pair.

## Known-origin controls

| Archive | Policy | Constructed records | Guide correct | Guide incorrect | Gene correct | Gene incorrect |
|---|---|---|---|---|---|---|
| ERR376998 | exact | 87,437 | 85,824 | 0 | 87,437 | 0 |
| ERR376998 | radius_k1 | 87,437 | 79,949 | 0 | 86,861 | 0 |
| ERR376998 | best_k1 | 87,437 | 85,824 | 0 | 87,437 | 0 |
| ERR376999 | exact | 87,437 | 85,824 | 0 | 87,437 | 0 |
| ERR376999 | radius_k1 | 87,437 | 79,949 | 0 | 86,861 | 0 |
| ERR376999 | best_k1 | 87,437 | 85,824 | 0 | 87,437 | 0 |
| SRR8297997 | exact | 77,441 | 77,410 | 0 | 77,412 | 0 |
| SRR8297997 | radius_k1 | 77,441 | 76,361 | 0 | 76,390 | 0 |
| SRR8297997 | best_k1 | 77,441 | 77,410 | 0 | 77,412 | 0 |

These are the previously archived balanced, error-free constructs, one per reference guide. They test implementation and information loss under known contexts; they do not estimate real sequencing accuracy or generalize to unknown assays. Ambiguous and unassigned results are fully retained in `known-origin-controls.tsv`.

## Validation and evidence boundaries

- **756,945** complete fixed-window guide-count cells reconciled with pinned native evidence, plus all fixed-window read states.
- **600** fixed-seed selected public records checked against every reference target at every permitted position; these are sampled checks, not an exhaustive all-target check on every archive record.
- **758,133** candidate classes independently reaggregated into every guide count, gene lower/upper bound and matched-state total.
- **2,068,092** representation calculations checked through ratio and log-difference arithmetic; old annotation sums independently reconciled through SQL and Python.

Full input SHA-256, MD5, byte count and record count were verified before execution and input SHA-256 checked again afterward. The report verifies each completed artifact checksum before analysis. Hash agreement proves byte identity, not biological accuracy. No unfinished local execution was used as evidence.

Conditional gene-count and effect ranges assume that each retained record originated from one of its candidates. They exclude unknown origins outside the model, including unmatched records. They are not confidence intervals, do not certify biological truth, and their marginal maxima need not be simultaneously achievable. No flanking-sequence validation, indel model, base-quality weighting, cell calling or UMI deduplication is implemented.

## Scientific contribution and next publication gate

This establishes an executed read-conserving, joint-position audit with separate guide and gene resolution, not a new biological mechanism. The scientifically testable contribution is whether position-aware candidate accounting preserves useful gene-level information while exposing recycled guide evidence. Context-aware comparators, outcome-blind replicated independent screens and accuracy/calibration evidence remain necessary before claims of superior screening performance or a publication-ready general method.

Grouping ambiguous reads is established prior art; elementary marginal bounds are not claimed as new mathematics. The strongest potential paper combines an assay-specific mechanism, complete empirical accounting, validated implementation and downstream robustness evidence rather than merely observing that counting programs differ.

## Reproduction

Protocol commit: `de0002a37259cf27000015f43dd335343cd01234`. Full-audit code commit: `68a406b41805b808e6428927c2b8ccd2a2e55f9f`. [Full audit run](https://github.com/dnncha/dotmatch/actions/runs/34032669590). [Prior independent replay](https://github.com/dnncha/dotmatch/actions/runs/34030217143). Source commands and hashes are recorded in every completion manifest. Primary source accessions remain the original input provenance; GitHub transport artifacts can expire.
