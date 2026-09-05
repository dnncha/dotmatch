# Scientific model

## The question, precisely

Given a finite collection of local alleles, a finite collection of diploid clonal
hypotheses, and explicitly configured assays, can those assays distinguish the
expected hypothesis from each alternative?

Model identifier: `original-sites-presence-v1`.

This is a question about a declared observation function. It is not a probability
model of CRISPR repair or of real assay performance.

## 1. Alleles and coordinates

A local reference is an uppercase A/C/G/T sequence. An allele is the reference
with sorted, nonoverlapping replacements. A replacement substitutes
`reference[start:end]` with `sequence`. All positions are local, zero-based,
half-open coordinates. Insertions have `start == end`; deletions have an empty
replacement sequence.

All edits use the **original reference coordinates**, not the coordinates after
previous edits. Ambiguous overlaps, repeated insertion boundaries, and no-op
replacements are rejected. Alleles are sequence configurations, not empirical
claims that those configurations occurred.

## 2. Original-site eligibility

An assay specifies an inward-facing left and right primer interval. The current
model emits a product only when both original sites remain pristine and the
reconstructed product meets the declared size bounds.

A nonempty replacement overlapping a primer interval disrupts its original site.
An insertion strictly inside the interval also disrupts it. An insertion exactly
at a site's boundary does not disrupt the site. Surviving original bases are
mapped through edit-induced coordinate shifts before extracting the insert.

**This is not a thermodynamic prediction.** A partially mismatched primer may
still amplify. A replacement can recreate the same binding sequence or create a
new one. The current original-site model does not rematch edited sequences,
model mismatch tolerance, or rescue those sites. Consequently, conclusions can
depend on how a complex replacement is represented. Do not use those cases as
experimentally established dropout predictions. Sequence-aware rematching and
representation-invariance tests are a priority before broader claims.

Only annotated sites are used. A local duplicate-site warning is not a full
specificity screen: other orientations, sites outside the supplied window,
paralogs, pseudogenes, and nonspecific amplification remain unmodeled.

## 3. The observation function

For allele `a` and assay `s`, define `O(a,s)` as either no sequence signal or one
of the following idealized observations:

* `full_insert`: the complete primer-trimmed insert sequence.
* `paired_end`: the ordered pair `(first k insert bases, reverse_complement(last k
  insert bases))`, where `k` is the number of usable **post-primer-trim** bases.

If the insert is shorter than `k`, the complete available insert is used. Empty
inserts are valid sequence observations, distinct from an allele emitting no
signal. Adapter handling and real read errors are not modeled.

A full-insert declaration means the full insert is genuinely observed. Merely
choosing a long-read instrument is not sufficient. Conversely, product length
appears in diagnostic metadata but **is not part of a paired-end observation**.
The tool does not smuggle knowledge of the unsequenced gap into the conclusion.

For a two-allele hypothesis `H = (a,b)`, the assay observes the **set**:

```text
O(H,s) = O(a,s) ∪ O(b,s), omitting absent signals
```

Equal sequences collapse. Multiplicity, read fractions, molecular counts, and
copy number are deliberately unobserved. `(intended,intended)` and
`(intended,unobserved)` can therefore have the same observation.

Every eligible signal is assumed detectable. Sampling, allele competition,
stochastic failure, contamination, and limits of detection are outside this
model. Even predicted differences need not be experimentally detectable.

## 4. What a witness proves

For existing assays `S`, a witness is a declared alternative `H_alt` such that:

```text
for every s in S: O(H_alt,s) == O(H_expected,s)
```

It proves observational equivalence **within this model**, conditional on all its
assumptions. It does not establish that the alternative occurred, is plausible at
some estimated frequency, or is the only alternative. No likelihood or posterior
probability is calculated.

If no witness exists, the conclusion is
`distinguishable_only_within_declared_model`. That wording is deliberate: the
provided hypotheses are not an exhaustive description of biology.

There is no inference from observed experimental data in this release. It is a
preflight/design and model-inspection tool. Incorporating real results requires
a separately validated adapter and a definition of observation uncertainty.

## 5. Candidate panels

For each currently equivalent alternative, a candidate covers it when its
modeled observation differs from the expected hypothesis. The planner covers all
alternatives distinguishable by at least one supplied candidate, while explicitly
retaining those that none can distinguish.

At most 18 useful candidates: enumerate all subsets, minimize the integer sum of
user-declared costs, then the number of assays, then lexicographic IDs. Optimality
is proven only for this objective, these candidates, and this model.

Above 18: greedily maximize newly covered alternatives per cost unit, with exact
rational comparisons and deterministic ties. Output says `not_proven`. The tool
does not claim globally optimal assay design or that selecting primers alone
resolves every alternative.

The objective distinguishes the expectation from alternatives. It does **not**
aim to distinguish every alternative from every other alternative.

## 6. Deletion geometry scan

The scanner enumerates single deletions `[start,end)` on the supplied reference
using configured inclusive endpoint ranges and a step. Length filtering follows
enumeration. Upper bounds are included only when reached by the step.

For each existing assay it counts disrupted original sites, eligible products,
and products excluded by hard size bounds. It stores only a bounded number of
blind examples. It does not compare read sequences or the declared genotypes.

Changing the grid changes the denominator. Counts are therefore **not** editing
outcome frequencies, biological risk, sensitivity, or a calibration curve.
Candidate assays and compound edits are not part of this scanner.

## 7. Orthogonal measurements

A missing local sequence signal may require a different kind of measurement.
This alpha reports that need without pretending to simulate a validated
copy-number assay. A future assay type needs its own response model, failure
modes, calibration, controls, and versioned contract.

## Scientific context and adjacent tools

Weisheit et al. reported hidden monoallelic deletions and loss of heterozygosity
in edited human stem cells, and described quantitative genotyping PCR and SNP
controls. This motivates inspecting what sequence-only assays cannot establish;
it does **not** validate the EditWitness observation function.

- Weisheit I et al. *Detection of Deleterious On-Target Effects after HDR-Mediated
  CRISPR Editing.* Cell Reports 31(8), 107689 (2020).
  DOI: [10.1016/j.celrep.2020.107689](https://doi.org/10.1016/j.celrep.2020.107689).
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/32460021/).
- [CleanFinder's Allelic Dropout SNP Analyzer](https://cleanfinder.org/allelic_dropout)
  analyzes long-read FASTQ data for heterozygous SNP evidence. Allelic-dropout
  analysis itself is not new. EditWitness instead compares explicit assay-model
  counterexamples and candidate panels without replacing a caller.
- [CRISPR-Analytics / CRISPR-A](https://pubmed.ncbi.nlm.nih.gov/37253059/)
  is relevant prior work in genome-editing analysis and simulation. Simulation
  and experiment-design support are not unique to this project.

References establish context, not package accuracy, experimental validation,
exhaustive novelty, or endorsement. A formal competitive review and independent
biological assessment remain release gates.
