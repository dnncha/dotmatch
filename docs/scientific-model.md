# Scientific model

## The question

For a finite set of local alleles, diploid clonal hypotheses, and configured
assays, which distinct genomic states cannot be distinguished from the expected
state by a specified observation function?

This is a deterministic design calculation, not inference from sequencing data
and not a probability model of DNA repair or assay performance.

## Coordinates and sequence identity

A reference is an explicit A/C/G/T sequence. Edits replace `reference[start:end]`
with `sequence`, in **original-reference, zero-based, half-open coordinates**.
Insertions use `start == end`; deletions use an empty replacement. Overlapping,
ambiguous and no-op edit declarations are rejected. Both expected and alternative
hypotheses contain exactly two declared alleles; other ploidies and mosaic
mixtures are unsupported.

In the new model, the final DNA is reconstructed before any primer matching.
Equivalent final sequences have the same observations regardless of edit
notation. A diploid state is an unordered pair of final sequences, with
multiplicity retained when comparing genomic states. An alternative with the
same pair as the expectation is explicitly labeled and excluded from witness
counts. The tool does not distinguish chromosome origin or haplotype context
outside the supplied local window.

## Exact local sequence model (v2)

Identifier: `exact-local-sequence-presence-v2`.

The reference specifies a forward oligo F and reverse oligo R, either explicitly
or derived from their annotated reference intervals. Each allele is reconstructed.
All exact F, R, reverse-complement(F) and reverse-complement(R) sites are found
inside that final local sequence. Two productive arrangements are enumerated:

```text
plus strand: F ... reverse_complement(R)   => forward orientation
plus strand: R ... reverse_complement(F)   => reverse orientation
```

Sites must be nonoverlapping and inward-facing. An adjacent pair with an empty
primer-trimmed insert is permitted. Product length includes both primer sites;
user-declared minimum/maximum bounds are hard exclusions. Surviving products are
normalized to F-to-R orientation. Multiple products are retained, not silently
reduced to the first match. The model includes new and recreated exact sites.

Identical F and R oligos are rejected because their read orientation is ambiguous
under this model. Products requiring the same primer at both ends (F/F or R/R)
are not represented. Mismatched binding, thermodynamics, competition, PCR cycle
count, sampling, limits of detection and off-window binding are also absent.
A lack of exact local product is **not** an empirical dropout diagnosis.

For a normalized primer-trimmed insert I, the measured sequence signal is:

```text
full_insert: (I,)
paired_end:  (first k available bases of I,
              reverse_complement(last k available bases of I))
```

Here k is the usable post-primer-trim read length. Short inserts use all available
bases. Empty observed inserts are signals, distinct from no product. Product
length and binding-site coordinates are diagnostic metadata, **not part of the
paired-end equivalence key**. The unsequenced gap must not leak into a decision.

All eligible exact local products are assumed detectable. A hypothesis observes
the union of its alleles' distinct sequence signals for each assay. Multiplicity,
fractions, counts and copy number are unobserved. Consequently, two copies of the
same intended allele and one intended allele plus an unobserved allele may
produce the same signal set. A predicted sequence difference need not be
detectable experimentally under these idealized assumptions.

## Original-site model (v1)

Identifier: `original-sites-presence-v1`.

This legacy model requires both annotated original primer intervals to survive
pristine and the reconstructed product to satisfy size bounds. Edits overlapping
a primer site disrupt it; an insertion at the boundary does not. It does not
search the edited DNA for rescued or new sites, so complex replacement notation
can affect its answer. This limitation is now a prominent machine-readable
notice. The model remains available for explicit legacy comparisons; it has not
been silently redefined.

Old manifests that omit `observation_model` continue to select v1. New `demo` and
`init` templates explicitly select v2 and schema 1.1. Use `compare-models` to
inspect sensitivity to these assumptions, not to decide which outcome occurred.

## What a witness proves

For every existing assay s, an alternative H is a witness when:

```text
observed_signal_set(H, s) == observed_signal_set(expected, s)
```

and its final local diploid sequence pair differs from the expected pair.

This establishes equivalence **within the declared model and hypotheses**. It
does not show that the alternative exists, estimate its likelihood, identify all
possible alternatives, or validate the expected state. With no witnesses, the
conclusion remains `distinguishable_only_within_declared_model`, never “safe.”

Every result records allele edit definitions and final sequence lengths/hashes,
including alleles with no signal. A focused `witness --include-sequences` command
also emits the final alternative DNA for inspection. Checksum integrity is not
an authenticated signature and does not establish scientific truth.

## Hypothesis generation versus geometry scanning

`expand-deletions` enumerates the valid single deletions in the declared reference
grid. Each is a deletion of the reference haplotype, paired with one selected
fixed expected allele. It is **not** a deletion applied on top of an intended
edited haplotype or a repair-outcome prediction. Physically identical local
states are deduplicated, keeping the first grid representation. Bounds, step,
input hash, filtering counts, and caps are recorded. Capacity excess is an error,
not subsampling. Generated results are then analyzed by the chosen full engine.

`scan` is different: it streams original-site deletion geometry and hard-size
eligibility counts without comparing read sequences or diploid states. It
always identifies itself as v1. Grid-dependent counts are not outcome frequencies,
risk, sensitivity or the proportion of all biologically possible deletions.

## Candidate panels

The objective separates the expectation from each currently equivalent
alternative—not all alternatives from each other. Candidate coverage is a
modeled difference from the expected observation. Alternatives no supplied
candidate covers remain explicitly unresolved.

Safe dominance filtering removes a candidate only when a no-worse deterministic
choice covers its entire witness set. With at most 18 remaining useful candidates,
all subsets are evaluated: minimize summed integer cost, then number of assays,
then lexicographic IDs. Above 18, a deterministic cost-aware greedy method is
labeled `not_proven`. None of this proves globally optimal experimental design,
experimental discriminability, or sufficient biological validation.

## Fail-closed resource limits

Exact matching permits at most 10,000 hits per primer sequence, 4,096 products
per allele/assay, 20,000 products and 20 million observed bases per analysis.
Manifest and generation limits also apply; inspect `capabilities`. Exceeding a
limit returns a structured error, not a truncated reassuring result. The HTML
may abbreviate previews without removing evidence from the full result.

## Scientific context, not package validation

Weisheit et al. describe on-target defects missed by routine genotyping and
quantitative/SNP-based validation. This motivates the problem but does not
validate this program or its assumptions:

- Weisheit I et al. *Detection of Deleterious On-Target Effects after HDR-Mediated
  CRISPR Editing.* Cell Reports 31(8), 107689 (2020).
  [DOI: 10.1016/j.celrep.2020.107689](https://doi.org/10.1016/j.celrep.2020.107689).
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/32460021/).
- [CleanFinder allelic-dropout analysis](https://cleanfinder.org/allelic_dropout)
  is relevant adjacent work. Allelic-dropout analysis itself is not new.
- [CRISPR-Analytics / CRISPR-A](https://pubmed.ncbi.nlm.nih.gov/37253059/)
  is relevant prior analysis/simulation work. Simulation and design support are
  not unique to EditWitness.

The release has no independently adjudicated experimental accuracy estimate.
Independent scientific review and a provenance-complete biological benchmark
remain explicit release-development gates.

### Haplotype choice during deletion generation

For a heterozygous expectation, the caller must explicitly choose the expected
allele to preserve. The other member of each generated challenge is a deletion
of the supplied reference haplotype, **not** a deletion composed onto the other
edited allele. To inspect both preserved-haplotype choices, run separate, clearly
identified generations. Equivalent allele-list orderings never choose different
haplotypes implicitly. The generated set is still finite and uncalibrated.
