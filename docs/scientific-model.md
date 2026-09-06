# Scientific model

EditWitness asks an identifiability question about a **specified measurement
model and finite set of local diploid states**. It does not infer how often an
editing outcome occurs or whether a submitted clone actually has that outcome.
The two shipped models are idealizations, not calibrated PCR simulators.

## 1. Local genomic states

A manifest supplies one forward-oriented local reference and named alleles.
An allele is a sorted, nonoverlapping list of replacements of that original
reference. `[start,end)` is zero-based and half-open; an empty replacement
sequence deletes that interval. An insertion has equal start and end.

A hypothesis contains exactly two allele IDs. Repeated IDs represent two copies.
The engine reconstructs the complete local sequence of each allele and compares
unordered pairs of those sequences **with multiplicity retained**. Different IDs,
allele ordering, or edit representations do not make sequence-identical pairs
into distinct local genomic alternatives. Such aliases remain in assessments,
but only one representative contributes a counterexample or planning constraint.

This is local sequence equivalence, not an assertion of genome-wide equivalence.
A whole-window deletion is represented by an empty local allele; its outside
breakpoints, distant rearrangements and chromosome state are not inferred.

## 2. Two explicit observation models

### `original-sites-presence-v1` — retained historical model

The original annotated primer intervals must remain pristine according to the
reference edit coordinates. An overlapping replacement disrupts the site even
when the replacement happens to preserve its bases. Boundary insertions have
the eligibility rules implemented and independently tested in `sequence.py`.
Only the product between those original sites is considered. No new, rescued,
relocated or alternative-orientation sites are searched.

This makes the model intentionally conservative about the original sites but
also **dependent on edit representation**. It is retained to make the original
assumption explicit and to support model comparison. Existing manifests that
omit `observation_model` continue to select it. It is not an exact model of which
PCR products a reconstructed sequence could yield.

### `exact-local-sites-presence-v2` — sequence-aware model

The reference intervals define two 5′→3′ oligos: the forward oligo is the left
interval sequence; the reverse oligo is the reverse complement of the right
interval sequence. On every reconstructed allele, the engine searches all exact
matches of these oligos and their required opposing binding sequences.

Both inward-facing heteroprimer configurations are considered: forward primer
on the left and reverse primer on the right, and the converse orientation. The
insert from a converse product is reverse-complemented so that all signals are
oriented forward-primer to reverse-primer. Products include both primer binding
sites for size filtering, but both primers are excluded from sequence readout.
Adjacent binding sites can produce an empty insert; an empty insert is a signal,
not an absence of product. Overlapping binding sites are not included.

Every in-bounds product contributes a signal, including products formed by
rescued or newly introduced sites. Product locations are expressed on the
**edited local allele**, not the input reference. Locations and full product
lengths are diagnostic metadata; they do not add unmeasured information to the
observational signature. Sequence-identical products collapse to the same signal.

This model still excludes primer mismatches, same-primer amplification, sites
outside the supplied window, secondary structure, competition, polymerase
failure, sequencing error and sampling. It assumes **every retained product is
detected**. In particular, an additional modeled signal can make an alternative
distinguishable even when a real experiment would fail to detect that signal.
That is a model-conditional distinction, never a claim of experimental sensitivity.

### Comparing assumptions

`editwitness compare-models experiment.json` runs the same hypotheses and assays
under both models and reports changes in witnesses and candidate panels. A
change demonstrates sensitivity to modeling assumptions, not validation of one
model against biology. New `demo` and `init` manifests explicitly select v2;
`demo --legacy-model` explicitly selects v1. No existing omission silently opts
into different scientific semantics.

## 3. What the assay observes

The readout of one retained product is either:

- `full_insert`: one complete primer-trimmed insert sequence; or
- `paired_end`: an ordered pair consisting of the first `read_bases` insert
  bases and the reverse complement of its last `read_bases` bases.

`read_bases` is usable insert sequence after trimming, not raw instrument cycles.
Reads may overlap. The paired-end signal does not include the unsequenced gap,
full insert length or product coordinates. Two products with different hidden
interiors can therefore have identical paired-end signals. Choosing a long-read
platform alone does not establish that a full insert was observed.

For assay a and allele x, let P(a,x) be the **set** of these sequence signals.
No eligible product means the empty set. A diploid hypothesis h=(x,y) has
S(a,h)=P(a,x) union P(a,y). This union deliberately discards allele multiplicity,
PCR yield, read fractions and copy number. Two identical copies and a single
contributing copy may therefore produce the same observed set.

Signal identifiers are SHA-256 hashes of canonical JSON read tuples. The actual
sequences are retained as evidence. These identifiers are for deterministic
comparison, not sequence anonymization or authentication.

## 4. Counterexamples and conclusion states

A distinct local genomic state is a counterexample when S(a,h) equals S(a,e)
for **every existing assay**, where e is the expected hypothesis. The output
includes its allele definitions and lists candidate assays whose modeled signal
sets differ from those of the expectation.

Four result states avoid treating lack of testing as evidence of success:

| Conclusion | Meaning |
|---|---|
| `no_distinct_alternatives` | The manifest contains no local genomic alternative distinct from the expectation. No discrimination claim is available. |
| `baseline_uninformative` | The expected state has no modeled signal in any existing assay. Absence of a product is not treated as reassuring validation. |
| `ambiguity_demonstrated` | At least one distinct declared alternative matches all baseline observations. |
| `distinguishable_only_within_declared_model` | Every distinct declared alternative differs under at least one baseline assay, conditional on all assumptions. |

The first two take precedence over the latter two. Witnesses may still be
reported when the baseline is uninformative; the state remains explicitly
qualified. A counterexample is **not evidence that it occurred**. No state
certifies safety, exhaustive coverage, experimental sensitivity or clone release.

## 5. Candidate-panel planning

For each candidate assay, the engine records which current counterexamples it
separates from the expectation. The objective is to cover all counterexamples
that any supplied candidate can separate, at minimum total declared integer
cost. It does not try to distinguish every alternative from every other one.
Alternatives no supplied candidate can separate remain unresolved in the result.

With at most 18 useful candidates, all subsets are searched. Ties favor fewer
assays, then lexicographic IDs. More useful candidates use deterministic weighted
set cover with exact rational gain/cost comparisons; optimality says `not_proven`.
This is conditional selection among **supplied assays**, not de novo oligo design
or proof of real-world assay performance. Costs have arbitrary user-defined units.

## 6. Bounded deletion hypotheses

`expand-deletions` converts a declared endpoint grid into explicit alternatives:
one expected allele plus a single deletion of the **reference haplotype**. It
does not superimpose a deletion on the intended edit. A homozygous local expected
sequence is required so that the program never silently chooses a phase.

Endpoint ranges are inclusive only at grid positions reached by `step`; invalid
or length-filtered deletions are excluded. Equivalent local diploid sequences
are collapsed against existing and generated hypotheses. Provenance records the
grid, source-manifest hash, valid-event count and duplicate count. Every unique
state must fit the declared limits; otherwise generation fails instead of
silently subsampling. A grid is not a distribution of repair outcomes.

The expanded manifest goes through the same readout-equivalence engine as
hand-authored hypotheses. It remains a finite, incomplete set of possibilities.
Genomic events outside that set remain untested.

## 7. Geometry scanner is separate

`scan` streams single reference deletions and counts disruption of the original
sites or exclusion by hard product bounds. It does not rematch exact sites,
compare read sequences, analyze hypotheses, or include candidate assays. Its
result always identifies `original-sites-presence-v1`, even when the manifest
selects v2 for the main analysis. Use `expand-deletions` followed by `analyze` for
sequence-aware hypothesis comparisons.

Changing the grid changes its denominator. Scan counts are never biological
frequencies, risk estimates, recall, precision, sensitivity or calibration.

## 8. Limits and missing measurements

The reference is at most 20 kb; a manifest supports 128 alleles, 1,000 diploid
hypotheses, 16 baseline assays and 24 candidates. Conservative sequence work is
bounded at 20 million bases. Exact matching additionally bounds each search at
512 matches, each allele/assay at 128 products, and total product reconstruction
at 20 million bases. Exceeding a bound fails the whole analysis, not part of it.

Hypothesis generation accepts at most 5,000 endpoint pairs and must fit the
manifest limits. Geometry scanning accepts at most 500,000 endpoint pairs.
These are computation limits, not statements about biological completeness.

An independent copy-number measurement might resolve a sequence-only ambiguity,
but no calibrated quantitative assay is implemented here. Such an extension
needs defined measurement error, controls, detection limits and independent
experimental assessment. An idealized perfect copy-number oracle must not be
presented as a validated laboratory method.


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
