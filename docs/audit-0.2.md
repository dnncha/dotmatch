# Audit of 0.1.0a1 and changes in 0.2.0a1

This is an engineering and model audit by the implementing assistant, not an independent scientific review. The original source archive and live staging branch were checked before changes. Test counts refer to collected tests; parameterized cases are not independent biological experiments.

## Findings and disposition

| Finding | Why it matters | Change |
|---|---|---|
| Original-site eligibility depends on edit representation. | A broad replacement can preserve an exact primer sequence but be labeled disrupted. A newly introduced binding site can be missed. | Added an explicitly selected exact local rematching model, with both inward heteroprimer orientations and bounded multiple-product enumeration. Retained the historical model for sensitivity analysis. |
| Genotype aliases could be counted as counterexamples. | Renaming an allele or describing the same sequence differently can inflate an apparent blind spot and panel coverage. | Group unordered diploid reconstructed sequences, exclude the expected genotype's aliases, and emit one representative witness per alternative genotype. |
| No alternatives looked like successful discrimination. | A zero-witness result could mean no scientific comparison was made. | New `no_distinct_alternatives` conclusion. Baselines with no expected positive signal use `baseline_uninformative`. |
| Public Python entry points trusted constructed model instances. | Pydantic's unchecked copy/construction helpers can bypass validation. | Analysis, scan, generation and model comparison revalidate input at their public boundary. |
| The original report did not carry the complete edit definitions. | A witness could be harder to inspect without locating the original manifest. | Full results include declared allele edits; reports show genomic changes; explicit sequence export reconstructs the local alleles. |
| Manual hypothesis authoring limited practical use. | The scientifically interesting part required repetitive hand-written alternatives. | Added bounded deletion-grid generation with explicit provenance, duplicate handling, collision checks and hard refusal rather than truncation. |
| Publication and CI configuration were conflated with distribution. | The previous package was staged, not independently published. | Public source/install instructions explicitly identify the hosting branch. Build/release status records observed outcomes only; standalone GitHub and PyPI remain separate operations. |

## Independent software checks

The new exact model is compared with a deliberately slow oracle that applies replacements using reverse list slicing and searches every substring pair. It does not reuse production sequence reconstruction, site search, interval filtering or coordinate mapping. Seeded cases vary replacements, insert lengths, orientations, product bounds and readouts.

Focused regressions cover representation invariance, an introduced primer site, reverse-orientation products, multiple signals, empty inserts, paired-end gaps, alias hypotheses, resource refusals, API revalidation, model disagreement, and generation provenance. Existing original-site and panel-search oracle tests remain in place.

These checks establish implementation behavior under the stated models. They do not establish experimental sensitivity or specificity.

## Important residual limitations

Exact matching is still not a PCR simulator. Partial matches may amplify, exact products may not be sampled, primer competition may alter readouts, and nonlocal or same-primer products may matter. The model comparison exposes only two selected assumptions; agreement between them is not a calibrated confidence level.

Generated deletions alter the reference haplotype and are paired with an expected allele. They are not deletions superimposed on the edited haplotype, an exhaustive repair-outcome model, or a plausible-frequency prior. The grid and phase restriction are explicit.

A checksum demonstrates that content matches a recorded digest, not authentic authorship. Replay is supported with the producing package/model version; old evidence is not silently reinterpreted by a new schema.

## Next scientific gate

Obtain a publicly redistributable case with primer geometry, reference build and locus, raw or adjudicated observations, and orthogonal genotype evidence. Freeze the expected result before running the package. Seek review of the observation model from a genome-engineering specialist. Do not mark this complete with another synthetic fixture.
