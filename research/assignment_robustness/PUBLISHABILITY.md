# Scientific impact: the claim worth earning

Decision note, 6 September 2026. This does not change AR003 endpoints or the completed audit code.

## The strongest question

**When apparent support from several CRISPR guide rows reuses the same sequencing records, what information is still identifiable, and does that reuse change replicated conclusions?**

This is more precise than comparing mapping rates or promising that a new counter discovers more genes. A gene's point estimate may be stable even when its apparently independent guide support is not. The distinction must be measured rather than assumed to inflate statistical significance.

AR003's specific contribution is a complete raw-read execution with joint target/position explanations and separately reconciled guide, gene-annotation and position resolution. A read compatible with several guides of one gene can retain gene-level information without pretending that it identifies one guide. Every retained record still has a visible one-count budget. The conditional ranges expose sensitivity to candidate allocation; they do not solve omitted-origin or error-model misspecification.

## Two deliverables, not one exaggerated claim

The near-term deliverable is a publicly reviewable technical research package: source, complete derived evidence, primary results, negative findings, known-origin controls, independent tests and executable reproduction. It should be cited as research code/evidence, not a peer-reviewed validation of a production default.

A stronger methods paper needs an independent replicated bulk-screen evaluation and matched assay-aware comparators. Its main result should explain *when and why* information is lost or duplicated, then test whether a conservative audit or correction improves a prespecified downstream endpoint. A biological paper would additionally need an independently corroborated gene/pathway finding; a newly crossing p-value is not that corroboration.

## Novelty exclusions from the literature review

- [ReCo, Bioinformatics 2023](https://doi.org/10.1093/bioinformatics/btad448) already automates CRISPR read counting with extraction/alignment context and supports single and combinatorial libraries. Staggered-position handling is not itself our invention. It is an appropriate comparator, not something to omit in favour of a weak fixed-window baseline.
- [CRISPR-Correct, official Pinello Lab repository](https://github.com/pinellolab/CRISPR-Correct) already maps imperfect protospacer observations with Hamming distance and supports additional construct components and several extraction specifications. Imperfect mapping or retaining ambiguity alone is not a novelty claim.
- [crispat, Bioinformatics 2024](https://doi.org/10.1093/bioinformatics/btae535) already demonstrates that guide-to-cell assignment strategies can change single-cell screen discoveries. AR003 is a different, bulk raw-read counting layer; it must not claim to originate the general assignment-to-conclusion question.
- [bcSeq, Bioinformatics 2018](https://doi.org/10.1093/bioinformatics/bty402) addresses sequencing-error-aware barcode counting. Accuracy comparisons cannot treat unweighted exact matching as a biological truth oracle.
- Ambiguous-read grouping and elementary per-candidate min/max allocation are established ideas. No new theorem or general-purpose invention is claimed for those components.

The literature search was focused, not exhaustive. The final manuscript requires a renewed review of current versions and closely overlapping methods before submitting a novelty claim.

## The decisive next experiment

Curate complete biological replicates with a compatible baseline, original guide reference, sample-level construct metadata and raw FASTQs. The [Sanson et al. primary study](https://doi.org/10.1038/s41467-018-07901-8) describes replicated screens and multiple CRISPR modalities, but the one Brunello plasmid archive in AR003 is not a replicated contrast. Different scaffolds and library versions must not be pooled or assigned the same baseline merely because they share the Brunello name. Also select at least one genuinely independent study/construct; its outcomes must remain unseen until the analysis plan is locked.

For each eligible contrast, compare matched-window policy arms first, then full native-workflow arms separately. Keep library annotation, normalization, filtering, contrasts and downstream implementation fixed within each comparison. Record read reuse, guide-support overlap, guide and gene representation, gene-effect estimates, rank stability, discovery overlap at a fixed threshold, and distance from that threshold. Biological replicates remain the inferential unit. Include downsampled/equal-depth analyses so a yield change is not mistaken for an assignment-specific benefit.

A prespecified known-origin error-channel evaluation should use complete real reference libraries and assay flanks, with balanced and imbalanced abundance and errors that violate the assumed radius. Retain wrong, ambiguous, unmatched and invalid outcomes. A narrow conditional range that excludes truth must count as a failure, not be hidden. Additional independent construct information or an orthogonal assay would provide stronger origin evidence than agreement between programs.

## Public communication rule

Lead with measured information preservation and transparent evidence accounting. Do not say that millions of reads are wrong, that published papers are invalid, that guide grouping is new, that all screens are vulnerable, or that improved mapping yield proves improved biological accuracy. Distinguish full-archive measurement, post-hoc mechanism witnesses, prespecified simulation, independent validation and unresolved hypotheses in every figure and abstract.

The technical package can invite narrowly scoped independent review now. Journal submission, author contact and a production behaviour change remain separate reviewed decisions. A null downstream result is retained and should narrow the paper rather than motivate selective reporting.
