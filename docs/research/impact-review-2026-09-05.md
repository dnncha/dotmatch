# DotMatch: scientific value, discovery, and adoption review

Date: 5 September 2026. Baseline: main commit ce5ade3dcae4e095d3a4695e12d7acb5a0f8cc70. This is a project review and proposed research programme, not an independent biological validation study.

## Decision

Make DotMatch the inspectable assignment and correction-risk layer between a known sequence library and downstream analysis. Lead with CRISPR guide counting. Extend through barcode QC and guide-capture workflows only where those workflows share the same assignment problem.

Do not make an autonomous agent, a new assay language, GPU support, or an all-in-one CRISPR platform the first explanation of the product. Those are implementation routes or research directions. The user-facing question is: which reads contributed to this count table, which did not, and how sensitive is the result to the assignment choices?

## What already exists

The project already includes a deterministic matcher, CRISPR count output, barcode demultiplexing, pair counting, feature matrices from extracted cell-labelled observations, library audits, diagnostic workflows, assay configuration, structured agent tools, public fixtures, package distribution, and citation metadata. It does not need another generic wrapper before those capabilities become understandable and useful.

Sources: [current README](https://github.com/dnncha/dotmatch/blob/ce5ade3dcae4e095d3a4695e12d7acb5a0f8cc70/README.md), [documentation](https://dotmatch.readthedocs.io/en/latest/), [release v0.4.1](https://github.com/dnncha/dotmatch/releases/tag/v0.4.1).

## Visibility findings and uncertainty

Searches surfaced DotMatch documentation and package entries, alongside the unrelated dotmatch component of GenomeMatcher. That older name appears in a 2008 scientific paper and also in CRISPR-related sequence-comparison contexts. This is a disambiguation problem, not proof that the current project is absent from every index.

The prior homepage led with a generic matching statement and an agent-first call to action. Its sitemap contained only the homepage. Existing structured data, llms files, documentation, package listings, and repository topics mean the answer is not simply to add more metadata files.

The connected Search Console application did not expose a DotMatch property. No claim is made here about its Google impression totals, exact positions, Google-selected canonical, crawl errors, or indexing exclusion reasons. Those need a verified property and URL inspection. A failed research-browser fetch is not evidence of a production outage.

Retain the name but pair it consistently with CRISPR guide counting and barcode QC. Keep one authoritative product identity across the homepage, GitHub, package metadata, documentation, and citations. Do not change domains without a redirect and canonical migration plan. Prioritise useful task pages with working tools and reproducible examples over a large set of thin keyword pages.

Google explicitly states that its AI search features do not need special AI files or special schema. Crawlability, indexability, useful textual content, internal links, and structured data that matches visible claims still matter. No ranking improvement is guaranteed by this patch.

Sources: [GenomeMatcher paper](https://doi.org/10.1186/1471-2105-9-376), [Google AI features guidance](https://developers.google.com/search/docs/appearance/ai-features).

## First implementation

The website now leads with the scientist's input/output task. It offers a dedicated CRISPR workflow and a local browser library checker before agent automation. Existing agent tools remain accessible lower on the homepage. Important existing section anchors are retained.

New task URLs have independent metadata and canonicals, internal navigation, sitemap entries, and static directory exports for GitHub Pages. Build time is no longer represented as an invented content modification time. Citation links use the stable concept DOI rather than implying a newly confirmed version-specific archive.

The library checker performs exhaustive Hamming radius-one neighbourhood ownership accounting over all supplied rows. It identifies duplicate target sequences, exact-read ambiguity under radius matching, and ambiguous single-substitution observations. It preserves target IDs, counts shared observations once, provides witnesses and an exportable report, and rejects unsupported input instead of silently sampling it. The interface is bounded to 2,000 targets and 8–32 ACGT bases. This is deliberately not a whole genome-scale browser pipeline.

The checker is a standalone transparent implementation, not the native matcher, a clinical safety tool, an off-target cleavage prediction, or an error-rate estimator. It does not collect sequencing reads. Only the displayed rows and examples are bounded; the aggregate audit uses every accepted target and JSON includes all per-target results.

## Scientific priority 1: assignment sensitivity and disagreement

The existing CRISPR comparison report is valuable but must not be interpreted as universal equivalence. Its Yusa exact-count comparison matches MAGeCK counts. Its Brunello exact-count row is marked non-comparable. Hamming comparisons marked ok still show non-zero count differences. A reported full-Brunello speed ratio is a one-run result, separate from repeated subsamples. These are different levels of evidence.

Add a separate machine-readable comparison interpretation rather than silently changing the existing status contract: execution_status, semantic_comparability, counts_identical, differing_guides, total_count_delta, validation_kind, and validation_scope. Regenerate reports from raw artifacts. A numeric total delta of zero also does not prove per-guide equality.

The strongest next native feature is a sensitivity ledger. Starting from a pinned baseline, compare exact versus one-mismatch assignment, radius versus best-distance ambiguity, fixed versus discovered offsets, and supported indel policies. For each transition, record which read or target changed, why it changed, whether it was newly counted or reassigned, and whether the change concentrates in a small set of guides. Preserve baseline counts; never silently replace them with a more permissive result.

Acceptance criteria: all input reads accounted for once under each declared counting policy; explicit unmatched/invalid states; fixed definitions and versions; bounded deterministic output; independent reference tests; discrepancy fixtures; and count-table conservation checks. A user should be able to see why two tools disagree without first learning the matcher internals.

Source: [baseline comparison report](https://github.com/dnncha/dotmatch/blob/ce5ade3dcae4e095d3a4695e12d7acb5a0f8cc70/docs/benchmarks/crispr_comparison/README.md).

## Scientific priority 2: independent, decision-relevant validation

Use three distinct evidence layers. Mathematical or exhaustive reference tests establish correctness under specified rules. Public real reads establish reproducibility and practical behaviour, but usually do not reveal the true origin of every read. Truth-labelled spike-ins or independently designed controlled experiments can estimate false assignments and missed assignments under that experimental setting.

Predeclare the task, offsets, reference library, ambiguity rules, strata, and evaluation metrics. Keep discovery and evaluation data disjoint. Include exact matches, substitutions, indels, duplicates, close neighbours, Ns, short reads, contamination, and targets absent from the whitelist. Include multiple library sizes and abundance distributions. Do not choose settings on the same held-out data used for the accuracy claim.

Report precision/recall only where origin is known; include abstention and the uncertainty of those estimates. Measure assignment changes by target and, with a fixed downstream method, whether those changes alter screen conclusions. Separate compressed FASTQ I/O, startup, indexing, steady-state assignment, and end-to-end runtime. Run repeated comparisons at matched resources and semantics. Publish checksums, commands, software versions, raw outcomes, and unsuccessful cases.

Do not use agreement with a reference matcher as biological ground truth. Do not interpret correlation as equality, assignment rate as accuracy, package downloads as unique users, or a decoy-derived error estimate as valid without its null-model assumptions.

## Scientific priority 3: selective, quality-aware decoding

The existing experimental calibration and assay layers are a starting point, not an excuse to introduce another overlapping command family. A defensible extension could use read quality and assay-specific error evidence to rank candidate assignments while retaining an abstain state and deterministic baseline outputs.

Train or calibrate against independently known origin, not against confident assignments made by the same system. Check calibration on held-out libraries and instruments. Ensure abundance priors do not absorb real low-abundance guides into common ones. Evaluate distribution shifts, absent targets, and imperfect quality calibration. Keep the experimental policy opt-in until false-assignment control is demonstrated. Do not advertise a confidence percentage before it has an empirical interpretation.

## Scientific priority 4: paired constructs and single-cell boundaries

Extend existing pair-count and assay infrastructure for combinatorial constructs only when expected linkage and read structure are specified. Distinguish independently assigned halves from a jointly supported construct. Report unexpected combinations; do not automatically equate them with genuine biological pairs or with a particular molecular artefact.

For Perturb-seq, preserve a clear boundary between sequence assignment and cell-level inference. Cell-barcode correction, UMI processing, cell calling, ambient-guide filtering, and perturbation-effect statistics are distinct tasks. Prefer tested adapters into established processing and analysis workflows over pretending a per-read counter has solved all of them.

CLEANSER provides a useful primary example: ambient gRNA handling can change guide-cell assignments and downstream differential-expression results, and its authors describe platform and training-data limits. Any DotMatch claim at that level needs its own corresponding evidence, not a per-read oracle test.

Source: [CLEANSER, Cell Genomics](https://doi.org/10.1016/j.xgen.2025.100766).

## Accessibility and ecosystem strategy

Support three entry routes: a no-install educational/diagnostic tool; a copyable, reviewed local workflow for an analyst; and a pinned workflow component for a core facility or pipeline maintainer. Keep data local by default. Explain failures with actual input or inference problems, not opaque generic errors. Provide concrete output examples and a methods/handoff artifact that can be reviewed without installing the application.

Prioritise one real upstream integration rather than many unaccepted wrappers. nf-core/crisprseq documents a count-table entry route, including sgRNA and gene columns. Validate exact header case, sample identity, gene identifiers, version, and normalization semantics in a tested adapter before calling it plug-and-play. A checked output adapter is not the same thing as an accepted nf-core module. MultiQC, Galaxy, and Snakemake are further candidate integrations; local examples alone are not proof of their acceptance.

Source: [nf-core/crisprseq screening](https://nf-co.re/crisprseq/usage/screening/). Comparator: [official guide-counter](https://github.com/fulcrumgenomics/guide-counter).

## Milestones, not popularity claims

First milestone: an external researcher reaches a useful output without maintainer intervention; record the actual friction. Proposed target: five completed external evaluations across at least three labs or teams, with their permission for any public attribution. These numbers are goals, not existing users.

Second milestone: three independently reproducible case studies, including one meaningful failure or disagreement analysis. Produce a compact methods figure or report grounded in those artifacts. Ask for concrete comparison feedback rather than generic endorsement.

Third milestone: one accepted upstream integration or documented repeated operational use by an independent core. Track repeat workflow use, accepted contributions, independent reproductions, citations, and non-branded search queries. Do not optimise on stars or aggregate package retrievals alone.

Publish focused task material only when it answers a real user problem: CRISPR counting; barcode mismatch collisions; high unmatched rate; tool disagreement; and the boundary between guide reads and guide-per-cell calls. Add each page to the sitemap and internal navigation. Verify the Search Console property, submit the deployed sitemap, inspect representative URLs, and measure changes after publication. Do not submit URLs that are only present on an unmerged branch.

## Merge and validation gates

The standalone audit has local reference tests covering exhaustive eight-base observations, seeded panels, duplicates, malformed inputs, multi-target overlaps, deterministic output, and supported input bounds. The pure TypeScript library was checked with strict TypeScript compilation. Metadata tests cover root hosting and GitHub Pages subpaths. A separate CI workflow builds the actual static export and checks rendered canonical URLs, H1s, links, and sitemap entries.

Before merging: run npm ci, npm run lint, npm run check:site, and a GitHub Pages-mode static export followed by node scripts/check_site_export.mjs. Run existing native/Python/evidence gates unchanged. Browser-review 390px, 768px, and 1440px widths, keyboard focus, table scrolling, error cases, cancellation, edits during a run, and JSON export; confirm no sequence-bearing network requests. Then verify deployed pages and Search Console state. This review did not run a new full biological benchmark or establish production ranking changes.
