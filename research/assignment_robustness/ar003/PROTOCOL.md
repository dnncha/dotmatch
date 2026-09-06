# AR003: full-read joint-position resolution audit

Status: prospective protocol amendment recorded before executing the new joint decoder on the complete archives in this amendment. This is a public versioned protocol, not an externally registered study.

## Prior knowledge and provenance

AR001/AR002 simulations and Yusa count-ratio analyses were already inspected. The independently completed replay at c33ff64a4a4cf931cb7b26e6ad3a17afe1cf9580 (workflow 34030217143) was inspected before this amendment. Its Yusa read reuse, Brunello plasmid comparison, selected positions, discovery-flank models and synthetic witnesses are known results, not held-out discoveries. The separate standalone continuation at e80619f2c963cdcbf12b998e3307a9a83388974f provides the tested joint-position semantics to validate here.

## Question

How much guide-level and gene-annotation-level information remains when all permitted positions are considered jointly and a sequencing record can contribute at most once to each resolution's table? Does merely resolving at gene level preserve information without pretending to identify an individual guide?

## Locked inputs and positions

Use full original compressed FASTQs and original reference tables from the successful archive acquisitions, verifying their archived SHA-256, ENA MD5, byte counts and record totals again. Transport ZIPs must match GitHub's published SHA-256. No synthetic fallback and no reconstitution of public reads from count tables.

Primary input pair: ERR376998 (10,093,905 records; plasmid) and ERR376999 (10,300,758 records; cellular ESC1), original Yusa 87,437-guide, 19-base reference with SHA-256 252e3b81b809c50f5cc347238a52926818027ad78a3ec98686e8012a8a46a896. Zero-based permitted starts [21,22,23,24], already selected in the earlier discovery-only replay; fixed-window comparison start 23.

Secondary descriptive design: SRR8297997 Brunello plasmid, 9,821,128 records, 20-base reference. Use exactly the archived discovery model's permitted positions and reference digest. Freeze those bytes before executing. This is not a held-out study or replicated phenotype comparison.

## Policies and outcomes

Enumerate every (position, target, Hamming distance <=1) explanation. Preserve duplicate target sequences as distinct IDs. Exact retains distance-zero explanations. Radius-one retains all explanations. Best-distance retains globally minimum-distance explanations, not separate votes from each window. Resolve guide identity and gene annotation independently. Multiple positions supporting the same target are not additional counts. Require all declared windows to be extractable; incomplete joint searches are invalid. Treat DNA/IUPAC symbols literally; no indels, Phred weighting, orientation inference or flank validation.

Primary descriptive outcomes for every sample/policy: complete unique/ambiguous/unmatched/invalid record budgets at guide, gene and position resolutions; guide-unique and gene-unique count tables; additional gene-identifiable records beyond guide-identifiable records; complete candidate-class sufficient statistics; marginal mandatory-allocation gene-count lower/upper bounds. Lower gene count must equal gene-unique count. Do not add guide-level and gene-level totals together.

Primary Yusa comparison: input-read-normalized cellular/plasmid log2 ratios with pseudocount 0.5; eligibility uses fixed-window exact plasmid gene sum >=50 and at least two guides with >=10 counts, identical across policies. Report number and proportion of eligible genes with absolute change >=0.5 log2 relative to fixed-window best-distance; all outliers retained. Report conditional effect ranges separately, not as confidence intervals or calibrated FDR. Threshold crossings are not biological discoveries. Alternative pseudocounts 0.1 and 1 and baseline thresholds 20 and 100 are labelled sensitivity analyses. No hit calling on this unreplicated pair.

## Validation

Run the original 127-test suite before changes. Validate optimized/cached execution against the uncached joint decoder and exhaustive all-target all-position enumeration on constructed panels, including ties, duplicates, literal N, missing windows, same-gene/cross-gene ambiguity and target-order invariance. Select independent real-record checks using fixed-seed random ordinals, not the most favourable outcomes. Reconcile per-position counts against pinned DotMatch where executed. Clearly separate sample checks from complete-count checks.

Also execute the decoder on the previously archived balanced error-free known-origin controls against each complete reference. These constructed controls test correct/incorrect/ambiguous guide and gene calls under the stated model; they are not new biological samples. Source origin outside the candidate set remains an explicit failure case. Repeat scientific outputs from a clean checkout or an independent implementation and report the exact scope; mere file hashing is integrity, not independent accuracy validation.

## Interpretation and publication

Bounds assume each matched record's true origin is among the permitted explanations and contributes exactly once. They do not cover missing references, errors outside the Hamming model, wrong positions or contamination. Marginal upper bounds need not be jointly attainable. Reads are not deduplicated molecules or biological replicates.

Ship reviewed research source, tests, complete aggregate outputs, input/source hashes and a clearly labelled technical report to a dedicated research branch/PR. Do not silently change DotMatch production behaviour, claim general accuracy superiority, claim new biological mechanisms, submit a manuscript, or contact authors without review. Maintain the external replicated-screen, comparator and biological-validation publication gates. Negative outcomes remain in the record.
