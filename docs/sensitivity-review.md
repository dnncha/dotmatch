# Interactive assignment evidence review

**Unreleased source candidate.** This is an upgrade to the existing sensitivity
report, not a new matching engine, a clinical product, or an automatic policy
selector. Published DotMatch 0.5.0 reports retain the old static interface.

## Start with a scientific question

The opening view asks whether a different matching rule changes individual guide
counts, even when the assigned-read total is identical. Reference and comparison
policies are explicit. The overview, guide list and transition matrix all follow
the same selection. Switching policies never reruns or changes an analysis.

A guide opens in a keyboard-accessible inspector with all three recorded counts.
The viewer distinguishes net assigned-count change, differing guides, read-state
transitions and the global number of reads changed across any policy. It never
substitutes one of those quantities for another. A zero baseline has no invented
percentage increase.

Search includes the full supplied library, including unchanged and zero-count
targets. Filtering, sorting and pagination do not remove rows from the underlying
data. Exported TSV contains all matching guides, not just the displayed page.
Formula-like text identifiers are spreadsheet-escaped in this convenience export;
original count artifacts are unchanged. SVG figures contain both selected policy
outcomes and the first 12 guides in the current filtered order, with explicit
subset, scale and interpretation labels. Print output labels page scope.

## Generate or rebuild a review

From a built checkout of this candidate, the existing command automatically
writes the interactive `report.html` when within viewer capacity:

```sh
dotmatch sensitivity \
  --targets examples/assignment_sensitivity/targets.tsv \
  --reads examples/assignment_sensitivity/reads.fastq \
  --target-start 0 --target-length 20 \
  --sample-label synthetic_example --write-read-changes \
  --out-dir sensitivity-review-example
```

Open the HTML in a current desktop browser. No web server, account, installation
or network is required to review the generated document. JavaScript-disabled
browsers receive a labelled static snapshot; they do not receive a misleading
empty interactive result. The no-JavaScript table contains the first 50 guides.

Existing completed v1 bundles can be reviewed without rerunning the reads. This
source-only rebuilding route uses the Python standard library, not the native
engine, and never replaces a file or modifies the original bundle:

```sh
python python/dotmatch/sensitivity_review.py \
  --bundle path/to/completed-sensitivity --out new-review.html
```

Rebuilding requires `summary.json`, `guide_deltas.tsv` and `transitions.tsv`.
Their sizes and hashes must match the supplied completion manifest. The optional
read-change file need not accompany the portable report until it is attached.
A re-rendered report is a new artifact; it is not the original manifest's
`report.html`. The renderer identity is included in the new document.

## Read decisions, without invented evidence

Select a transition or guide, then optionally attach the run's original
`read_changes.tsv`. Native browser Web Crypto checks its SHA-256 against the
report's recorded artifact identity. The parser checks columns, record ordinals,
states, target IDs, Hamming policy-call consistency, changed-record count,
off-diagonal transitions and per-guide deltas before displaying any records.
A mismatch or unsupported hash context leaves no accepted attachment.

Read IDs may repeat; the 1-based record ordinal distinguishes occurrences. The
file contains only records changed across any policy, not every record in a
selected transition cell. At most 100 matching records are displayed, with the
actual matching count stated. Clearing, replacing with an invalid file, or
cancelling a pending attachment removes stale read-ID views, including a closed
inspector's DOM. The original HTML never embeds read IDs. Read-ID tables are
excluded from the viewer's TSV, SVG and print exports. Do not save the rendered
page after attaching private records.

Recorded-policy explanations describe consequences of the fixed matching rules.
They do **not** fabricate a base mismatch, candidate sequence, quality score or
biological conclusion. These artifacts do not contain raw reads, qualities or
enumerated candidate targets. A sequence-level inspector requires a separate,
explicitly consented evidence export; it is not simulated here.

## Evidence and privacy boundaries

The renderer reconciles all guide counts and all 48 transition cells against the
summary. Unsupported schemas, missing cells, duplicate targets and contradictory
statistics fail closed. Generated reports distinguish a producer's staged
snapshot from a rebuild checked against a supplied completion manifest. Normal
producer publication remains no-clobber, with `summary.json` published last.

A matching checksum proves consistency with the supplied identity, not that the
identity or report is authentic. Original input bytes are not re-read by the
viewer. There is no signature verification, authenticated human approval, sample
authentication, assay-validity score, clinical conclusion or claim of biological
accuracy. More uniquely assigned reads are not inherently a better result.

The report has no external dependencies, network requests, telemetry, browser
storage, external fonts or automatic policy selection. Its content-security
policy denies connections and permits only the generated script hash. Input
identifiers are inserted as text; JSON cannot terminate its script element.
The HTML, figures and exported count tables still contain sample/target
identifiers and may describe unpublished research. Review before sharing.

## Capacity and architecture

The portable view supports up to 250,000 guides, 64 MiB per required TSV and
16,384 characters per displayed identifier. Optional attached read decisions
are limited to 32 MiB / 100,000 records; the viewer does not silently sample a
larger file. The producer retains a labelled static fallback and all full count
artifacts when the interactive viewer exceeds its capacity. Corrupt evidence
is not treated as a capacity fallback.

`dotmatch.sensitivity.v1` remains the scientific output contract. The embedded
presentation snapshot is `dotmatch.review.v1`; it is not a new source of scientific
truth. Python reconciles/serializes; dependency-free JavaScript filters and renders
recorded data. Assets are Python resources so ordinary source/wheel packaging
includes them without adding a bundler or duplicating the matcher in JavaScript.
Large optional read-change artifacts are hashed once by the producer, not read
again merely to create an attachment binding.

This PR does not replace the Community Workbench or couple Sequence Doctor to
DotMatch. Workbench embedding and independent Sequence Doctor change-review
adapters remain separate integration work. Execution/approval controls are not
introduced into this read-only report.

## Validation and release gates

```sh
make shared
python -m pytest -q python/tests/test_sensitivity.py python/tests/test_sensitivity_review.py
python -m pip install pytest==9.0.2 playwright==1.57.0
python -m playwright install chromium
python scripts/check_sensitivity_review.py --file-mode --screenshots review-checks/screenshots --result review-checks/result.json
```

The dedicated `sensitivity-review` workflow builds the native engine, runs the
existing sensitivity tests and new renderer tests, then requires real local-file
Chromium navigation and native browser Web Crypto. Its result must be checked on
the actual candidate commit; configuration is not evidence of a successful run.

Local validation on 15 September 2026: **53 renderer/contract cases passed; three
native-integration cases were skipped** because the native package was not
available in the restricted working container. **26 Node checks** exercised the
actual JavaScript validators and native Node Web Crypto. **38 Chromium DOM checks**
passed; eight accepted-evidence DOM checks used an explicitly identified Python
digest test double because managed browser policy blocked local-file and localhost
navigation. That is not browser secure-context acceptance. No outgoing requests
or JavaScript exceptions were observed during those DOM checks.

A single synthetic 100,005-guide exercise produced a 4,160,141-byte HTML document:
0.822 s Python rendering, 1.215 s Chromium DOM load, 0.104 s showing all guides and
0.142 s searching the final guide (including debounce), in that container. This
is a recorded local fixture measurement, not a universal performance guarantee.

Fresh native/installed-wheel checks, complete repository CI, local-file browser
acceptance, and Safari/macOS acceptance remain release gates until actually run.
No merge, package release, deployment, customer pilot, scientific validation or
institutional security certification is implied by this implementation.
