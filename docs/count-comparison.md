# Compare your current counts before switching

You should not have to replace a working analysis pipeline to find out whether
DotMatch helps. Keep your current result, run DotMatch on the same reads and
library, and compare the two raw-count tables locally.

The useful question is not simply “do the totals agree?” Ten counts lost from
one guide and ten gained by another leave the total unchanged. This comparison
shows that difference rather than hiding it behind a correlation coefficient.
It does not tell you which result is biologically correct.

**Availability:** `dotmatch compare-counts` is included in DotMatch 0.6.0.
Install the release with `python3 -m pip install dotmatch==0.6.0`. The
`dotmatch-compare-counts` entrypoint and `python3 -m dotmatch.count_compare`
module are also available.

## Try a small, deliberately revealing example

From an installed DotMatch 0.6.0 package, run the checked synthetic first run:

```bash
dotmatch demo --out-dir first-run/
```

Open `first-run/comparison/report.html`. Exact and radius-one matching have
equal assigned totals but different guide counts. The source tables are in
`first-run/sensitivity/`. From a source checkout, the separate
`examples/count_comparison/` fixture exercises reordered rows and columns.

This fixture is synthetic. It tests the comparison mechanics, not a biological
claim, a speed advantage, or evidence of an independent user.

## Compare a real run

Use the **same biological samples and input reads**, the intended library
revision, and explicit extraction and preprocessing settings. A sample name
matching in two files does not prove that the underlying sample is the same.
Keep those records with the [lab evaluation and handoff](lab-evaluation.md).

```bash
dotmatch compare-counts \
  --baseline existing-workflow/counts.tsv \
  --candidate dotmatch-run/counts.mageck.tsv \
  --out-dir comparison/
```

Inputs are headered, tab-separated raw-count tables, optionally gzip-compressed
with a `.gz` suffix. MAGeCK-style `sgRNA`, `Gene`, and sample columns are
supported, as are DotMatch detailed outputs: only their `*_count_total`
columns become samples, not their exact/corrected components or QC fields.
The established count-table reader retains textual guide IDs such as `NA`
and `001`, and accepts exact nonnegative integer values, including `10.0`.

Do not supply normalized abundances, log counts, gene-level scores, or a
cell-by-feature matrix. The two tables must use the same counting unit and
represent the same biological samples. This tool does not aggregate genes,
correct cell barcodes, or deduplicate UMIs for you.

### Different sample column names (next release)

The source checkout now supports an explicit two-column, tab-separated mapping
for the same biological samples. This option is **not in the published 0.6.0
package**. Create `samples.tsv` with this header and one row for each renamed
candidate column:

```bash
printf 'baseline\tcandidate\ncontrol_day0\tL001\ntreated_day7\tL002\n' > samples.tsv
```

Then use `--sample-map samples.tsv` with the ordinary `--baseline`,
`--candidate`, and `--out-dir` arguments. Unmapped columns retain their
original names. Missing names, duplicates, and collisions are refused before
any output is written. The JSON report records the exact pairs and a hash of
the mapping file; HTML marks biological identity as a user assertion. Confirm
identity from your sample sheet first. Counts and source files are not rewritten.

## What you get

| File | Purpose |
| --- | --- |
| `report.html` | An offline, readable summary and the 20 largest count changes per sample. |
| `report.json` | Exact integer summaries, input hashes, identity checks, all exclusions, and bounded change previews. |
| `changes.tsv` | Every changed guide/sample cell, with baseline, candidate, and candidate-minus-baseline counts. |
| `manifest.json` | SHA-256 and byte size for each output; written last to mark completed output. |

Each sample shows compared totals, increased/decreased guides, guides becoming
zero or nonzero, and the sum of absolute count deltas. Equal totals with
changed guide counts are called out explicitly. The absolute-delta sum is
**not** a number of reassigned reads: aggregate count tables do not contain
the read identities needed to establish that.

Count arithmetic is exact; a one-count difference above `2**53` is not rounded
away. Missing, negative, fractional, non-finite, and duplicate data are errors,
not values to repair silently. The HTML has no scripts or remote resources.

## Missing guides are not zero counts

By default, differing guide or sample sets stop the comparison. Resolve a
library or sample-sheet mismatch first. For an intentionally partial comparison:

```bash
dotmatch compare-counts \
  --baseline baseline.tsv --candidate candidate.tsv \
  --shared-only --out-dir partial-comparison/
```

The report is marked **partial**, lists every excluded guide and sample in JSON,
and reports excluded count totals for each compared sample. It never invents
zero counts for an absent guide. No overlap means no comparison.

Identically named gene or sequence annotation fields are checked for conflicts
on shared guides. A conflict stops the comparison. Missing annotations and
differently named fields are reported as unverified, not assumed equivalent.
The tool does not infer that `gene_id` and `gene_symbol` use the same namespace.

## Use in a regression check

Add `--fail-on-difference` to write the report and return exit status 1 when
counts differ **or** guides/samples were excluded. Status 0 means the compared
counts agree with complete axis coverage; it does not certify scientific
accuracy. Input or filesystem errors return status 2.

Use a new output directory for every comparison. Existing directories and
files are never replaced. Validation happens before output creation. A disk
failure may leave a partial new directory; without a completed manifest,
do not treat it as a finished comparison.

## Decide what a difference means

Inspect the affected guides, not just a headline statistic. Check library
identity and sequence annotations, extraction window and orientation,
read-quality filtering, mismatch radius, ambiguity handling, and counting
unit. Keep downstream analysis settings fixed when testing whether changed
counts alter a conclusion. Do not choose whichever policy yields more counts
or more significant hits simply because the number is larger.

A useful evaluation ends with one of three conclusions: the results agree
under the tested conditions; the differences have an understood explanation;
or a discrepancy remains unresolved and blocks changing the workflow.

## Keep private data private

The command makes no network requests and sends no telemetry. It hashes and
parses the same temporary local copy of each input, then removes those copies.
It omits input paths and filenames from the report. **Reports still contain
guide identifiers, sample names, and counts.** Review them before sharing and
follow your institution's data-handling rules. Import TSV identifier columns
as text when opening them in spreadsheet software.

For optional feedback, use the existing pilot-feedback issue form with a
small synthetic or approved shareable example. A useful report says what you
expected, what happened, whether the result was useful, and what blocked the
next step. Never post private reads or unpublished study identifiers. A
synthetic demo, download, or maintainer-run test is not evidence of an
independent laboratory adopting DotMatch.
