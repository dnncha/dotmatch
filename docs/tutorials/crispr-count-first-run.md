# Your first CRISPR guide-counting run

Start with a checked example, then prepare your own reads. Keep the counting
workflow you already use while you evaluate DotMatch. No study data needs to
leave your computer.

## Check which version you are using

The released package can be installed with:

```bash
python3 -m pip install dotmatch==0.6.0
dotmatch --version
```

DotMatch 0.6.0 includes the checked offline demo, `dotmatch compare-counts`,
and the validated `crispr quickstart --link-reads` option. Use the published
package and confirm `dotmatch --version` reports 0.6.0 before following this
tutorial.

For a no-install view using the public synthetic fixture, open the
[assignment review example](https://dnncha.github.io/dotmatch/assignment-sensitivity/).

## 1. Get a checked result without assembling files

From an installed DotMatch 0.6.0 package:

```bash
dotmatch demo --out-dir first-run/
```

Open `first-run/index.html`. The package includes the tiny synthetic inputs;
the command does not fetch them or upload anything. It runs the existing native
assignment engine, checks every guide count and the four read-outcome totals
against the committed expectations, and creates an offline review bundle.

The example deliberately includes a near-neighbour guide pair, duplicate target
sequences, a short read, an unmatched read, and a literal `N`. Exact and
radius-one assignment each count three reads, but their counts differ for three
guides. Best-distance counts five reads. More assigned reads are not proof of
better biological accuracy. This is a software example, not a biological study.

The start page links to the interactive assignment review, the count-table
comparison, the complete tables and checksums. A completed `manifest.json`
marks a successful demo. If verification fails, the command returns status 2
and does not write that completion marker. Existing output paths are refused.

## 2. Prepare your own guide library and FASTQs

A guide table must identify the guide, its sequence, and, for a useful
MAGeCK handoff, its gene annotation. For example:

```text
target_id	sequence	gene
guide_a	GACTAGCTACGATCGTACGA	GENE_A
guide_b	TTGGCACCTTAGGACCTGAA	GENE_B
```

Use a TSV or CSV exported from the intended library revision. Do not replace
missing annotations by guessing. Retain the library's source and checksum.

Choose the read that contains the guide sequence. **Do not give quickstart both
R1 and R2 just because they are in the same directory.** Filename-derived sample
names are not biological replicate declarations, and technical lanes are not
silently combined into one biological sample.

```bash
dotmatch crispr quickstart \
  --library guides.csv \
  --fastq 'fastqs/*_R1*.fastq.gz' \
  --out crispr-screen/
```

Use paths that match your actual naming convention. Repeat `--fastq` for more
input patterns. DotMatch 0.6.0 requires every argument to resolve: one valid
file does not excuse a second missing pattern. Duplicate inputs, colliding
basenames, non-files and invalid resource limits are rejected before a project
is created. Existing output paths are never overwritten.

The default copies reads into the project. For large FASTQs, the evaluation
build supports:

```bash
dotmatch crispr quickstart \
  --library guides.csv \
  --fastq 'fastqs/*_R1*.fastq.gz' \
  --link-reads --out crispr-screen/
```

Linked reads stay in their original location, avoiding another full-size copy.
**They must remain available and unchanged. The project is not self-contained.**
The library is still copied, and `inference_report.json` records the original
read paths and whether reads were linked or copied.

The established directory-based route also supports linking:

```bash
dotmatch assay new crispr \
  --library guides.csv --reads-dir guide-fastqs/ \
  --link-reads --out crispr-screen/
```

Put only the intended guide-read inputs in that directory. Inspect the generated
status and settings before starting a run; do not assume this older command
uses the new quickstart's always-draft default.

## 3. Review the inferred settings before running

Open `crispr-screen/inference_report.json`, `samples.generated.tsv`, and
`assay.toml`. Check the sample mapping, library revision, zero-based guide
window, read orientation, matching rule and warnings against the assay protocol.
Inference proposes a configuration; it does not establish biological validity.

The quickstart stays in draft by default. After confirming the settings, change
the top-level `status = "draft"` to `status = "ready"` in `assay.toml`, then run:

```bash
dotmatch assay start crispr-screen/assay.toml
```

DotMatch 0.6.0's explicit `--accept-inference` option only starts a run
when inference itself reports ready. It no longer promotes an uncertain
inference automatically. `--no-run` always leaves the project in draft.

A setup failure can leave diagnostic files in the newly created project. Do
not treat an incomplete scaffold as ready; inspect the error and use a new
output directory for the corrected attempt. Existing projects are not deleted.

## 4. Review the outputs, not only the assigned percentage

Start with `assay_out/reliability_report.html`, then inspect the sample QC,
CRISPR QC and `counts.mageck.tsv`. A software QC verdict is not a substitute for
sample identity, positive and negative controls, or the study's analysis plan.

The matrix contains `sgRNA`, `Gene`, then one raw integer count column per
sample. Only unique assignments add target counts. Ambiguous, unmatched and
invalid-window reads remain visible in the QC outputs.

For a side-by-side evaluation with DotMatch 0.6.0:

```bash
dotmatch compare-counts \
  --baseline existing-workflow/counts.tsv \
  --candidate crispr-screen/assay_out/counts.mageck.tsv \
  --out-dir comparison/
```

Both tables must refer to the same samples, guides and counting unit. Open
`comparison/report.html`. Equal totals do not guarantee equal guide counts,
and equal guide counts do not prove biological accuracy. Do not select a policy
just because it produces more counts or significant hits. See
[Count comparison](../count-comparison.md) for exclusions and interpretation.

## 5. Keep a reviewable record

After reviewing a completed run:

```bash
dotmatch assay handoff crispr-screen/assay.toml
```

This packages the configuration, counts, reports, methods, citation and checksums
without copying raw FASTQs. Reports can still contain private sample names,
guide identifiers and counts. Review before sharing. Use the
[lab-evaluation guide](../lab-evaluation.md) to record the technical decision.

Useful optional feedback is specific: what task you tried, what became useful,
what blocked you, and whether you chose to use DotMatch again. The
[evaluation form](https://github.com/dnncha/dotmatch/issues/new?template=pilot_feedback.yml)
is public; do not attach private reads or unpublished study identifiers.

## Direct counting with an explicit sample sheet

The direct command also supports an explicit sample sheet:

```bash
dotmatch crispr-count \
  --library guides.csv --samples samples.tsv \
  --guide-start 0 --guide-length 20 \
  --k 1 --metric hamming --ambiguity-policy radius \
  --out counts.mageck.tsv --summary summary.json \
  --sample-qc sample_qc.tsv --ambiguous discard
```

These extraction and assignment values are examples, not recommendations for
your assay. Supply your confirmed settings. Sample sheets use `sample_id` and
`fastq` columns; absolute FASTQ paths avoid ambiguity about the working directory.
Keep biological sample identities explicit and do not treat separate lanes as
independent biological replicates. DotMatch writes MAGeCK-compatible counts;
it does not perform downstream screen statistics.
