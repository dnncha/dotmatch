# CRISPR Count First Run

This tutorial uses the tiny checked workflow fixtures under
`examples/workflows/fixtures/`. It does not download public data. The goal is to
show the production assay path, the MAGeCK-compatible count matrix, and the
sample QC table in a few commands.

## 1. Build DotMatch

```bash
make
```

## 2. Recommended: scaffold or start an assay project

The production path is `dotmatch assay new` followed by `dotmatch assay start`
(or `./run.sh` inside a scaffolded project). That runs preflight `check`, counts
guides, runs CRISPR QC, and writes a reliability report with suggested
`assay_fixes.tsv` edits when thresholds fail.

From the checked fixture:

```bash
cd examples/workflows/fixtures
../../dotmatch assay start crispr_assay.toml
```

To scaffold a fresh project from your own FASTQs:

```bash
dotmatch assay new crispr \
  --library guides.csv \
  --reads-dir fastqs/ \
  --out crispr_screen/

cd crispr_screen
./run.sh
```

Key outputs under the configured `out_dir`:

- `counts.mageck.tsv` — MAGeCK-style count matrix
- `sample_qc.tsv` — per-sample assignment and representation QC
- `summary.json` — run metadata and assignment rates
- `reliability_report.html` — evidence-bounded preflight/postrun review
- `crispr_qc.json` — guide-level QC summary

CPU remains the assignment authority. GPU Metal is opt-in via `[backend]` in the
assay spec and requires `--metal-validate` when enabled.

## 3. Direct `crispr-count` (single command)

For a minimal single command without the full assay wrapper:

```bash
mkdir -p tmp/crispr-first-run

cat > tmp/crispr-first-run/samples.tsv <<'EOF'
sample_id	fastq
sample_a	examples/workflows/fixtures/sample_a.fastq
sample_b	examples/workflows/fixtures/sample_b.fastq
EOF
```

`crispr-count --samples` accepts TSV or CSV sample sheets. With a header, the
sample column may be named `sample_id`, `sample`, or `label`; the FASTQ column
may be named `fastq`, `fastq_path`, `reads`, `path`, or `file`. Without a
header, DotMatch treats the first column as the sample name and the second as
the FASTQ path.

The fixture library contains three guides:

```bash
cat examples/workflows/fixtures/crispr_library.csv
```

Guide libraries may be TSV or CSV. DotMatch detects common CRISPR headers such
as `sgRNA`, `sgRNAID`, `guide_id`, `gRNA.sequence`, `sgRNA_sequence`,
`guide_seq`, `sequence`, `Gene`, and `gene_symbol`.

```bash
./dotmatch crispr-count \
  --library examples/workflows/fixtures/crispr_library.csv \
  --samples tmp/crispr-first-run/samples.tsv \
  --guide-start 0 \
  --guide-length 4 \
  --k 1 \
  --metric hamming \
  --ambiguity-policy radius \
  --out tmp/crispr-first-run/counts.mageck.tsv \
  --summary tmp/crispr-first-run/qc.json \
  --ambiguous discard
```

`sample_qc.tsv` is written automatically beside `--out`. Progress and QC review
warnings go to stderr on long runs.

## 4. Inspect the count matrix

```bash
cat tmp/crispr-first-run/counts.mageck.tsv
```

Expected output:

```text
sgRNA	Gene	sample_a	sample_b
guide_a	GENEA	0	0
guide_b	GENEB	0	0
guide_c	GENEC	0	1
```

The output is ready for MAGeCK-style downstream analysis: `sgRNA`, `Gene`, then
one integer count column per sample. DotMatch does not run MAGeCK statistics; it
only writes the count matrix expected by those tools.

## 5. Inspect sample QC

```bash
cat tmp/crispr-first-run/sample_qc.tsv
```

The key columns are:

- `total_reads`: input reads observed for the sample.
- `assigned_reads`: reads assigned uniquely to one guide.
- `exact_reads`: exact guide-window matches.
- `k1_rescued_reads`: one-edit rescued reads.
- `ambiguous_reads`: reads matching multiple guides within the allowed radius.
- `no_match_reads`: valid guide windows that matched no guide.
- `invalid_reads`: reads too short for the configured guide window.
- `assignment_rate`, `ambiguous_rate`, and `no_match_rate`: the same outcomes
  divided by valid extracted reads.
- `targets_observed`, `zero_count_targets`, `gini_index`, and
  `top_1pct_read_fraction`: guide-representation checks for quick review.
- `candidates_verified`: native target candidates checked after indexing.

In `sample_a`, the fixture deliberately includes one exact read that is withheld
because another guide is inside the one-edit radius, one ambiguous one-edit
read, one unmatched read, and one invalid short read. That is the behavior
DotMatch is designed to expose rather than hide.

## 6. Verify against the checked fixture outputs

```bash
diff -u examples/workflows/fixtures/expected_counts.mageck.tsv \
  tmp/crispr-first-run/counts.mageck.tsv
```

No diff means the tutorial count matrix matches the repository fixture.

For a public-data CRISPR example, use `examples/crispr_guides/run.sh` and the
checked evidence reports under `docs/benchmarks/public_crispr/`.