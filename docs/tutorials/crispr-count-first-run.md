# CRISPR Count First Run

This tutorial uses the tiny checked workflow fixtures under
`examples/workflows/fixtures/`. It does not download public data. The goal is to
show the full CRISPR-facing command, the MAGeCK-compatible count matrix, and the
sample QC table in a few commands.

## 1. Build DotMatch

```bash
make
```

## 2. Create a sample sheet

```bash
mkdir -p tmp/crispr-first-run

cat > tmp/crispr-first-run/samples.tsv <<'EOF'
sample_id	fastq
sample_a	examples/workflows/fixtures/sample_a.fastq
sample_b	examples/workflows/fixtures/sample_b.fastq
EOF
```

The fixture library contains three guides:

```bash
cat examples/workflows/fixtures/crispr_library.csv
```

## 3. Count guides

```bash
./dotmatch crispr-count \
  --library examples/workflows/fixtures/crispr_library.csv \
  --samples tmp/crispr-first-run/samples.tsv \
  --guide-start 0 \
  --guide-length 4 \
  --k 1 \
  --metric hamming \
  --out tmp/crispr-first-run/counts.mageck.tsv \
  --summary tmp/crispr-first-run/qc.json \
  --sample-qc tmp/crispr-first-run/sample_qc.tsv \
  --ambiguous discard
```

## 4. Inspect the count matrix

```bash
cat tmp/crispr-first-run/counts.mageck.tsv
```

Expected output:

```text
sgRNA	Gene	sample_a	sample_b
guide_a	GENEA	1	0
guide_b	GENEB	0	0
guide_c	GENEC	0	1
```

The output is MAGeCK-compatible: one row per guide, one count column per sample.

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
- `candidates_verified`: native target candidates checked after indexing.

In `sample_a`, the fixture deliberately includes one exact guide assignment, one
ambiguous one-edit read, one unmatched read, and one invalid short read. That is
the behavior DotMatch is designed to expose rather than hide.

## 6. Verify against the checked fixture outputs

```bash
diff -u examples/workflows/fixtures/expected_counts.mageck.tsv \
  tmp/crispr-first-run/counts.mageck.tsv
```

No diff means the tutorial count matrix matches the repository fixture.

For a public-data CRISPR example, use `examples/crispr_guides/run.sh` and the
checked evidence reports under `docs/benchmarks/public_crispr/`.
