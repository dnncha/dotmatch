# Sanson/Brunello Guide-Counter-Style Example

This example mirrors the public `guide-counter` Sanson/Brunello paper-data lane: four single-end FASTQ samples from PRJNA508200/SRP172473, the Broad Brunello corrected guide library, automatic guide-offset detection, exact plus one-substitution guide matching, and no indels.

Create the local public-data subset and run DotMatch:

```bash
DOTMATCH_SANSON_SUBSAMPLE=100000 DOTMATCH_COUNT_THREADS=4 ./run.sh
```

Run the full FASTQs instead:

```bash
DOTMATCH_EXAMPLE_FULL=1 DOTMATCH_COUNT_THREADS=4 ./run.sh
```

Record the backend optimizer recommendation for the full workload:

```bash
PYTHONPATH=../../python python3 -m dotmatch.cli assay optimize assay.full.toml
```

Use the GuideCounter-compatible command shape on the fetched subset:

```bash
../../dotmatch guide-counter count \
  --input data/plasmid.subsample100000.fastq.gz data/RepA.subsample100000.fastq.gz \
  --samples plasmid RepA \
  --library data/broadgpp-brunello-library-corrected.txt \
  --output output/guide_counts
```

Outputs are written under `output/`:

- `counts.hamming.mageck.tsv`: MAGeCK-ready `sgRNA`, `Gene`, then one count
  column per sample;
- `summary.hamming.json`: selected offsets, assignment counts, QC rates, and
  timing;
- `sample_qc.tsv`: per-sample totals, assignment/no-match/ambiguity rates,
  guide coverage, zero-count guide counts, Gini index, and top-guide fraction;
- `assay_full/backend_optimization.json`: CPU-authoritative optimizer plan for
  the full AssaySpec workload;
- `guide_counts.counts.txt`, `guide_counts.extended-counts.txt`, and
  `guide_counts.stats.txt`: GuideCounter-compatible wrapper outputs when the
  compatibility command is used;
- `guide_counter.counts.txt` and `guide_counter.stats.txt`: optional comparator outputs when `DOTMATCH_RUN_GUIDE_COUNTER=1` and `guide-counter` is on `PATH`.

The DotMatch lane uses `--k 1 --metric hamming --ambiguity-policy best --offset-mode multi --offset-min-fraction 0.0025`, matching guide-counter's one-mismatch/no-indel counting mode and current offset threshold default. The compatibility wrapper applies the same defaults, samples up to 100000 reads for offset detection by default, infers sample names when `--samples` is omitted, and switches to exact-only counting with `--exact-match`. Assignment remains DotMatch's deterministic CPU path; GPU/backend optimizer evidence is advisory and does not change the counts. The broader benchmark and count-agreement gate live in `docs/benchmarks/crispr_comparison/README.md`.
