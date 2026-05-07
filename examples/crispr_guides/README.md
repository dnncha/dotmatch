# CRISPR Guide Counting Example

This example mirrors the MAGeCK/Yusa demo workflow: two FASTQ samples, a guide library, 23 bases of 5-prime sequence before the guide, and 19 nt guide targets.

Create a small local fixture:

```bash
python3 ../../scripts/fetch_mageck_demo.py --small --out data
./run.sh
```

Fetch the full public MAGeCK demo data instead:

```bash
python3 ../../scripts/fetch_mageck_demo.py --out data
DOTMATCH_EXAMPLE_FULL=1 ./run.sh
```

Outputs are written under `output/`:

- `counts.tsv`: detailed DotMatch count table;
- `counts.mageck.tsv`: MAGeCK-compatible `sgRNA Gene sample...` table;
- `summary.json`: assignment summary;
- `assignments.tsv`: per-read diagnostics for the fixture/full run.

The full public workflow uses `ERR376998.fastq.gz`, `ERR376999.fastq.gz`, `yusa_library.csv`, `--target-start 23`, and `--target-length 19`.
