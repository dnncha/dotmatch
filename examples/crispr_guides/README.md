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
- `counts.mageck.tsv`: MAGeCK-ready `sgRNA`, `Gene`, then one count column per
  sample;
- `summary.json`: selected counting settings and per-sample assignment summary;
- `assignments.tsv`: per-read diagnostics for the fixture/full run.

DotMatch reads common guide-library headers including `sgRNA`, `sgRNAID`,
`guide_id`, `gRNA.sequence`, `sgRNA_sequence`, `guide_seq`, `sequence`, `Gene`,
and `gene_symbol`.

The tiny fixture outputs are tracked under `expected_output/`, and `make cli-test` checks them. Set `DOTMATCH_EXAMPLE_DATA_DIR` and `DOTMATCH_EXAMPLE_OUT_DIR` to run the example against alternate local paths.

The full public workflow uses `ERR376998.fastq.gz`, `ERR376999.fastq.gz`, `yusa_library.csv`, `--target-start 23`, and `--target-length 19`.

## Run it with the published package

The example can use the installed command instead of a locally built binary:

```bash
git clone https://github.com/dnncha/dotmatch.git
cd dotmatch
python3 -m pip install dotmatch
DOTMATCH_BIN=dotmatch ./examples/crispr_guides/run.sh
```

The default command fetches the small public fixture and writes the count,
assignment, and summary artifacts under `examples/crispr_guides/output/`.
To fetch the full public MAGeCK/Yusa demo instead, run:

```bash
python3 scripts/fetch_mageck_demo.py --out examples/crispr_guides/data
DOTMATCH_EXAMPLE_FULL=1 DOTMATCH_BIN=dotmatch ./examples/crispr_guides/run.sh
```

The checked comparison and count-agreement results are summarized in the
[public CRISPR benchmark report](../../docs/benchmarks/public_crispr/README.md).
