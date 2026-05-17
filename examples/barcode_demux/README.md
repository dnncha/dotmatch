# Barcode Demultiplexing Example

This example is the real-data path for proving barcode-demultiplexing performance.

The first target dataset is SRP009896, a public maize GBS dataset described in public Cutadapt demultiplexing examples as 5-prime inline barcode reads with 96 demultiplexed FASTQ outputs. The example page links `ExampleDataset.zip` on Google Drive with FASTQ files and the barcode file; the fetch script uses ENA for FASTQ access and can extract the first ZIP member barcode sheet with a ranged request.

Fetch metadata and the public barcode sheet without downloading the 7.4 GB Google Drive archive:

```bash
python3 scripts/fetch_srp009896_barcode_demo.py \
  --metadata-only \
  --use-public-example-barcodes \
  --require-barcodes
```

Fetch a real FASTQ subsample:

```bash
python3 scripts/fetch_srp009896_barcode_demo.py \
  --accession SRR391079 \
  --subsample 100000 \
  --use-public-example-barcodes \
  --require-barcodes
```

Run DotMatch versus Cutadapt once the matching barcode file is installed:

```bash
make barcode-competitor-env

PATH="$PWD/build/barcode-competitors/bin:$PATH" python3 scripts/bench_barcode_demux.py \
  --reads examples/barcode_demux/data/SRR391079.subsample100000.fastq.gz \
  --barcodes examples/barcode_demux/data/barcodes.tsv \
  --barcode-start 1 \
  --barcode-length auto \
  --k 0 \
  --workflow-name srp009896_srr391079_real_subsample \
  --run-cutadapt \
  --run-hash-splitter

python3 scripts/generate_barcode_demux_report.py
```

For a workflow-facing diagnostic report, run the barcode autopsy demo:

```bash
make barcode-autopsy-demo
```

This writes `examples/barcode_autopsy/results/report.html` plus offset-scan,
collision-audit, assignment, unmatched-read, provenance, and MultiQC custom
content artifacts. The report is diagnostic; comparator throughput
claims still come from the checked comparison gate above.

For the broader fixed-window barcode validation checks, run:

```bash
make barcode-validation-ready
```

That command checks the public fixed-window matrix in
`docs/barcode-science-readiness.json` and the failure-mode fixtures under
`examples/barcode_autopsy/failure_modes/`.

The SRP009896 barcode sheet contains variable-length barcodes (`4-8 bp`) and separate run blocks with reused barcode sequences; the fetcher filters to the requested accession when that run column is present. The SRP009896 reads include a leading `N`, so use `--barcode-start 1` for this public example. DotMatch supports this with `--barcode-length auto` and conservative ambiguity handling for prefix-overlapping barcodes. A full comparison claim still requires repeated real-data rows and comparator evidence that pass `make barcode-comparison-gate`.

The built-in benchmark fixture is only a smoke test. Comparative barcode claims require real public FASTQ data, the matching sample/barcode sheet, repeated runs, and fair comparator rows.
