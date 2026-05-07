# Barcode Demultiplexing Example

This example is the real-data path for proving barcode-demultiplexing performance.

The first target dataset is SRP009896, a public maize GBS dataset described in public Cutadapt demultiplexing examples as 5-prime inline barcode reads with 96 demultiplexed FASTQ outputs. The example page links `ExampleDataset.zip` on Google Drive with FASTQ files and the barcode file; the fetch script uses ENA for FASTQ access and expects the barcode sheet via `--barcodes-file` or `--barcodes-url`.

Fetch metadata or a real FASTQ subsample:

```bash
python3 scripts/fetch_srp009896_barcode_demo.py --metadata-only

python3 scripts/fetch_srp009896_barcode_demo.py \
  --accession SRR391079 \
  --subsample 100000 \
  --barcodes-file /path/to/barcodes.tsv
```

Run DotMatch versus Cutadapt once the matching barcode file is installed:

```bash
make barcode-competitor-env

PATH="$PWD/build/barcode-competitors/bin:$PATH" python3 scripts/bench_barcode_demux.py \
  --reads examples/barcode_demux/data/SRR391079.subsample100000.fastq.gz \
  --barcodes examples/barcode_demux/data/barcodes.tsv \
  --barcode-start 0 \
  --barcode-length 8 \
  --k 0 \
  --workflow-name srp009896_srr391079_real_subsample \
  --run-cutadapt

python3 scripts/generate_barcode_demux_report.py
```

The built-in benchmark fixture is only a smoke test. State-of-the-art barcode claims require real public FASTQ data, the matching sample/barcode sheet, repeated runs, and fair comparator rows.
