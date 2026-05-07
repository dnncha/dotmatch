#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
PREFIX="${DOTMATCH_BARCODE_COMPETITOR_PREFIX:-$ROOT/build/barcode-competitors}"

python3 -m venv "$PREFIX"
"$PREFIX/bin/python" -m pip install --upgrade pip >/dev/null
"$PREFIX/bin/python" -m pip install \
  "cutadapt==5.2" \
  "ultraplex==0.8.0"

cat <<EOF
DOTMATCH_BARCODE_COMPETITOR_BIN=$PREFIX/bin

Run barcode demux benchmark rows with:

PATH="$PREFIX/bin:\$PATH" python3 scripts/bench_barcode_demux.py --run-cutadapt

For real SRP009896-style data:

python3 scripts/fetch_srp009896_barcode_demo.py --accession SRR391079 --subsample 100000 --barcodes-file /path/to/barcodes.tsv
PATH="$PREFIX/bin:\$PATH" python3 scripts/bench_barcode_demux.py \\
  --reads examples/barcode_demux/data/SRR391079.subsample100000.fastq.gz \\
  --barcodes examples/barcode_demux/data/barcodes.tsv \\
  --barcode-start 0 \\
  --barcode-length 8 \\
  --k 0 \\
  --workflow-name srp009896_srr391079_real_subsample \\
  --run-cutadapt
EOF
