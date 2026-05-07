#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

mkdir -p benchmarks/raw/atlas docs/benchmarks/barcode_demux benchmarks/figures
LOG="benchmarks/raw/atlas/barcode_sota_$(date -u +%Y%m%dT%H%M%SZ).log"
exec > >(tee "$LOG") 2>&1

echo "dotmatch_commit,$(git rev-parse HEAD 2>/dev/null || echo unknown)" > benchmarks/raw/atlas/barcode_sota_environment.csv
if command -v sha256sum >/dev/null 2>&1; then
  find Makefile README.md include src tests scripts python docs -type f -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    | sha256sum \
    | awk '{print "source_tree_sha256," $1}' >> benchmarks/raw/atlas/barcode_sota_environment.csv
fi
echo "uname,$(uname -a)" >> benchmarks/raw/atlas/barcode_sota_environment.csv
if command -v lscpu >/dev/null 2>&1; then
  lscpu | sed 's/,/;/g' > benchmarks/raw/atlas/barcode_sota_lscpu.txt
fi
python3 --version | sed 's/^/python,/' >> benchmarks/raw/atlas/barcode_sota_environment.csv

if [[ -z "${DOTMATCH_BARCODE_SOTA_BARCODES:-}" && -z "${DOTMATCH_BARCODE_SOTA_BARCODES_URL:-}" ]]; then
  echo "Set DOTMATCH_BARCODE_SOTA_BARCODES or DOTMATCH_BARCODE_SOTA_BARCODES_URL to a real barcode sheet." >&2
  exit 2
fi

FETCH_ARGS=(--out examples/barcode_demux/data --require-barcodes --subsample "${DOTMATCH_BARCODE_SOTA_SUBSAMPLE:-0}" --barcode-start "${DOTMATCH_BARCODE_START:-0}")
if [[ -n "${DOTMATCH_BARCODE_LENGTH:-}" ]]; then
  FETCH_ARGS+=(--barcode-length "$DOTMATCH_BARCODE_LENGTH")
fi
if [[ -n "${DOTMATCH_BARCODE_SOTA_BARCODES:-}" ]]; then
  FETCH_ARGS+=(--barcodes-file "$DOTMATCH_BARCODE_SOTA_BARCODES")
else
  FETCH_ARGS+=(--barcodes-url "$DOTMATCH_BARCODE_SOTA_BARCODES_URL")
fi
python3 scripts/fetch_srp009896_barcode_demo.py "${FETCH_ARGS[@]}"

READS="$(python3 - <<'PY'
import json
from pathlib import Path
meta=json.loads(Path("examples/barcode_demux/data/metadata.json").read_text())
print(meta["runs"][0]["local_fastq"])
PY
)"
BARCODES="$(python3 - <<'PY'
import json
from pathlib import Path
meta=json.loads(Path("examples/barcode_demux/data/metadata.json").read_text())
print(meta["barcodes"])
PY
)"
BARCODE_LENGTH="$(python3 - <<'PY'
import json
from pathlib import Path
meta=json.loads(Path("examples/barcode_demux/data/metadata.json").read_text())
print(meta.get("barcode_length") or 8)
PY
)"

make clean
make
make test cli-test python-test coverage
make barcode-competitor-env

PATH="$ROOT/build/barcode-competitors/bin:$PATH" \
  python3 scripts/bench_barcode_demux.py \
    --reads "$READS" \
    --barcodes "$BARCODES" \
    --barcode-start "${DOTMATCH_BARCODE_START:-0}" \
    --barcode-length "$BARCODE_LENGTH" \
    --k "${DOTMATCH_BARCODE_K:-1}" \
    --metric "${DOTMATCH_BARCODE_METRIC:-hamming}" \
    --workflow-name "real_srp009896_inline_barcode" \
    --run-cutadapt \
    --repeats "${DOTMATCH_BARCODE_REPEATS:-5}"

python3 scripts/generate_barcode_demux_report.py
python3 scripts/check_barcode_sota_gate.py --no-second-comparator

echo "Barcode SOTA atlas run complete: $LOG"
