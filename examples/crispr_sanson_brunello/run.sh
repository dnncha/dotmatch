#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
HERE=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
DATA="${DOTMATCH_EXAMPLE_DATA_DIR:-$HERE/data}"
OUT="${DOTMATCH_EXAMPLE_OUT_DIR:-$HERE/output}"
THREADS="${DOTMATCH_COUNT_THREADS:-4}"

mkdir -p "$DATA" "$OUT"

if [ "${DOTMATCH_EXAMPLE_FULL:-0}" = "1" ]; then
  SUBSAMPLE=0
  SUFFIX=".fastq.gz"
else
  SUBSAMPLE="${DOTMATCH_SANSON_SUBSAMPLE:-100000}"
  SUFFIX=".subsample${SUBSAMPLE}.fastq.gz"
fi

python3 "$ROOT/scripts/fetch_sanson_brunello_demo.py" --out "$DATA" --subsample "$SUBSAMPLE"

if [ "$SUBSAMPLE" = "0" ]; then
  (cd "$HERE" && PYTHONPATH="$ROOT/python${PYTHONPATH:+:$PYTHONPATH}" python3 -m dotmatch.cli assay optimize assay.full.toml)
fi

LIBRARY="$DATA/broadgpp-brunello-library-corrected.txt"
PLASMID="$DATA/plasmid$SUFFIX"
REPA="$DATA/RepA$SUFFIX"
REPB="$DATA/RepB$SUFFIX"
REPC="$DATA/RepC$SUFFIX"

"$ROOT/dotmatch" count \
  --targets "$LIBRARY" \
  --reads "$PLASMID" \
  --reads "$REPA" \
  --reads "$REPB" \
  --reads "$REPC" \
  --sample-label plasmid,RepA,RepB,RepC \
  --target-start 20 \
  --target-length 20 \
  --k 1 \
  --metric hamming \
  --ambiguity-policy best \
  --auto-offset 20 \
  --auto-offset-sample 100000 \
  --offset-mode multi \
  --offset-min-fraction 0.0025 \
  --format mageck \
  --threads "$THREADS" \
  --out "$OUT/counts.hamming.mageck.tsv" \
  --summary "$OUT/summary.hamming.json"

if [ "${DOTMATCH_RUN_GUIDE_COUNTER:-0}" = "1" ] && command -v guide-counter >/dev/null 2>&1; then
  guide-counter count \
    --input "$PLASMID" "$REPA" "$REPB" "$REPC" \
    --samples plasmid RepA RepB RepC \
    --library "$LIBRARY" \
    --output "$OUT/guide_counter"
else
  printf '%s\n' "guide-counter comparator skipped; set DOTMATCH_RUN_GUIDE_COUNTER=1 with guide-counter on PATH to run it." > "$OUT/guide_counter_skipped.txt"
fi

echo "wrote $OUT"
