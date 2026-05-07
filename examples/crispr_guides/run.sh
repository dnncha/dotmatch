#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
HERE=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
DATA="${DOTMATCH_EXAMPLE_DATA_DIR:-$HERE/data}"
OUT="${DOTMATCH_EXAMPLE_OUT_DIR:-$HERE/output}"
mkdir -p "$OUT"

if [ ! -f "$DATA/yusa_library.csv" ]; then
  python3 "$ROOT/scripts/fetch_mageck_demo.py" --small --out "$DATA"
fi

READ1="$DATA/ERR376998.fastq.gz"
READ2="$DATA/ERR376999.fastq.gz"
TARGETS="$DATA/yusa_library.csv"

"$ROOT/dotmatch" count \
  --targets "$TARGETS" \
  --reads "$READ1" \
  --reads "$READ2" \
  --sample-label plasmid,ESC1 \
  --target-start 23 \
  --target-length 19 \
  --k 1 \
  --metric levenshtein \
  --indel-window 1 \
  --out "$OUT/counts.tsv" \
  --assignments "$OUT/assignments.tsv" \
  --summary "$OUT/summary.json" \
  --ambiguous report \
  --ambiguous-out "$OUT/ambiguous.tsv" \
  --unmatched-out "$OUT/unmatched.tsv"

"$ROOT/dotmatch" count \
  --targets "$TARGETS" \
  --reads "$READ1" \
  --reads "$READ2" \
  --sample-label plasmid,ESC1 \
  --target-start 23 \
  --target-length 19 \
  --k 1 \
  --metric levenshtein \
  --indel-window 1 \
  --format mageck \
  --out "$OUT/counts.mageck.tsv"

if command -v mageck >/dev/null 2>&1 && [ "${DOTMATCH_EXAMPLE_FULL:-0}" = "1" ]; then
  (cd "$OUT" && mageck count -l "$TARGETS" -n mageck_exact --sample-label "plasmid,ESC1" --trim-5 23 --fastq "$READ1" "$READ2")
else
  printf '%s\n' "MAGeCK comparator skipped; install mageck and set DOTMATCH_EXAMPLE_FULL=1 to run it." > "$OUT/mageck_skipped.txt"
fi

echo "wrote $OUT"
