#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
DOTMATCH_BIN=${DOTMATCH_BIN:-"$ROOT/dotmatch"}

cd "$ROOT"
if [ ! -x "$DOTMATCH_BIN" ]; then
  make dotmatch
fi
python3 scripts/run_perturb_seq_gse146194.py public --dotmatch "$DOTMATCH_BIN"
python3 scripts/check_perturb_seq_gse146194.py --require-public --require-work
