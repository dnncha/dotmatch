#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

mkdir -p benchmarks/raw/atlas docs/benchmarks/crispr_sota benchmarks/figures
LOG="benchmarks/raw/atlas/crispr_sota_$(date -u +%Y%m%dT%H%M%SZ).log"
exec > >(tee "$LOG") 2>&1

echo "dotmatch_commit,$(git rev-parse HEAD 2>/dev/null || echo unknown)" > benchmarks/raw/atlas/crispr_sota_environment.csv
if command -v sha256sum >/dev/null 2>&1; then
  find Makefile README.md include src tests scripts python docs -type f -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    | sha256sum \
    | awk '{print "source_tree_sha256," $1}' >> benchmarks/raw/atlas/crispr_sota_environment.csv
fi
echo "uname,$(uname -a)" >> benchmarks/raw/atlas/crispr_sota_environment.csv
if command -v lscpu >/dev/null 2>&1; then
  lscpu | sed 's/,/;/g' > benchmarks/raw/atlas/crispr_sota_lscpu.txt
fi
if command -v cc >/dev/null 2>&1; then
  cc --version | head -n 1 | sed 's/^/cc,/' >> benchmarks/raw/atlas/crispr_sota_environment.csv
fi
python3 --version | sed 's/^/python,/' >> benchmarks/raw/atlas/crispr_sota_environment.csv

make clean
make
make test cli-test python-test coverage
make competitor-env edlib-tools

PATH="$ROOT/build/guide-counter/bin:$ROOT/build/competitor-env/bin:$PATH" \
  DOTMATCH_SOTA_READ_SIZES="${DOTMATCH_SOTA_READ_SIZES:-10000,100000}" \
  DOTMATCH_SOTA_REPEATS="${DOTMATCH_SOTA_REPEATS:-5}" \
  DOTMATCH_COUNT_THREADS="${DOTMATCH_COUNT_THREADS:-1}" \
  python3 scripts/run_crispr_sota_repeated.py --run-mageck --run-guide-counter --full

DOTMATCH_SOTA_VALIDATION_RECORDS="${DOTMATCH_SOTA_VALIDATION_RECORDS:-100000}" \
DOTMATCH_SOTA_VALIDATION_SAMPLE="${DOTMATCH_SOTA_VALIDATION_SAMPLE:-100000}" \
  python3 scripts/validate_crispr_sota_edlib.py

python3 scripts/compare_crispr_sota_counts.py
python3 scripts/generate_crispr_sota_report.py
python3 scripts/check_crispr_sota_gate.py

echo "CRISPR SOTA atlas run complete: $LOG"
