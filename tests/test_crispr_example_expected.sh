#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
TMP_BASE="${TMPDIR:-/tmp}"
TMP=$(mktemp -d "$TMP_BASE/dotmatch-crispr-example.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

DOTMATCH_EXAMPLE_DATA_DIR="$TMP/data" \
DOTMATCH_EXAMPLE_OUT_DIR="$TMP/output" \
  "$ROOT/examples/crispr_guides/run.sh" >/dev/null

EXPECTED="$ROOT/examples/crispr_guides/expected_output"
for name in counts.tsv counts.mageck.tsv assignments.tsv ambiguous.tsv unmatched.tsv mageck_skipped.txt; do
  diff -u "$EXPECTED/$name" "$TMP/output/$name"
done

python3 - "$TMP/output/summary.json" "$EXPECTED/summary.stable.json" <<'PY'
import json
import sys
from pathlib import Path

summary_path = Path(sys.argv[1])
expected_path = Path(sys.argv[2])

top_keys = [
    "k",
    "metric",
    "ambiguity_policy",
    "indel_window",
    "target_start",
    "auto_offset",
    "offset_mode",
    "offset_detection_strategy",
    "count_engine",
    "hamming_index",
    "target_length",
    "n_targets",
    "read_threads",
]
sample_keys = [
    "sample",
    "selected_target_start",
    "selected_target_starts",
    "total_reads",
    "assigned_unique",
    "assigned_exact",
    "assigned_corrected",
    "k1_rescued_reads",
    "percent_rescued_by_k1",
    "ambiguous",
    "percent_ambiguous",
    "unmatched",
    "percent_unmatched",
    "invalid",
    "library_covered_targets",
    "library_coverage_fraction",
    "top_target_id",
    "top_target_count",
    "candidates_considered",
    "candidates_verified",
]

def stable_projection(summary):
    projected = {key: summary[key] for key in top_keys}
    projected["samples"] = [
        {key: sample[key] for key in sample_keys}
        for sample in summary["samples"]
    ]
    return projected

actual = stable_projection(json.loads(summary_path.read_text(encoding="utf-8")))
expected = json.loads(expected_path.read_text(encoding="utf-8"))
if actual != expected:
    print("Stable summary fields differ from expected output.", file=sys.stderr)
    print("Actual:", json.dumps(actual, indent=2, sort_keys=True), file=sys.stderr)
    print("Expected:", json.dumps(expected, indent=2, sort_keys=True), file=sys.stderr)
    raise SystemExit(1)
PY
