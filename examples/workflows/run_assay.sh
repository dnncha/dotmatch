#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FIXTURES="${ROOT}/examples/workflows/fixtures"

if [[ ! -x "${DOTMATCH_NATIVE_CLI:-${ROOT}/dotmatch}" ]]; then
  make -C "${ROOT}" dotmatch
fi

export DOTMATCH_NATIVE_CLI="${DOTMATCH_NATIVE_CLI:-${ROOT}/dotmatch}"
export PYTHONPATH="${ROOT}/python${PYTHONPATH:+:${PYTHONPATH}}"

cd "${FIXTURES}"
set +e
if command -v dotmatch >/dev/null 2>&1; then
  dotmatch assay start crispr_assay.toml "$@"
else
  python3 -m dotmatch.cli assay start crispr_assay.toml "$@"
fi
rc=$?
set -e

if [[ "$rc" -eq 0 ]]; then
  echo "done: assay_out/reliability_report.html (passed)" >&2
elif [[ "$rc" -eq 1 ]]; then
  echo "done: assay_out/reliability_report.html (needs review)" >&2
else
  echo "done: assay_out/reliability_report.html (failed or blocked)" >&2
fi
exit "$rc"