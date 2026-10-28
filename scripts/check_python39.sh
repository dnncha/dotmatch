#!/usr/bin/env bash
set -euo pipefail

PYTHON="${PYTHON39:-python3}"
VERSION="$("$PYTHON" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"

if [[ "$VERSION" != "3.9" ]]; then
  echo "compat39: skipped because $PYTHON is Python $VERSION"
  exit 0
fi

"$PYTHON" -c 'import dotmatch, quickdna; assert dotmatch.distance("ACGT", "AGGT") == 1; assert quickdna.distance_leq("ACGT", "AGGT", 1)'
"$PYTHON" -m pytest -q python/tests/test_examples_typed.py python/tests/test_doctests.py
