"""Inline snapshot expectations for deterministic example tests.

The interface mirrors Rust's ``expect-test``: a test writes the value it wants
to pin, the first run (with ``UPDATE_EXPECT=1``) records it under
``python/tests/expected/``, and later runs must match byte-for-byte.
"""

from __future__ import annotations

import difflib
import os
from pathlib import Path

EXPECTED_DIR = Path(__file__).with_name("expected")


def expect_snapshot(actual: str, *, name: str) -> None:
    """Assert that ``actual`` matches the named snapshot, or record it when updating."""
    path = EXPECTED_DIR / f"{name}.txt"
    normalized = actual.rstrip("\n") + "\n"
    if os.environ.get("UPDATE_EXPECT") == "1":
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(normalized, encoding="utf-8")
        return
    if not path.exists():
        raise AssertionError(
            f"missing expected snapshot {path}; run UPDATE_EXPECT=1 pytest python/tests/test_examples_typed.py"
        )
    expected = path.read_text(encoding="utf-8")
    if expected != normalized:
        diff = "".join(
            difflib.unified_diff(
                expected.splitlines(keepends=True),
                normalized.splitlines(keepends=True),
                fromfile=f"expected/{name}.txt",
                tofile="actual",
            )
        )
        raise AssertionError(f"snapshot mismatch for {name}:\n{diff}")
