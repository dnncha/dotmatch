"""Dependency-free source hygiene: syntax, UTF-8, indentation, and trailing spaces.

This intentionally does not claim to replace a full static type checker or Ruff.
"""
from __future__ import annotations

import ast
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    errors: list[str] = []
    paths = sorted(path for directory in ("src", "tests", "scripts", "benchmarks")
                   for path in (root / directory).rglob("*.py"))
    for path in paths:
        text = path.read_text(encoding="utf-8")
        label = str(path.relative_to(root))
        try:
            ast.parse(text, filename=label, feature_version=(3, 11))
        except SyntaxError as error:
            errors.append(f"{label}: {error}")
        for number, line in enumerate(text.splitlines(), 1):
            if line.rstrip() != line:
                errors.append(f"{label}:{number}: trailing whitespace")
            if "\t" in line[:len(line) - len(line.lstrip())]:
                errors.append(f"{label}:{number}: tab indentation")
        if text and not text.endswith("\n"):
            errors.append(f"{label}: missing final newline")
    if errors:
        print("\n".join(errors))
        return 1
    print(f"Source hygiene passed for {len(paths)} Python files (not a type check).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
