#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCES = [ROOT / "src" / "qda.c", ROOT / "src" / "qdalign.c", ROOT / "src" / "target_table.h"]

BANNED_CALLS = {
    "fgets": "fgets splits long lines and hides NUL bytes; use read_text_line or the FASTQ reader",
    "gets": "gets has no bound",
    "sprintf": "use snprintf",
    "vsprintf": "use vsnprintf",
    "strcpy": "use memcpy with a checked length",
    "strcat": "use memcpy with a checked length",
    "strtok": "strtok keeps hidden global state; use split_fields",
    "atoi": "atoi cannot report errors; use parse_int_value",
    "atol": "atol cannot report errors; use parse_size_value",
    "atof": "atof cannot report errors; use parse_double_value",
}

OUTPUT_ASSIGN = re.compile(r"(?:FILE\s*\*\s*)?([A-Za-z_][A-Za-z0-9_\[\]\.>-]*)\s*=\s*open_output_file\(")


def functions(text: str) -> list[tuple[int, str]]:
    found: list[tuple[int, str]] = []
    depth = 0
    start = -1
    for i, ch in enumerate(text):
        if ch == "{":
            if depth == 0:
                start = i
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0 and start >= 0:
                found.append((start, text[start : i + 1]))
                start = -1
    return found


def line_of(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def lint(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    rel = path.relative_to(ROOT)
    problems: list[str] = []
    for name, why in BANNED_CALLS.items():
        for match in re.finditer(rf"(?<![A-Za-z0-9_]){name}\s*\(", text):
            problems.append(f"{rel}:{line_of(text, match.start())}: banned call {name}(): {why}")
    for start, body in functions(text):
        names = {m.group(1) for m in OUTPUT_ASSIGN.finditer(body)}
        for name in sorted(names):
            escaped = re.escape(name)
            if re.search(rf"return\s+{escaped}\s*;", body):
                continue
            for match in re.finditer(rf"(?<![A-Za-z0-9_])fclose\(\s*{escaped}\s*\)", body):
                problems.append(
                    f"{rel}:{line_of(text, start + match.start())}: fclose({name}) on an output file drops write "
                    "errors; close it with close_output_file"
                )
            if not re.search(rf"close_output_file\(\s*{escaped}\s*[,)]", body):
                problems.append(
                    f"{rel}:{line_of(text, start)}: output file {name} is never closed with close_output_file, so "
                    "a failed write can still exit 0"
                )
    return problems


def main() -> int:
    paths = [Path(arg).resolve() for arg in sys.argv[1:]] or SOURCES
    problems = [problem for path in paths for problem in lint(path)]
    for problem in problems:
        print(problem, file=sys.stderr)
    if problems:
        print(f"lint_c: {len(problems)} problem(s)", file=sys.stderr)
        return 1
    print("lint_c: ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
