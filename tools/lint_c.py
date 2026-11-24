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


def call_args(text: str, open_paren: int) -> str:
    depth = 0
    for index in range(open_paren, len(text)):
        if text[index] == "(":
            depth += 1
        elif text[index] == ")":
            depth -= 1
            if depth == 0:
                return text[open_paren + 1 : index]
    return text[open_paren + 1 :]


def top_level_has(expr: str, needle: str) -> bool:
    depth = 0
    for ch in expr:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == needle and depth == 0:
            return True
    return False


def lint(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    rel = path.relative_to(ROOT)
    problems: list[str] = []
    for match in re.finditer(r"(?<![A-Za-z0-9_])(malloc|realloc)\s*\(", text):
        enclosing = text.rfind("\nstatic ", 0, match.start())
        if re.match(r"\nstatic void \*(malloc|realloc)_array\(", text[enclosing:]):
            continue
        args = call_args(text, match.end() - 1)
        size_arg = args.split(",")[-1] if match.group(1) == "realloc" else args
        if top_level_has(size_arg, "*") or re.search(r"\)\s*\*|\*\s*\(", size_arg):
            problems.append(
                f"{rel}:{line_of(text, match.start())}: {match.group(1)} with a multiplied size can overflow; use "
                f"{match.group(1)}_array(count, size)"
            )
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
    problems.extend(duplicate_fastq_validation(text, rel))
    return problems


FASTQ_ERROR_REASON = re.compile(r'fastq_error\(\s*reader\s*,\s*\d+\s*,\s*"([^"]+)"')


def duplicate_fastq_validation(text: str, rel: Path) -> list[str]:
    first_seen: dict[str, int] = {}
    problems: list[str] = []
    for match in FASTQ_ERROR_REASON.finditer(text):
        reason = match.group(1)
        line = line_of(text, match.start())
        if reason in first_seen:
            problems.append(
                f"{rel}:{line}: FASTQ check \"{reason}\" is already written at line {first_seen[reason]}; a second "
                "copy of a record check lets two readers disagree, so call the shared fastq_validate_* helper"
            )
        else:
            first_seen[reason] = line
    return problems


def lint_makefile(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    if re.search(r"^\.SHELLFLAGS\s*:=\s*-o pipefail", text, flags=re.MULTILINE):
        return []
    return [f"{path.relative_to(ROOT)}:1: recipes must run with .SHELLFLAGS := -o pipefail so a piped check cannot hide its exit code"]


def main() -> int:
    paths = [Path(arg).resolve() for arg in sys.argv[1:]] or SOURCES
    problems = [problem for path in paths for problem in lint(path)]
    if not sys.argv[1:]:
        problems.extend(lint_makefile(ROOT / "Makefile"))
    for problem in problems:
        print(problem, file=sys.stderr)
    if problems:
        print(f"lint_c: {len(problems)} problem(s)", file=sys.stderr)
        return 1
    print("lint_c: ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
