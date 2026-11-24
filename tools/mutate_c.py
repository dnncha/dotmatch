#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TARGET = ROOT / "src" / "qdalign.c"
SWAPS = {"<": "<=", "<=": "<", ">": ">=", ">=": ">", "==": "!=", "!=": "=="}
OPERATOR = re.compile(r"(?<![<>=!\-])(<=|>=|==|!=|<|>)(?![<>=])")
CC = os.environ.get("CC", "cc")
EXTRA_CFLAGS = os.environ.get("MUTATE_CFLAGS", "").split()
COMMENT = re.compile(r"/\*.*?\*/|//[^\n]*", re.S)


@dataclass(frozen=True)
class Mutant:
    offset: int
    line: int
    original: str
    replacement: str
    source_line: str


def without_comments(text: str) -> str:
    return COMMENT.sub(lambda m: re.sub(r"[^\n]", " ", m.group(0)), text)


def mutants(text: str) -> list[Mutant]:
    found: list[Mutant] = []
    code_text = without_comments(text)
    line_starts = [0] + [m.end() for m in re.finditer("\n", text)]
    for line_index, start in enumerate(line_starts):
        end = text.find("\n", start)
        stop = end if end >= 0 else len(text)
        line = text[start:stop]
        code = code_text[start:stop]
        if code.lstrip().startswith("#"):
            continue
        code = re.sub(r'"(?:[^"\\]|\\.)*"|\'(?:[^\'\\]|\\.)*\'', lambda m: " " * len(m.group(0)), code)
        for match in OPERATOR.finditer(code):
            op = match.group(1)
            found.append(Mutant(start + match.start(), line_index + 1, op, SWAPS[op], line.strip()))
    return found


def run(argv: list[str], cwd: Path, timeout: float) -> bool:
    try:
        return subprocess.run(argv, cwd=cwd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=timeout).returncode == 0
    except subprocess.TimeoutExpired:
        return False


def compile_object(source: Path, output: Path) -> bool:
    return run([CC, "-O1", "-std=c11", *EXTRA_CFLAGS, f"-I{ROOT / 'include'}", "-c", str(source), "-o", str(output)],
               source.parent, 120)


def baseline_object(text: str, work: Path) -> bytes:
    (work / "baseline.c").write_text(text, encoding="utf-8")
    if not compile_object(work / "baseline.c", work / "baseline.o"):
        raise SystemExit("mutate_c: unmutated src/qdalign.c does not compile")
    return (work / "baseline.o").read_bytes()


def outcome(mutant: Mutant, text: str, baseline: bytes, fuzz_iters: int) -> str:
    mutated = text[: mutant.offset] + mutant.replacement + text[mutant.offset + len(mutant.original) :]
    with tempfile.TemporaryDirectory(prefix="dotmatch-mutant-") as raw:
        work = Path(raw)
        (work / "qdalign.c").write_text(mutated, encoding="utf-8")
        if not compile_object(work / "qdalign.c", work / "qdalign.o"):
            return "killed"
        if (work / "qdalign.o").read_bytes() == baseline:
            return "not built"
        include = f"-I{ROOT / 'include'}"
        for name, harness in (("unit", "test_qdalign.c"), ("fuzz", "fuzz_qdalign_diff.c")):
            argv = [CC, "-O1", "-std=c11", *EXTRA_CFLAGS, include, str(ROOT / "tests" / harness), "qdalign.o", "-o", name]
            if not run(argv, work, 120):
                return "killed"
        if not run(["./unit"], work, 60):
            return "killed"
        return "survived" if run(["./fuzz", str(fuzz_iters), "1"], work, 60) else "killed"


EQUIVALENT = ROOT / "tools" / "mutate_c_equivalent.tsv"


def equivalent_mutants() -> set[tuple[str, str, str]]:
    if not EQUIVALENT.exists():
        return set()
    rows = set()
    for line in EQUIVALENT.read_text(encoding="utf-8").splitlines()[1:]:
        fields = line.split("\t")
        if len(fields) == 4 and fields[3].strip():
            rows.add((fields[0], fields[1], fields[2]))
    return rows


def base_exists(base: str) -> bool:
    return subprocess.run(["git", "rev-parse", "--verify", "--quiet", base + "^{commit}"], cwd=ROOT,
                          capture_output=True).returncode == 0


def changed_lines(base: str) -> set[int]:
    diff = subprocess.run(
        ["git", "diff", "-U0", base, "--", str(TARGET.relative_to(ROOT))], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout
    lines: set[int] = set()
    for match in re.finditer(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@", diff, re.M):
        first, count = int(match.group(1)), int(match.group(2) or "1")
        lines.update(range(first, first + count))
    return lines


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", type=int, default=0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--changed-since")
    parser.add_argument("--fuzz-iters", type=int, default=3000)
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = parser.parse_args()
    if shutil.which(CC) is None:
        print(f"mutate_c: compiler {CC} not found", file=sys.stderr)
        return 2
    text = TARGET.read_text(encoding="utf-8")
    pool = mutants(text)
    if args.changed_since:
        if not base_exists(args.changed_since):
            print(f"mutate_c: SKIPPED, base {args.changed_since} is not in this checkout; changed-line mutants not run")
            return 0
        lines = changed_lines(args.changed_since)
        pool = [m for m in pool if m.line in lines]
    if args.sample and len(pool) > args.sample:
        pool = sorted(random.Random(args.seed).sample(pool, args.sample), key=lambda m: m.offset)
    with tempfile.TemporaryDirectory(prefix="dotmatch-mutant-base-") as raw:
        baseline = baseline_object(text, Path(raw))
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        outcomes = list(executor.map(lambda m: (m, outcome(m, text, baseline, args.fuzz_iters)), pool))
    results = [(m, state != "survived") for m, state in outcomes if state != "not built"]
    not_built = len(outcomes) - len(results)
    known = equivalent_mutants()
    survivors = [m for m, dead in results if not dead]
    equivalent = [m for m in survivors if (m.source_line, m.original, m.replacement) in known]
    missed = [m for m in survivors if m not in equivalent]
    rel = TARGET.relative_to(ROOT)
    for m in missed:
        print(f"{rel}:{m.line}: mutant survived: '{m.original}' -> '{m.replacement}' in: {m.source_line}")
    print(
        f"mutate_c: {len(results) - len(survivors)}/{len(results)} mutants killed by test_qdalign + fuzz_qdalign_diff, "
        f"{len(equivalent)} listed as equivalent in {EQUIVALENT.relative_to(ROOT)}, {len(missed)} missing assertions, "
        f"{not_built} not compiled into this build"
    )
    return 1 if missed and not args.report_only else 0


if __name__ == "__main__":
    raise SystemExit(main())
