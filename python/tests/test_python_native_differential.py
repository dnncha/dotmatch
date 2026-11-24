from __future__ import annotations

import json
import os
import random
import subprocess
from collections import Counter
from pathlib import Path

import pytest

from dotmatch import MATCH_UNIQUE, assignment_summary, stream_assign
from dotmatch._types import AssignmentPolicy, Metric

ROOT = Path(__file__).resolve().parents[2]

COMBINATIONS: list[tuple[Metric, int, AssignmentPolicy]] = [
    ("hamming", 0, "best"),
    ("hamming", 0, "radius"),
    ("hamming", 1, "best"),
    ("hamming", 1, "radius"),
    ("hamming", 2, "best"),
    ("hamming", 2, "radius"),
    ("hamming", 3, "best"),
    ("hamming", 3, "radius"),
    ("levenshtein", 0, "best"),
    ("levenshtein", 0, "radius"),
    ("levenshtein", 1, "best"),
    ("levenshtein", 1, "radius"),
    ("levenshtein", 2, "best"),
    ("levenshtein", 2, "radius"),
]


def _native_binary() -> Path | None:
    env = os.environ.get("DOTMATCH_BIN")
    if env:
        candidate = Path(env)
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
    candidate = ROOT / "dotmatch"
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return candidate
    return None


def _mutate_sequence(rng: random.Random, seq: str, edits: int, allow_indel: bool) -> str:
    chars = list(seq)
    for _ in range(edits):
        op = rng.choice(["sub", "sub", "ins", "del"] if allow_indel else ["sub"])
        if op == "sub" and chars:
            pos = rng.randrange(len(chars))
            alphabet = [b for b in "ACGTN" if b != chars[pos].upper()]
            base = rng.choice(alphabet)
            chars[pos] = base.lower() if rng.random() < 0.15 else base
        elif op == "ins":
            pos = rng.randrange(len(chars) + 1)
            chars.insert(pos, rng.choice("ACGT"))
        elif op == "del" and len(chars) > 1:
            pos = rng.randrange(len(chars))
            chars.pop(pos)
    return "".join(chars)


def _generate_case(
    tmp_path: Path, seed: int
) -> tuple[Path, Path, int, int, list[tuple[str, str]]]:
    rng = random.Random(seed)
    target_len = rng.choice([8, 12, 16, 20])
    target_start = rng.choice([0, 3, 5])
    n_targets = rng.randint(6, 14)

    base_targets: list[str] = []
    for idx in range(n_targets):
        if base_targets and rng.random() < 0.45:
            parent = rng.choice(base_targets)
            seq = _mutate_sequence(rng, parent, rng.choice([1, 2]), allow_indel=False)
        else:
            seq = "".join(rng.choices("ACGT", k=target_len))
        base_targets.append(seq.upper())

    targets = [(f"T{idx}", seq) for idx, seq in enumerate(base_targets)]
    targets_path = tmp_path / f"targets_{seed}.tsv"
    targets_path.write_text(
        "id\tsequence\n" + "".join(f"{tid}\t{seq}\n" for tid, seq in targets),
        encoding="utf-8",
    )

    n_reads = 40
    fq_lines: list[str] = []
    for r_idx in range(n_reads):
        kind = rng.randrange(10)
        if kind == 0:
            # Too short -> invalid window
            read_seq = "".join(rng.choices("ACGT", k=max(1, target_start + target_len - 2)))
        else:
            prefix = "".join(rng.choices("ACGT", k=target_start))
            suffix = "".join(rng.choices("ACGT", k=rng.randint(0, 6)))
            if kind <= 3:
                window = rng.choice(base_targets)
                if rng.random() < 0.2:
                    window = window.lower()
            elif kind <= 7:
                window = _mutate_sequence(
                    rng,
                    rng.choice(base_targets),
                    rng.choice([1, 2, 3]),
                    allow_indel=(kind == 7),
                )
            else:
                window = "".join(rng.choices("ACGTN", k=target_len))
            read_seq = prefix + window + suffix
        qual = "I" * len(read_seq)
        fq_lines.append(f"@read_{seed}_{r_idx}\n{read_seq}\n+\n{qual}\n")

    fastq_path = tmp_path / f"reads_{seed}.fastq"
    fastq_path.write_text("".join(fq_lines), encoding="utf-8")
    return targets_path, fastq_path, target_start, target_len, targets


def _parse_native_counts(path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    id_col = 0
    count_col = header.index("count_total") if "count_total" in header else len(header) - 1
    for line in lines[1:]:
        parts = line.split("\t")
        counts[parts[id_col]] = int(parts[count_col])
    return counts


@pytest.mark.parametrize("seed", range(1, 31))
def test_python_stream_assign_matches_native_count(tmp_path: Path, seed: int) -> None:
    binary = _native_binary()
    if binary is None:
        pytest.skip("native dotmatch binary not built")

    targets_path, fastq_path, target_start, target_len, targets = _generate_case(tmp_path, seed)
    metric, k, policy = COMBINATIONS[(seed - 1) % len(COMBINATIONS)]

    assignments = list(
        stream_assign(
            fastq_path,
            targets_path,
            target_start=target_start,
            target_length=target_len,
            k=k,
            metric=metric,
            policy=policy,
        )
    )
    py_summary = assignment_summary(assignments)
    py_counts: Counter[str] = Counter()
    for row in assignments:
        if row.status == MATCH_UNIQUE:
            py_counts[row.target_name] += 1
    expected_counts = {tid: py_counts.get(tid, 0) for tid, _ in targets}

    out_counts = tmp_path / f"counts_{seed}.tsv"
    out_summary = tmp_path / f"summary_{seed}.json"
    cmd = [
        str(binary),
        "count",
        "--targets",
        str(targets_path),
        "--reads",
        str(fastq_path),
        "--target-start",
        str(target_start),
        "--target-length",
        str(target_len),
        "--k",
        str(k),
        "--metric",
        metric,
        "--ambiguity-policy",
        policy,
        "--threads",
        "1",
        "--out",
        str(out_counts),
        "--summary",
        str(out_summary),
    ]
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    assert proc.returncode == 0, f"native count failed (seed={seed}, {metric}, k={k}, {policy}): {proc.stderr}"

    native_summary = json.loads(out_summary.read_text(encoding="utf-8"))
    native_sample = native_summary["samples"][0]
    native_counts = _parse_native_counts(out_counts)

    assert native_counts == expected_counts, (
        f"per-target count mismatch (seed={seed}, metric={metric}, k={k}, policy={policy})"
    )
    for field in (
        "total_reads",
        "assigned_unique",
        "assigned_exact",
        "assigned_corrected",
        "ambiguous",
        "unmatched",
        "invalid",
    ):
        assert native_sample[field] == py_summary[field], (
            f"summary.{field} mismatch (seed={seed}, metric={metric}, k={k}, policy={policy}): "
            f"native={native_sample[field]} vs python={py_summary[field]}"
        )
