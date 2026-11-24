"""Unified deterministic run triage for known-target short-DNA FASTQ assays.

Combines target-library collision checks, three-policy sensitivity analysis
(exact, radius k=1, best k=1), and unmatched-read diagnostics (offset shifts
and reverse-complemented windows) into one local, structured verdict.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal, Sequence

from .core import MATCH_UNIQUE, Matcher, iter_fastq
from .sensitivity import run_sensitivity
from .target_io import read_target_table

TriageSeverity = Literal["info", "warning", "critical"]
TriageVerdict = Literal["healthy", "review_required", "likely_misconfigured"]

SCHEMA_VERSION = "dotmatch.triage.v1"
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


@dataclass(frozen=True)
class TriageFinding:
    code: str
    severity: TriageSeverity
    summary: str
    evidence: dict[str, Any]
    recommended_action: str


@dataclass(frozen=True)
class TriageReport:
    schema_version: str
    verdict: TriageVerdict
    read_count: int
    target_count: int
    parameters: dict[str, Any]
    exact_assignment_rate: float
    radius_k1_assignment_rate: float
    best_k1_assignment_rate: float
    changed_reads: int
    findings: tuple[TriageFinding, ...]
    artifacts: dict[str, str]

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "verdict": self.verdict,
            "read_count": self.read_count,
            "target_count": self.target_count,
            "parameters": dict(self.parameters),
            "exact_assignment_rate": self.exact_assignment_rate,
            "radius_k1_assignment_rate": self.radius_k1_assignment_rate,
            "best_k1_assignment_rate": self.best_k1_assignment_rate,
            "changed_reads": self.changed_reads,
            "findings": [asdict(finding) for finding in self.findings],
            "artifacts": dict(self.artifacts),
        }


def _derive_verdict(findings: tuple[TriageFinding, ...]) -> TriageVerdict:
    if any(finding.severity == "critical" for finding in findings):
        return "likely_misconfigured"
    if any(finding.severity == "warning" for finding in findings):
        return "review_required"
    return "healthy"


def _rate(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator > 0 else 0.0


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _audit_library_safety(matcher: Matcher, sequences: Sequence[str]) -> list[TriageFinding]:
    findings: list[TriageFinding] = []
    counts = Counter(seq.upper() for seq in sequences)
    dup_seqs = [seq for seq, count in counts.items() if count > 1]
    if dup_seqs:
        dup_target_count = sum(counts[seq] for seq in dup_seqs)
        findings.append(
            TriageFinding(
                code="duplicate_targets",
                severity="critical",
                summary=(
                    f"{dup_target_count} targets share {len(dup_seqs)} duplicate sequence(s); "
                    "reads matching them can never be uniquely assigned."
                ),
                evidence={
                    "duplicate_target_count": dup_target_count,
                    "distinct_duplicate_sequences": len(dup_seqs),
                    "target_count": len(sequences),
                },
                recommended_action="Deduplicate or collapse identical target sequences before counting.",
            )
        )

    if sequences:
        k1_results = matcher.assign_hamming(list(sequences), k=1, policy="radius")
        k1_neighbors = sum(
            1
            for seq, res in zip(sequences, k1_results)
            if counts[seq.upper()] == 1 and res.match_count > 1
        )
        if k1_neighbors > 0:
            findings.append(
                TriageFinding(
                    code="library_collision_k1",
                    severity="warning",
                    summary=(
                        f"{k1_neighbors} of {len(sequences)} unique targets have at least one Hamming k=1 neighbor."
                    ),
                    evidence={
                        "targets_with_k1_neighbor": k1_neighbors,
                        "target_count": len(sequences),
                    },
                    recommended_action=(
                        "Use radius ambiguity policy or exact matching; inspect sensitivity deltas before trusting k>=1 counts."
                    ),
                )
            )
    return findings


def _diagnose_read_orientation_and_offset(
    *,
    matcher: Matcher,
    reads_path: Path,
    target_start: int,
    target_length: int,
    offset_window: int,
    sample_limit: int,
) -> list[TriageFinding]:
    sampled: list[str] = []
    for record in iter_fastq(reads_path):
        sampled.append(record.seq)
        if len(sampled) >= sample_limit:
            break
    if not sampled:
        return []

    n = len(sampled)
    current_windows = [
        seq[target_start : target_start + target_length]
        for seq in sampled
        if target_start + target_length <= len(seq)
    ]
    current_hits = (
        sum(1 for r in matcher.assign_exact(current_windows) if r.status == MATCH_UNIQUE)
        if current_windows
        else 0
    )
    current_rate = _rate(current_hits, n)

    findings: list[TriageFinding] = []

    best_alt_offset: int | None = None
    best_alt_hits = 0
    low = max(0, target_start - offset_window)
    high = target_start + offset_window
    for candidate in range(low, high + 1):
        if candidate == target_start:
            continue
        windows = [
            seq[candidate : candidate + target_length]
            for seq in sampled
            if candidate + target_length <= len(seq)
        ]
        if not windows:
            continue
        hits = sum(1 for r in matcher.assign_exact(windows) if r.status == MATCH_UNIQUE)
        if hits > best_alt_hits:
            best_alt_hits = hits
            best_alt_offset = candidate

    best_alt_rate = _rate(best_alt_hits, n)
    if best_alt_offset is not None and best_alt_rate >= 0.25 and best_alt_rate > current_rate + 0.15:
        findings.append(
            TriageFinding(
                code="offset_shift_detected",
                severity="critical",
                summary=(
                    f"Exact match rate rises from {current_rate:.1%} at target_start={target_start} "
                    f"to {best_alt_rate:.1%} at target_start={best_alt_offset} in {n} sampled reads."
                ),
                evidence={
                    "configured_target_start": target_start,
                    "suggested_target_start": best_alt_offset,
                    "configured_exact_rate": round(current_rate, 6),
                    "suggested_exact_rate": round(best_alt_rate, 6),
                    "sampled_reads": n,
                },
                recommended_action=f"Re-run with --target-start {best_alt_offset} or inspect primer/spacer offsets.",
            )
        )

    if current_windows:
        rc_windows = [_reverse_complement(window) for window in current_windows]
        rc_hits = sum(1 for r in matcher.assign_exact(rc_windows) if r.status == MATCH_UNIQUE)
        rc_rate = _rate(rc_hits, n)
        if rc_rate >= 0.25 and rc_rate > current_rate + 0.15:
            findings.append(
                TriageFinding(
                    code="reverse_complement_detected",
                    severity="critical",
                    summary=(
                        f"Reverse-complementing windows at target_start={target_start} raises exact match rate "
                        f"from {current_rate:.1%} to {rc_rate:.1%} in {n} sampled reads."
                    ),
                    evidence={
                        "target_start": target_start,
                        "configured_exact_rate": round(current_rate, 6),
                        "reverse_complement_exact_rate": round(rc_rate, 6),
                        "sampled_reads": n,
                    },
                    recommended_action="Reverse-complement the target library or confirm read strand orientation.",
                )
            )

    return findings


def triage_run(
    *,
    targets: str | Path,
    reads: str | Path,
    target_start: int,
    target_length: int,
    out_dir: str | Path,
    sample_label: str = "sample",
    offset_window: int = 8,
    sample_limit: int = 10000,
) -> TriageReport:
    """Run deterministic library, policy-sensitivity, and failure-mode triage."""
    if offset_window < 0:
        raise ValueError("offset_window must be non-negative")
    if sample_limit < 1:
        raise ValueError("sample_limit must be positive")

    output = Path(out_dir)
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f"out_dir must be a new or empty directory: {output}")
    output.mkdir(parents=True, exist_ok=True)

    sensitivity_dir = output / "sensitivity"
    sensitivity = run_sensitivity(
        targets=targets,
        reads=reads,
        target_start=target_start,
        target_length=target_length,
        out_dir=sensitivity_dir,
        sample_label=sample_label,
    )

    read_count = int(sensitivity["read_count"])
    target_count = int(sensitivity["target_count"])
    outcomes = sensitivity["outcomes"]
    exact_rate = _rate(int(outcomes["exact"]["unique"]), read_count)
    radius_rate = _rate(int(outcomes["radius_k1"]["unique"]), read_count)
    best_rate = _rate(int(outcomes["best_k1"]["unique"]), read_count)
    invalid_reads = int(outcomes["exact"]["invalid"])
    invalid_rate = _rate(invalid_reads, read_count)
    changed_reads = int(sensitivity["changed_reads"])

    found: list[TriageFinding] = []
    target_rows = read_target_table(Path(targets))
    sequences = [row.sequence for row in target_rows]
    with Matcher(sequences) as matcher:
        found.extend(_audit_library_safety(matcher, sequences))
        if int(outcomes["exact"]["none"]) > 0:
            found.extend(
                _diagnose_read_orientation_and_offset(
                    matcher=matcher,
                    reads_path=Path(reads),
                    target_start=target_start,
                    target_length=target_length,
                    offset_window=offset_window,
                    sample_limit=sample_limit,
                )
            )

    if read_count == 0:
        found.append(
            TriageFinding(
                code="empty_fastq",
                severity="critical",
                summary="FASTQ input contains zero reads.",
                evidence={"read_count": 0},
                recommended_action="Verify the FASTQ file path and upstream demultiplexing output.",
            )
        )
    elif invalid_rate >= 0.10:
        found.append(
            TriageFinding(
                code="high_invalid_window_rate",
                severity="critical",
                summary=(
                    f"{invalid_reads} of {read_count} reads ({invalid_rate:.1%}) are too short for "
                    f"target_start={target_start}, target_length={target_length}."
                ),
                evidence={
                    "invalid_reads": invalid_reads,
                    "read_count": read_count,
                    "invalid_rate": round(invalid_rate, 6),
                    "required_min_read_length": target_start + target_length,
                },
                recommended_action="Check read trimming and verify target_start and target_length.",
            )
        )

    differing_pairs = [
        {
            "left": pair["left"],
            "right": pair["right"],
            "differing_guides": pair["differing_guides"],
            "total_count_delta": pair["total_count_delta"],
        }
        for pair in sensitivity["count_comparisons"]
        if not pair["counts_identical"]
    ]
    if changed_reads > 0 or differing_pairs:
        found.append(
            TriageFinding(
                code="policy_sensitive_counts",
                severity="warning",
                summary=(
                    f"{changed_reads} of {read_count} reads ({_rate(changed_reads, read_count):.1%}) "
                    "change assignment status or target across exact, radius k=1, and best k=1 policies."
                ),
                evidence={
                    "changed_reads": changed_reads,
                    "changed_rate": round(_rate(changed_reads, read_count), 6),
                    "differing_policy_pairs": differing_pairs,
                },
                recommended_action="Review sensitivity/report.html before fixing a mismatch policy for downstream inference.",
            )
        )

    diagnosed_codes = {finding.code for finding in found}
    if (
        read_count > 0
        and invalid_rate < 0.10
        and exact_rate < 0.50
        and "offset_shift_detected" not in diagnosed_codes
        and "reverse_complement_detected" not in diagnosed_codes
    ):
        found.append(
            TriageFinding(
                code="low_exact_assignment_rate",
                severity="warning",
                summary=f"Exact assignment rate is {exact_rate:.1%} with no offset-shift or strand explanation found.",
                evidence={
                    "exact_assignment_rate": round(exact_rate, 6),
                    "radius_k1_assignment_rate": round(radius_rate, 6),
                    "best_k1_assignment_rate": round(best_rate, 6),
                },
                recommended_action="Inspect unmatched reads for adapter contamination, wrong library, or quality drop.",
            )
        )

    findings = tuple(found)
    report_path = output / "triage_report.json"
    report = TriageReport(
        schema_version=SCHEMA_VERSION,
        verdict=_derive_verdict(findings),
        read_count=read_count,
        target_count=target_count,
        parameters={
            "target_start": target_start,
            "target_length": target_length,
            "offset_window": offset_window,
            "sample_limit": sample_limit,
            "sample_label": sample_label,
        },
        exact_assignment_rate=exact_rate,
        radius_k1_assignment_rate=radius_rate,
        best_k1_assignment_rate=best_rate,
        changed_reads=changed_reads,
        findings=findings,
        artifacts={
            "triage_report": str(report_path),
            "sensitivity_summary": str(sensitivity_dir / "summary.json"),
            "sensitivity_report": str(sensitivity_dir / "report.html"),
        },
    )
    report_path.write_text(json.dumps(report.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report
