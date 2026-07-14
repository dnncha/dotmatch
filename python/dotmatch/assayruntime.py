"""Experimental streaming runtime for compiled AssayScript v2 plans.

The runtime is deliberately fail-closed: it verifies compiled input
fingerprints, synchronizes all supplied FASTQs by read id, preserves segment
ambiguity, and writes outputs only after a complete successful run.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import os
import tempfile
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterator, Mapping, Sequence

from .assayscript import AssayScriptError, CompiledAssay, CompiledSegment
from .core import FastqRecord, iter_fastq, load_targets


STATUSES = {"unique", "ambiguous", "none", "invalid", "missing"}


@dataclass(frozen=True)
class SegmentCall:
    status: str
    reason: str
    observed: str
    target_id: str | None
    target_sequence: str | None
    best_distance: int | None
    candidates: tuple[str, ...]
    extraction_start: int | None
    orientation: str | None


@dataclass(frozen=True)
class JointCall:
    status: str
    reason: str
    combination: Mapping[str, str] | None
    compatible_combinations: int


@dataclass(frozen=True)
class RuntimeResult:
    output_dir: Path
    assignments: Path
    counts: Path
    events: Path
    summary: Path
    total_reads: int
    status_counts: Mapping[str, int]


@dataclass(frozen=True)
class _RuntimeSegment:
    compiled: CompiledSegment
    target_ids: tuple[str, ...]
    target_sequences: tuple[str, ...]


def load_compiled_plan(path: str | Path) -> CompiledAssay:
    plan_path = Path(path).expanduser().resolve()
    try:
        data = json.loads(plan_path.read_text(encoding="utf-8"))
        if data.get("compiler_schema_version") != 1:
            raise ValueError("unsupported compiler_schema_version")
        segments = []
        for item in data["segments"]:
            normalized = dict(item)
            if "length" not in normalized:
                target_lengths = normalized.get("target_lengths") or []
                if len(target_lengths) != 1:
                    raise ValueError(f"segment {normalized.get('name')} has no executable length")
                normalized["length"] = target_lengths[0]
            segments.append(CompiledSegment(**normalized))
        return CompiledAssay(
            compiler_schema_version=data["compiler_schema_version"],
            source_schema_version=data["source_schema_version"],
            name=data["name"],
            assay_type=data["assay_type"],
            source=data["source"],
            source_sha256=data["source_sha256"],
            segments=segments,
            allowed_combinations=data.get("allowed_combinations"),
            allowed_combinations_sha256=data.get("allowed_combinations_sha256"),
            execution_order=list(data["execution_order"]),
            findings=list(data.get("findings", [])),
        )
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise AssayScriptError(f"invalid compiled AssayScript plan {plan_path}: {exc}") from exc


def verify_plan_inputs(plan: CompiledAssay) -> dict[str, str]:
    fingerprints: dict[str, str] = {}
    _verify_fingerprint(Path(plan.source), plan.source_sha256, "AssayScript source")
    fingerprints[plan.source] = plan.source_sha256
    for segment in plan.segments:
        _verify_fingerprint(Path(segment.library), segment.library_sha256, f"{segment.name} library")
        fingerprints[segment.library] = segment.library_sha256
    if plan.allowed_combinations is not None:
        if not plan.allowed_combinations_sha256:
            raise AssayScriptError("compiled plan omits allowed-combinations fingerprint")
        _verify_fingerprint(
            Path(plan.allowed_combinations),
            plan.allowed_combinations_sha256,
            "allowed combinations",
        )
        fingerprints[plan.allowed_combinations] = plan.allowed_combinations_sha256
    return fingerprints


def run_compiled_plan(
    plan_path: str | Path,
    read_paths: Mapping[str, str | Path],
    output_dir: str | Path,
    *,
    max_reads: int | None = None,
    accept_findings: bool = False,
) -> RuntimeResult:
    plan_file = Path(plan_path).expanduser().resolve()
    plan = load_compiled_plan(plan_file)
    input_fingerprints = verify_plan_inputs(plan)
    if plan.findings and not accept_findings:
        raise AssayScriptError(
            "compiled plan has review findings; inspect them and explicitly accept_findings to execute"
        )
    if max_reads is not None and (
        isinstance(max_reads, bool) or not isinstance(max_reads, int) or max_reads <= 0
    ):
        raise ValueError("max_reads must be a positive integer")

    normalized_reads = {name.upper(): Path(path).expanduser().resolve() for name, path in read_paths.items()}
    unknown = sorted(set(normalized_reads) - {segment.read for segment in plan.segments})
    if unknown:
        raise AssayScriptError(f"read inputs are not used by the plan: {unknown}")
    required_reads = {segment.read for segment in plan.segments if segment.required}
    missing = sorted(required_reads - set(normalized_reads))
    if missing:
        raise AssayScriptError(f"missing FASTQ inputs for required reads: {missing}")
    for read_name, path in normalized_reads.items():
        if not path.is_file():
            raise AssayScriptError(f"{read_name} FASTQ does not exist: {path}")
        input_fingerprints[str(path)] = _sha256(path)

    order = {name: index for index, name in enumerate(plan.execution_order)}
    runtime_segments = [
        _load_runtime_segment(segment)
        for segment in sorted(plan.segments, key=lambda segment: order.get(segment.name, len(order)))
    ]
    combinations = _load_allowed_combinations(plan, runtime_segments)
    output = Path(output_dir).expanduser().resolve()
    if output.exists():
        if not output.is_dir():
            raise AssayScriptError(f"output path is not a directory: {output}")
        if any(output.iterdir()):
            raise AssayScriptError(f"output directory is not empty: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f".{output.name}-", dir=output.parent) as temporary:
        temp = Path(temporary)
        assignments_path = temp / "assignments.tsv"
        counts_path = temp / "counts.tsv"
        events_path = temp / "events.jsonl"
        summary_path = temp / "summary.json"
        status_counts: Counter[str] = Counter()
        segment_counts = {segment.compiled.name: Counter() for segment in runtime_segments}
        combination_counts: Counter[tuple[tuple[str, str], ...]] = Counter()
        total_reads = 0

        with assignments_path.open("w", encoding="utf-8", newline="") as assignments_handle, events_path.open(
            "w", encoding="utf-8"
        ) as events_handle:
            writer = csv.DictWriter(
                assignments_handle,
                fieldnames=_assignment_fields(runtime_segments),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for read_id, records in _synchronized_records(normalized_reads):
                if max_reads is not None and total_reads >= max_reads:
                    break
                calls = {
                    segment.compiled.name: _call_segment(segment, records.get(segment.compiled.read))
                    for segment in runtime_segments
                }
                joint = _joint_call(runtime_segments, calls, combinations)
                total_reads += 1
                status_counts[joint.status] += 1
                for name, call in calls.items():
                    segment_counts[name][call.status] += 1
                key = _combination_key(joint.combination)
                if joint.status == "unique" and joint.combination is not None:
                    combination_counts[tuple(sorted(joint.combination.items()))] += 1
                writer.writerow(_assignment_row(read_id, joint, calls, runtime_segments))
                events_handle.write(
                    json.dumps(
                        {"read_id": read_id, "status": joint.status, "target": key},
                        sort_keys=True,
                    )
                    + "\n"
                )

        with counts_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["combination", "count"])
            for items, count in sorted(combination_counts.items()):
                writer.writerow([_combination_key(dict(items)), count])

        summary = {
            "runtime_schema_version": 1,
            "status": "experimental",
            "assay": plan.name,
            "assay_type": plan.assay_type,
            "compiled_plan": str(plan_file),
            "compiled_plan_sha256": _sha256(plan_file),
            "inputs": dict(sorted(input_fingerprints.items())),
            "total_reads": total_reads,
            "status_counts": {status: status_counts.get(status, 0) for status in ("unique", "ambiguous", "none", "invalid")},
            "rates": {
                status: (status_counts.get(status, 0) / total_reads if total_reads else 0.0)
                for status in ("unique", "ambiguous", "none", "invalid")
            },
            "segment_status_counts": {
                name: {status: counter.get(status, 0) for status in sorted(STATUSES)}
                for name, counter in segment_counts.items()
            },
            "unique_combinations": len(combination_counts),
            "allowed_combinations": len(combinations) if combinations is not None else None,
            "plan_findings": plan.findings,
        }
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        output.mkdir(parents=True, exist_ok=True)
        final_paths = {}
        for name, source in {
            "assignments": assignments_path,
            "counts": counts_path,
            "events": events_path,
            "summary": summary_path,
        }.items():
            destination = output / source.name
            os.replace(source, destination)
            final_paths[name] = destination

    return RuntimeResult(
        output_dir=output,
        assignments=final_paths["assignments"],
        counts=final_paths["counts"],
        events=final_paths["events"],
        summary=final_paths["summary"],
        total_reads=total_reads,
        status_counts=dict(status_counts),
    )


def _load_runtime_segment(segment: CompiledSegment) -> _RuntimeSegment:
    targets = load_targets(segment.library)
    identifiers = tuple(target_id for target_id, _sequence in targets)
    if len(set(identifiers)) != len(identifiers):
        raise AssayScriptError(f"{segment.name} library contains duplicate target identifiers")
    return _RuntimeSegment(
        compiled=segment,
        target_ids=identifiers,
        target_sequences=tuple(sequence.upper() for _target_id, sequence in targets),
    )


def _synchronized_records(read_paths: Mapping[str, Path]) -> Iterator[tuple[str, Mapping[str, FastqRecord]]]:
    if not read_paths:
        raise AssayScriptError("at least one FASTQ input is required")
    names = sorted(read_paths)
    iterators = [iter_fastq(read_paths[name]) for name in names]
    sentinel = object()
    for index, group in enumerate(itertools.zip_longest(*iterators, fillvalue=sentinel), start=1):
        if sentinel in group:
            exhausted = [name for name, record in zip(names, group) if record is sentinel]
            raise AssayScriptError(
                f"FASTQ record counts differ at record {index}; exhausted inputs: {exhausted}"
            )
        records = {name: record for name, record in zip(names, group) if isinstance(record, FastqRecord)}
        canonical = {_canonical_read_id(record.read_id) for record in records.values()}
        if len(canonical) != 1:
            observed = {name: record.read_id for name, record in records.items()}
            raise AssayScriptError(f"FASTQ read ids are not synchronized at record {index}: {observed}")
        yield next(iter(canonical)), records


def _call_segment(segment: _RuntimeSegment, record: FastqRecord | None) -> SegmentCall:
    compiled = segment.compiled
    if record is None:
        return SegmentCall("missing", "read_not_supplied", "", None, None, None, (), None, None)
    windows, absent_reason = _candidate_windows(record.seq, compiled)
    if not windows:
        status = "none" if absent_reason == "anchor_not_found" else "invalid"
        return SegmentCall(status, absent_reason, "", None, None, None, (), None, None)

    evidence: dict[str, tuple[int, str, str, int, str]] = {}
    for window, start, orientation in windows:
        for target_id, target_sequence in zip(segment.target_ids, segment.target_sequences):
            distance = _distance(window, target_sequence, compiled.metric, compiled.k)
            if distance <= compiled.k:
                candidate = (distance, target_id, window, start, orientation)
                previous = evidence.get(target_id)
                if previous is None or candidate < previous:
                    evidence[target_id] = candidate
    if not evidence:
        first = min(windows, key=lambda item: (item[1], item[2], item[0]))
        return SegmentCall("none", "no_target_within_radius", first[0], None, None, None, (), first[1], first[2])

    candidates = sorted(evidence)
    if compiled.ambiguity_policy == "best":
        best_distance = min(item[0] for item in evidence.values())
        candidates = sorted(target for target, item in evidence.items() if item[0] == best_distance)
    ranked = sorted((evidence[target] for target in candidates), key=lambda item: item)
    best = ranked[0]
    if len(candidates) == 1:
        index = segment.target_ids.index(candidates[0])
        return SegmentCall(
            "unique",
            "single_compatible_target",
            best[2],
            candidates[0],
            segment.target_sequences[index],
            best[0],
            tuple(candidates),
            best[3],
            best[4],
        )
    return SegmentCall(
        "ambiguous",
        "multiple_compatible_targets",
        best[2],
        None,
        None,
        best[0],
        tuple(candidates),
        best[3],
        best[4],
    )


def _candidate_windows(sequence: str, segment: CompiledSegment) -> tuple[list[tuple[str, int, str]], str]:
    oriented = [("forward", sequence.upper())]
    if segment.orientation == "reverse_complement":
        oriented = [("reverse_complement", _reverse_complement(sequence))]
    elif segment.orientation == "auto":
        oriented.append(("reverse_complement", _reverse_complement(sequence)))

    windows: set[tuple[str, int, str]] = set()
    saw_anchor = False
    for orientation, read in oriented:
        bases: list[int] = []
        if segment.anchor is None:
            assert segment.start is not None
            bases = [segment.start]
        else:
            offset = 0
            while True:
                position = read.find(segment.anchor, offset)
                if position < 0:
                    break
                saw_anchor = True
                bases.append(position + len(segment.anchor) + segment.after_anchor)
                offset = position + 1
        for base in bases:
            for start in range(base - segment.jitter, base + segment.jitter + 1):
                end = start + segment.length
                if start >= 0 and end <= len(read):
                    windows.add((read[start:end], start, orientation))
    if windows:
        return sorted(windows, key=lambda item: (item[1], item[2], item[0])), ""
    if segment.anchor is not None and not saw_anchor:
        return [], "anchor_not_found"
    return [], "window_out_of_bounds"


def _joint_call(
    segments: Sequence[_RuntimeSegment],
    calls: Mapping[str, SegmentCall],
    combinations: list[dict[str, str]] | None,
) -> JointCall:
    required = [segment.compiled for segment in segments if segment.compiled.required]
    if any(calls[segment.name].status == "invalid" for segment in required):
        return JointCall("invalid", "required_segment_invalid", None, 0)
    if any(not calls[segment.name].candidates for segment in required):
        return JointCall("none", "required_segment_unassigned", None, 0)

    if combinations is None:
        if any(len(calls[segment.name].candidates) > 1 for segment in required):
            return JointCall("ambiguous", "required_segment_ambiguous", None, 0)
        combination = {
            segment.compiled.name: calls[segment.compiled.name].candidates[0]
            for segment in segments
            if len(calls[segment.compiled.name].candidates) == 1
        }
        return JointCall("unique", "all_required_segments_unique", combination, 1)

    constrained = set(combinations[0]) if combinations else set()
    required_by_name = {segment.compiled.name: segment.compiled.required for segment in segments}
    compatible = []
    for combination in combinations:
        if all(
            segment_name not in calls
            or (calls[segment_name].status == "missing" and not required_by_name[segment_name])
            or combination[segment_name] in calls[segment_name].candidates
            for segment_name in constrained
        ):
            compatible.append(combination)
    if not compatible:
        return JointCall("none", "disallowed_combination", None, 0)
    if len(compatible) > 1:
        return JointCall("ambiguous", "multiple_allowed_combinations", None, len(compatible))

    resolved = dict(compatible[0])
    for segment in segments:
        name = segment.compiled.name
        if name not in resolved and len(calls[name].candidates) == 1:
            resolved[name] = calls[name].candidates[0]
        elif name not in resolved and segment.compiled.required:
            return JointCall("ambiguous", "unconstrained_required_segment_ambiguous", None, 1)
    return JointCall("unique", "allowed_combination_rescue", resolved, 1)


def _load_allowed_combinations(
    plan: CompiledAssay,
    segments: Sequence[_RuntimeSegment],
) -> list[dict[str, str]] | None:
    if plan.allowed_combinations is None:
        return None
    path = Path(plan.allowed_combinations)
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise AssayScriptError("allowed combinations have no header")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise AssayScriptError("allowed combinations contain duplicate segment columns")
        target_ids = {segment.compiled.name: set(segment.target_ids) for segment in segments}
        combinations = []
        for row_index, row in enumerate(reader, start=2):
            normalized = {name: str(value or "").strip() for name, value in row.items()}
            if any(not value for value in normalized.values()):
                raise AssayScriptError(f"allowed combinations row {row_index} contains an empty target")
            for segment_name, target_id in normalized.items():
                if segment_name not in target_ids:
                    raise AssayScriptError(
                        f"allowed combinations row {row_index} references unknown segment {segment_name}"
                    )
                if target_id not in target_ids[segment_name]:
                    raise AssayScriptError(
                        f"allowed combinations row {row_index} references unknown "
                        f"{segment_name} target {target_id}"
                    )
            combinations.append(normalized)
    return combinations


def _distance(left: str, right: str, metric: str, limit: int) -> int:
    if metric == "hamming":
        if len(left) != len(right):
            return limit + 1
        return sum(a != b for a, b in zip(left, right))
    if abs(len(left) - len(right)) > limit:
        return limit + 1
    previous = list(range(len(right) + 1))
    for i, left_base in enumerate(left, start=1):
        current = [i]
        row_min = i
        for j, right_base in enumerate(right, start=1):
            value = min(current[j - 1] + 1, previous[j] + 1, previous[j - 1] + (left_base != right_base))
            current.append(value)
            row_min = min(row_min, value)
        if row_min > limit:
            return limit + 1
        previous = current
    return previous[-1]


def _assignment_fields(segments: Sequence[_RuntimeSegment]) -> list[str]:
    fields = ["read_id", "status", "reason", "combination", "compatible_combinations"]
    for segment in segments:
        name = segment.compiled.name
        fields.extend(
            [
                f"{name}_status",
                f"{name}_reason",
                f"{name}_target",
                f"{name}_observed",
                f"{name}_distance",
                f"{name}_candidates",
                f"{name}_start",
                f"{name}_orientation",
            ]
        )
    return fields


def _assignment_row(
    read_id: str,
    joint: JointCall,
    calls: Mapping[str, SegmentCall],
    segments: Sequence[_RuntimeSegment],
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "read_id": read_id,
        "status": joint.status,
        "reason": joint.reason,
        "combination": _combination_key(joint.combination),
        "compatible_combinations": joint.compatible_combinations,
    }
    for segment in segments:
        name = segment.compiled.name
        call = calls[name]
        resolved_target = joint.combination.get(name) if joint.combination else call.target_id
        row.update(
            {
                f"{name}_status": call.status,
                f"{name}_reason": call.reason,
                f"{name}_target": resolved_target or "",
                f"{name}_observed": call.observed,
                f"{name}_distance": "" if call.best_distance is None else call.best_distance,
                f"{name}_candidates": ";".join(call.candidates),
                f"{name}_start": "" if call.extraction_start is None else call.extraction_start,
                f"{name}_orientation": call.orientation or "",
            }
        )
    return row


def _combination_key(combination: Mapping[str, str] | None) -> str:
    if not combination:
        return ""
    return "|".join(f"{name}={target}" for name, target in sorted(combination.items()))


def _canonical_read_id(read_id: str) -> str:
    return read_id[:-2] if read_id.endswith(("/1", "/2")) else read_id


def _reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _verify_fingerprint(path: Path, expected: str, label: str) -> None:
    if not path.is_file():
        raise AssayScriptError(f"{label} does not exist: {path}")
    observed = _sha256(path)
    if observed != expected:
        raise AssayScriptError(f"{label} fingerprint changed: expected {expected}, observed {observed}")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
