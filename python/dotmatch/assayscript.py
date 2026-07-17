"""AssayScript v2 validation and deterministic compilation.

The existing AssaySpec schema_version=1 remains supported by dotmatch.assayspec.
This module introduces a multi-read, multi-segment representation without
changing current execution semantics. Compilation produces a portable JSON plan
that records strategy selection, input fingerprints, and safety findings.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore[no-redef]


READS = {"R1", "R2", "I1", "I2"}
ORIENTATIONS = {"forward", "reverse_complement", "auto"}
METRICS = {"hamming", "levenshtein"}
POLICIES = {"radius", "best"}
DNA = frozenset("ACGTN")


class AssayScriptError(ValueError):
    pass


@dataclass(frozen=True)
class Segment:
    name: str
    read: str
    library: Path
    length: int
    start: int | None
    anchor: str | None
    after_anchor: int
    jitter: int
    orientation: str
    metric: str
    k: int
    ambiguity_policy: str
    required: bool


@dataclass(frozen=True)
class CompiledSegment:
    name: str
    read: str
    library: str
    library_sha256: str
    target_count: int
    target_lengths: list[int]
    strategy: str
    start: int | None
    anchor: str | None
    after_anchor: int
    jitter: int
    orientation: str
    metric: str
    k: int
    ambiguity_policy: str
    required: bool
    safety_status: str
    findings: list[str]


@dataclass(frozen=True)
class CompiledAssay:
    compiler_schema_version: int
    source_schema_version: int
    name: str
    assay_type: str
    source: str
    source_sha256: str
    segments: list[CompiledSegment]
    allowed_combinations: str | None
    allowed_combinations_sha256: str | None
    execution_order: list[str]
    findings: list[str]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def load_and_compile(path: str | Path) -> CompiledAssay:
    source = Path(path).expanduser().resolve()
    try:
        data = tomllib.loads(source.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise AssayScriptError(f"AssayScript does not exist: {source}") from exc
    except tomllib.TOMLDecodeError as exc:
        raise AssayScriptError(f"{source}: invalid TOML: {exc}") from exc
    if not isinstance(data, dict):
        raise AssayScriptError("AssayScript top level must be a table")
    return compile_assayscript(data, source=source)


def compile_assayscript(data: Mapping[str, Any], *, source: Path) -> CompiledAssay:
    if data.get("schema_version") != 2:
        raise AssayScriptError("schema_version must be 2 for AssayScript")
    name = _required_text(data.get("name"), "name")
    assay_type = _required_text(data.get("assay_type"), "assay_type")
    raw_segments = data.get("segments")
    if not isinstance(raw_segments, list) or not raw_segments:
        raise AssayScriptError("at least one [[segments]] entry is required")

    segments: list[Segment] = []
    names: set[str] = set()
    for index, raw in enumerate(raw_segments):
        if not isinstance(raw, dict):
            raise AssayScriptError(f"segments[{index}] must be a table")
        segment = _parse_segment(raw, source=source, index=index)
        if segment.name in names:
            raise AssayScriptError(f"duplicate segment name: {segment.name}")
        names.add(segment.name)
        segments.append(segment)

    constraints = data.get("constraints", {})
    if not isinstance(constraints, dict):
        raise AssayScriptError("constraints must be a table")
    allowed_path: Path | None = None
    if constraints.get("allowed_combinations") is not None:
        allowed_path = _resolve_input(source, constraints["allowed_combinations"], "constraints.allowed_combinations")

    compiled_segments = [_compile_segment(segment) for segment in segments]
    findings = [finding for segment in compiled_segments for finding in segment.findings]
    if allowed_path is not None:
        findings.extend(_validate_combinations(allowed_path, names))

    execution_order = [
        segment.name
        for segment in sorted(
            compiled_segments,
            key=lambda item: (
                item.anchor is not None,
                item.jitter,
                item.k,
                item.target_count,
                item.name,
            ),
        )
    ]
    return CompiledAssay(
        compiler_schema_version=1,
        source_schema_version=2,
        name=name,
        assay_type=assay_type,
        source=str(source),
        source_sha256=_sha256(source),
        segments=compiled_segments,
        allowed_combinations=str(allowed_path) if allowed_path else None,
        allowed_combinations_sha256=_sha256(allowed_path) if allowed_path else None,
        execution_order=execution_order,
        findings=sorted(set(findings)),
    )


def write_compiled_plan(plan: CompiledAssay, path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(plan.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return output


def _parse_segment(raw: Mapping[str, Any], *, source: Path, index: int) -> Segment:
    prefix = f"segments[{index}]"
    name = _required_text(raw.get("name"), f"{prefix}.name")
    read = _enum(raw.get("read"), READS, f"{prefix}.read")
    library = _resolve_input(source, raw.get("library"), f"{prefix}.library")
    length = _integer(raw.get("length"), 1, 512, f"{prefix}.length")
    start_value = raw.get("start")
    anchor_value = raw.get("anchor")
    if (start_value is None) == (anchor_value is None):
        raise AssayScriptError(f"{prefix} requires exactly one of start or anchor")
    start = None if start_value is None else _integer(start_value, 0, 100000, f"{prefix}.start")
    anchor = None
    if anchor_value is not None:
        anchor = _required_text(anchor_value, f"{prefix}.anchor").upper()
        if any(base not in DNA for base in anchor):
            raise AssayScriptError(f"{prefix}.anchor contains unsupported symbols")
    after_anchor = _integer(raw.get("after_anchor", 0), -512, 512, f"{prefix}.after_anchor")
    jitter = _integer(raw.get("jitter", 0), 0, 32, f"{prefix}.jitter")
    orientation = _enum(raw.get("orientation", "forward"), ORIENTATIONS, f"{prefix}.orientation")
    metric = _enum(raw.get("metric", "hamming"), METRICS, f"{prefix}.metric")
    k = _integer(raw.get("k", 0), 0, 3, f"{prefix}.k")
    if metric == "levenshtein" and k > 2:
        raise AssayScriptError(f"{prefix}: levenshtein currently supports k <= 2")
    policy = _enum(raw.get("ambiguity_policy", "radius"), POLICIES, f"{prefix}.ambiguity_policy")
    required = raw.get("required", True)
    if not isinstance(required, bool):
        raise AssayScriptError(f"{prefix}.required must be true or false")
    return Segment(
        name=name,
        read=read,
        library=library,
        length=length,
        start=start,
        anchor=anchor,
        after_anchor=after_anchor,
        jitter=jitter,
        orientation=orientation,
        metric=metric,
        k=k,
        ambiguity_policy=policy,
        required=required,
    )


def _compile_segment(segment: Segment) -> CompiledSegment:
    targets = _read_library(segment.library)
    lengths = sorted({len(sequence) for sequence in targets})
    findings: list[str] = []
    if lengths != [segment.length]:
        findings.append(
            f"{segment.name}: declared length {segment.length} differs from target lengths {lengths}"
        )
    duplicates = len(targets) - len(set(targets))
    if duplicates:
        findings.append(f"{segment.name}: {duplicates} duplicate target sequences")
    safety_status = "safe"
    if duplicates:
        safety_status = "unsafe"
    elif segment.k > 0:
        close_pairs = _bounded_close_pairs(targets, segment.k, segment.metric)
        if close_pairs is None:
            safety_status = "not_computed"
            findings.append(
                f"{segment.name}: pairwise safety not computed above 5000 targets; run native audit"
            )
        elif close_pairs:
            safety_status = "unsafe"
            findings.append(
                f"{segment.name}: {close_pairs} target pairs overlap the configured correction radius"
            )

    strategy = _select_strategy(segment, targets)
    if segment.orientation == "auto":
        findings.append(f"{segment.name}: orientation auto requires runtime evidence and review")
    if segment.jitter:
        findings.append(f"{segment.name}: jitter search spans ±{segment.jitter} bases")
    return CompiledSegment(
        name=segment.name,
        read=segment.read,
        library=str(segment.library),
        library_sha256=_sha256(segment.library),
        target_count=len(targets),
        target_lengths=lengths,
        strategy=strategy,
        start=segment.start,
        anchor=segment.anchor,
        after_anchor=segment.after_anchor,
        jitter=segment.jitter,
        orientation=segment.orientation,
        metric=segment.metric,
        k=segment.k,
        ambiguity_policy=segment.ambiguity_policy,
        required=segment.required,
        safety_status=safety_status,
        findings=findings,
    )


def _select_strategy(segment: Segment, targets: Sequence[str]) -> str:
    uniform = len({len(target) for target in targets}) == 1
    if segment.k == 0 and uniform:
        base = "exact_hash"
    elif segment.metric == "hamming" and uniform and segment.k == 1 and segment.length <= 32:
        base = "packed_hamming_neighborhood"
    elif segment.metric == "hamming" and uniform:
        base = "seeded_hamming_verify"
    elif segment.metric == "levenshtein" and segment.length <= 32 and segment.k <= 2:
        base = "packed_levenshtein_neighborhood"
    else:
        base = "seeded_levenshtein_verify"
    if segment.anchor is not None:
        base = "anchor_scan+" + base
    if segment.jitter:
        base = "jitter_scan+" + base
    if segment.orientation != "forward":
        base = "orientation_dispatch+" + base
    return base


def _read_library(path: Path) -> list[str]:
    rows: list[str] = []
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        columns = [column.strip() for column in line.split(delimiter)]
        sequence = columns[1] if len(columns) > 1 else columns[0]
        normalized = sequence.upper()
        if normalized.lower() in {"sequence", "seq", "target_seq", "guide_seq", "barcode_seq"}:
            continue
        if not normalized:
            raise AssayScriptError(f"empty target sequence in {path}")
        rows.append(normalized)
    if not rows:
        raise AssayScriptError(f"no targets found in {path}")
    return rows


def _bounded_close_pairs(targets: Sequence[str], k: int, metric: str) -> int | None:
    if len(targets) > 5000:
        return None
    threshold = 2 * k
    count = 0
    for left_index, left in enumerate(targets):
        for right in targets[left_index + 1 :]:
            if metric == "hamming":
                if len(left) == len(right) and sum(a != b for a, b in zip(left, right)) <= threshold:
                    count += 1
            elif _levenshtein_bounded(left, right, threshold) <= threshold:
                count += 1
    return count


def _levenshtein_bounded(left: str, right: str, limit: int) -> int:
    if abs(len(left) - len(right)) > limit:
        return limit + 1
    previous = list(range(len(right) + 1))
    for i, left_base in enumerate(left, start=1):
        current = [i]
        row_min = i
        for j, right_base in enumerate(right, start=1):
            value = min(
                current[j - 1] + 1,
                previous[j] + 1,
                previous[j - 1] + (left_base != right_base),
            )
            current.append(value)
            row_min = min(row_min, value)
        if row_min > limit:
            return limit + 1
        previous = current
    return previous[-1]


def _validate_combinations(path: Path, names: set[str]) -> list[str]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise AssayScriptError("allowed combinations file is empty")
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    header = {column.strip() for column in lines[0].split(delimiter)}
    unknown = sorted(header - names)
    if unknown:
        raise AssayScriptError(f"allowed combinations reference unknown segments: {unknown}")
    missing = sorted(names - header)
    findings: list[str] = []
    if missing:
        findings.append(f"allowed combinations do not constrain optional segments: {missing}")
    findings.append(f"allowed combinations contain {max(0, len(lines) - 1)} tuples")
    return findings


def _resolve_input(source: Path, value: Any, name: str) -> Path:
    text = _required_text(value, name)
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = source.parent / path
    path = path.resolve()
    if not path.is_file():
        raise AssayScriptError(f"{name} does not exist: {path}")
    return path


def _required_text(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise AssayScriptError(f"{name} must be a non-empty string")
    return value.strip()


def _enum(value: Any, choices: set[str], name: str) -> str:
    text = _required_text(value, name)
    if text not in choices:
        raise AssayScriptError(f"{name} must be one of {sorted(choices)}")
    return text


def _integer(value: Any, minimum: int, maximum: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise AssayScriptError(f"{name} must be an integer")
    if value < minimum or value > maximum:
        raise AssayScriptError(f"{name} must be between {minimum} and {maximum}")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
