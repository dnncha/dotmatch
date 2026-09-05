"""Bounded JSON I/O, canonical checksums and explicit, atomic output replacement."""
from __future__ import annotations

import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeVar

from pydantic import BaseModel

if TYPE_CHECKING:
    from .models import Analysis, Manifest, ScanResult

MAX_INPUT_BYTES = 8 * 1024 * 1024
MAX_RESULT_BYTES = 64 * 1024 * 1024
M = TypeVar("M", bound=BaseModel)


class InputError(ValueError):
    pass


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise InputError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    raise InputError(f"nonfinite JSON number is not supported: {value}")


def read_json(path: str | Path, *, max_bytes: int = MAX_INPUT_BYTES) -> Any:
    if str(path) == "-":
        data = sys.stdin.buffer.read(max_bytes + 1)
    else:
        with Path(path).open("rb") as handle:
            data = handle.read(max_bytes + 1)
    if len(data) > max_bytes:
        raise InputError(f"input exceeds {max_bytes} byte limit")
    try:
        return json.loads(data.decode("utf-8"), object_pairs_hook=_unique_object, parse_constant=_reject_constant)
    except InputError:
        raise
    except (ValueError, UnicodeDecodeError, RecursionError) as error:
        raise InputError("invalid or excessively nested UTF-8 JSON") from error


def load_manifest(path: str | Path) -> Manifest:
    """Load strictly typed JSON. Unknown fields, duplicate keys and ambiguous coordinates fail."""
    from .models import Manifest
    return Manifest.model_validate(read_json(path))


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode()).hexdigest()


def seal(result: M) -> M:
    data = result.model_dump(mode="json", exclude={"result_sha256"})
    return result.model_copy(update={"result_sha256": digest(data)})


def verify_result(path: str | Path) -> Analysis | ScanResult:
    from .models import Analysis, ScanResult
    data = read_json(path, max_bytes=MAX_RESULT_BYTES)
    if not isinstance(data, dict):
        raise InputError("result must be a JSON object")
    result: Analysis | ScanResult
    kind = data.get("kind")
    if kind == "editwitness.analysis":
        result = Analysis.model_validate(data)
    elif kind == "editwitness.deletion_scan":
        result = ScanResult.model_validate(data)
    else:
        raise InputError("unsupported result kind; compact summaries cannot be verified as full results")
    if not result.result_sha256 or seal(result).result_sha256 != result.result_sha256:
        raise InputError("result checksum mismatch")
    return result


def check_destinations(paths: list[str | Path | None], *, force: bool = False,
                       input_path: str | Path | None = None) -> None:
    files = [Path(p) for p in paths if p is not None and str(p) != "-"]
    resolved = [p.resolve() for p in files]
    if len(resolved) != len(set(resolved)):
        raise InputError("output destinations must be different files")
    if input_path is not None and str(input_path) != "-" and Path(input_path).resolve() in resolved:
        raise InputError("refusing to replace the input file, even with --force")
    for path in files:
        if not path.parent.is_dir():
            raise FileNotFoundError(f"output directory does not exist: {path.parent}")
        if path.exists() and (not force or path.is_dir()):
            raise FileExistsError(f"output exists: {path}; use --force to replace a file")


def atomic_write(path: str | Path, text: str, *, force: bool = False) -> None:
    """Write one file atomically. No silent overwrite; permissions default to owner-only."""
    path = Path(path)
    descriptor, name = tempfile.mkstemp(prefix=".editwitness-", dir=path.parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        if force:
            os.replace(temporary, path)
        else:
            os.link(temporary, path)  # Atomic create-if-absent, including cross-process races.
    finally:
        temporary.unlink(missing_ok=True)
