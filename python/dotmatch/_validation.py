"""Small, dependency-free validators for the public Python/native boundary."""
from __future__ import annotations

import ctypes
import operator

C_INT_MAX = (1 << (ctypes.sizeof(ctypes.c_int) * 8 - 1)) - 1


def integer(value, name: str, *, minimum: int = 0, maximum: int = C_INT_MAX) -> int:
    """Accept integer scalars, never booleans, rounding, or C integer wrapping."""
    if isinstance(value, bool) or type(value).__name__ == "bool_":
        raise TypeError(f"{name} must be an integer, not a boolean")
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc
    if not minimum <= result <= maximum:
        raise ValueError(f"{name} must be between {minimum} and {maximum}")
    return result


def identifier(value, name: str = "identifier") -> str:
    """Preserve textual identifiers (including NA and leading zeros)."""
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if isinstance(value, str):
        result = value.strip()
    else:
        # In-memory integer labels are useful, but None/NaN/float coercion is not.
        try:
            if isinstance(value, bool) or type(value).__name__ == "bool_":
                raise TypeError
            result = str(operator.index(value))
        except TypeError as exc:
            raise ValueError(f"{name} must be nonempty text or an integer label") from exc
    if not result or any(ord(ch) < 32 or ord(ch) == 127 for ch in result):
        raise ValueError(f"{name} must be nonempty and contain no control characters")
    return result


def sequence(value, name: str = "sequence") -> str:
    """Normalize high-level sequence inputs without stringifying missing data."""
    if isinstance(value, bytes):
        value = value.decode("ascii")
    if not isinstance(value, str) or not value or not value.isascii():
        raise ValueError(f"{name} must be nonempty ASCII text or bytes")
    if any(ord(ch) <= 32 or ord(ch) >= 127 for ch in value):
        raise ValueError(f"{name} must contain no whitespace or control characters")
    return value.upper()


def named_targets(rows) -> list[tuple[str, str]]:
    result, seen = [], set()
    for raw_id, raw_seq in rows:
        target_id = identifier(raw_id, "target ID")
        if target_id in seen:
            raise ValueError(f"duplicate target ID: {target_id!r}; use distinct IDs even for duplicate sequences")
        seen.add(target_id)
        result.append((target_id, sequence(raw_seq, f"sequence for {target_id}")))
    if not result:
        raise ValueError("targets must not be empty")
    return result
