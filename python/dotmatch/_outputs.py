"""Private staging for count artifacts and streamed assignment files."""
from __future__ import annotations

import os
import stat
import tempfile
import warnings
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Mapping, Sequence


def _identity(path: Path) -> tuple[int, int] | None:
    try:
        info = path.lstat()
    except FileNotFoundError:
        return None
    return info.st_dev, info.st_ino


def _version(path: Path) -> tuple[int, ...] | None:
    try:
        info = path.lstat()
    except FileNotFoundError:
        return None
    if not stat.S_ISREG(info.st_mode):
        raise ValueError(f"output must be a regular file, not a symlink or special file: {path}")
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns


def _sibling(path: Path, suffix: str = '') -> Path:
    descriptor, name = tempfile.mkstemp(prefix='.dotmatch-', suffix=suffix, dir=path.parent)
    os.close(descriptor)
    return Path(name)


def _unlink_owned(path: Path, identity: tuple[int, int] | None) -> bool:
    current = _identity(path)
    if current is None:
        return True
    if current != identity:
        return False
    path.unlink()
    return True


@dataclass
class _Output:
    final: Path
    original: tuple[int, ...] | None
    temporary: Path | None = None
    temporary_identity: tuple[int, int] | None = None
    backup: Path | None = None
    keep_backup: bool = False


def _check_unchanged(entry: _Output) -> None:
    if _version(entry.final) != entry.original:
        raise FileExistsError(f"output changed while counting; refusing to replace it: {entry.final}")


@contextmanager
def staged_outputs(destinations: Mapping[str, str | Path], *,
                   inputs: Sequence[str | Path] = ()) -> Iterator[dict[str, Path]]:
    """Publish closed artifacts after success; restore old files if replacement fails.

    Individual replacements are atomic. The set is not atomic across a crash.
    """
    input_ids = set()
    for source in inputs:
        info = Path(source).stat()
        if not stat.S_ISREG(info.st_mode):
            raise ValueError(f"input must be a regular file: {source}")
        input_ids.add((info.st_dev, info.st_ino))
    entries: dict[str, _Output] = {}
    paths, output_ids = set(), set()
    for role, destination in destinations.items():
        path = Path(destination).absolute()
        path = path.parent.resolve(strict=True) / path.name
        original = _version(path)
        identity = original[:2] if original else None
        if identity in input_ids:
            raise ValueError(f"output aliases an input file: {destination}")
        if path in paths or identity is not None and identity in output_ids:
            raise ValueError(f"output paths alias each other: {destination}")
        paths.add(path)
        if identity is not None:
            output_ids.add(identity)
        entries[role] = _Output(path, original)
    published: list[_Output] = []
    try:
        for entry in entries.values():
            entry.temporary = _sibling(entry.final, entry.final.suffix)
            entry.temporary_identity = _identity(entry.temporary)
        yield {role: entry.temporary for role, entry in entries.items()}
        for entry in entries.values():
            _check_unchanged(entry)
        for entry in entries.values():
            _check_unchanged(entry)
            if entry.original is not None:
                entry.backup = _sibling(entry.final)
                entry.backup.unlink()
                os.link(entry.final, entry.backup)
                linked = _version(entry.final)
                if linked is None or linked[:4] != entry.original[:4]:
                    raise FileExistsError(f"output changed during backup: {entry.final}")
                entry.original = linked
        for entry in entries.values():
            _check_unchanged(entry)
            if _identity(entry.temporary) != entry.temporary_identity:
                raise OSError(f"staged output changed before publication: {entry.temporary}")
            os.replace(entry.temporary, entry.final)
            published.append(entry)
    except BaseException as exc:
        recovery_errors = []
        for entry in reversed(published):
            try:
                if _identity(entry.final) != entry.temporary_identity:
                    raise OSError(f"published output changed during recovery: {entry.final}")
                if entry.backup is not None:
                    os.replace(entry.backup, entry.final)
                else:
                    entry.final.unlink()
            except OSError as recovery:
                entry.keep_backup = True
                retained = f"; previous file retained at {entry.backup}" if entry.backup else ''
                recovery_errors.append(f"{recovery}{retained}")
        if recovery_errors:
            raise OSError('; '.join(recovery_errors)) from exc
        raise
    finally:
        for entry in entries.values():
            owned = [(entry.temporary, entry.temporary_identity)]
            if entry.backup is not None and not entry.keep_backup:
                owned.append((entry.backup, entry.original[:2]))
            for path, identity in owned:
                if path is not None:
                    try:
                        _unlink_owned(path, identity)
                    except OSError as exc:
                        warnings.warn(f"could not remove temporary output {path}: {exc}", RuntimeWarning)
