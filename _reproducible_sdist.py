"""Deterministic source-distribution archive support.

Setuptools does not currently normalize the release-tree or gzip timestamps
when SOURCE_DATE_EPOCH is set.  Keep this helper dependency-free because
setup.py imports it while building from both a checkout and an sdist.
"""
from __future__ import annotations

import gzip
import os
import tarfile
from pathlib import Path


_MAX_GZIP_MTIME = (1 << 32) - 1


def parse_source_date_epoch(value: str) -> int:
    try:
        epoch = int(value)
    except ValueError as exc:
        raise ValueError("SOURCE_DATE_EPOCH must be an integer") from exc
    if not 0 <= epoch <= _MAX_GZIP_MTIME:
        raise ValueError(f"SOURCE_DATE_EPOCH must be between 0 and {_MAX_GZIP_MTIME}")
    return epoch


def make_reproducible_gztar(
    base_name: str | os.PathLike[str],
    *,
    root_dir: str | os.PathLike[str] | None,
    base_dir: str,
    epoch: int,
) -> str:
    """Archive base_dir with stable ordering, ownership and timestamps."""
    if not 0 <= epoch <= _MAX_GZIP_MTIME:
        raise ValueError(f"epoch must be between 0 and {_MAX_GZIP_MTIME}")

    root = Path(root_dir or os.curdir).resolve()
    source = root / base_dir
    if not source.is_dir():
        raise FileNotFoundError(f"sdist release tree does not exist: {source}")

    archive_path = Path(f"{base_name}.tar.gz")
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    members = [source, *sorted(source.rglob("*"), key=lambda path: path.relative_to(root).as_posix())]

    with archive_path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=epoch) as compressed:
            with tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as archive:
                for path in members:
                    arcname = path.relative_to(root).as_posix()
                    info = archive.gettarinfo(str(path), arcname=arcname)
                    info.uid = 0
                    info.gid = 0
                    info.uname = ""
                    info.gname = ""
                    info.mtime = epoch
                    info.pax_headers = {}
                    if info.isfile():
                        with path.open("rb") as stream:
                            archive.addfile(info, stream)
                    else:
                        archive.addfile(info)

    return str(archive_path)
