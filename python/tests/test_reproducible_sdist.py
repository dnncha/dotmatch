from __future__ import annotations

import hashlib
import os
import tarfile
from pathlib import Path

from _reproducible_sdist import make_reproducible_gztar, parse_source_date_epoch


def _write_release_tree(root: Path, *, mtime: int) -> None:
    release = root / "dotmatch-0.1.0"
    nested = release / "nested"
    nested.mkdir(parents=True)
    (release / "a.txt").write_text("alpha\n", encoding="utf-8")
    (nested / "b.txt").write_text("beta\n", encoding="utf-8")
    for path in (release, nested, release / "a.txt", nested / "b.txt"):
        os.utime(path, (mtime, mtime))


def test_reproducible_gztar_normalizes_filesystem_metadata(tmp_path):
    epoch = 1_700_000_000
    first_root = tmp_path / "first"
    second_root = tmp_path / "second"
    _write_release_tree(first_root, mtime=epoch - 100)
    _write_release_tree(second_root, mtime=epoch + 100)

    first = Path(
        make_reproducible_gztar(
            tmp_path / "one" / "dotmatch-0.1.0",
            root_dir=first_root,
            base_dir="dotmatch-0.1.0",
            epoch=epoch,
        )
    )
    second = Path(
        make_reproducible_gztar(
            tmp_path / "two" / "dotmatch-0.1.0",
            root_dir=second_root,
            base_dir="dotmatch-0.1.0",
            epoch=epoch,
        )
    )

    assert hashlib.sha256(first.read_bytes()).digest() == hashlib.sha256(second.read_bytes()).digest()
    assert int.from_bytes(first.read_bytes()[4:8], "little") == epoch
    with tarfile.open(first, "r:gz") as archive:
        members = archive.getmembers()
    assert [member.name for member in members] == sorted(member.name for member in members)
    assert all(member.mtime == epoch for member in members)
    assert all(member.uid == 0 and member.gid == 0 for member in members)
    assert all(member.uname == "" and member.gname == "" for member in members)


def test_source_date_epoch_parser_rejects_invalid_values():
    assert parse_source_date_epoch("0") == 0
    assert parse_source_date_epoch("1700000000") == 1_700_000_000
    for value in ("not-an-integer", "-1", str(1 << 32)):
        try:
            parse_source_date_epoch(value)
        except ValueError:
            pass
        else:
            raise AssertionError(f"expected invalid SOURCE_DATE_EPOCH to fail: {value}")
