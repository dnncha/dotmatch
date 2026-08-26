import importlib.util
import io
import tarfile
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
FETCHER = ROOT / "scripts" / "fetch_10x_tiny_bcl.py"


def _load_fetcher():
    spec = importlib.util.spec_from_file_location("fetch_10x_tiny_bcl", FETCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _tar_payload(members: dict[str, bytes]) -> bytes:
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w:gz") as archive:
        for name, content in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))
    return payload.getvalue()


def test_extract_tiny_bcl_writes_expected_regular_files(tmp_path):
    fetcher = _load_fetcher()
    archive = tmp_path / "tiny.tar.gz"
    archive.write_bytes(_tar_payload({
        "cellranger-tiny-bcl-1.2.0/RunInfo.xml": b"<RunInfo />",
    }))

    fetcher.extract_tiny_bcl(archive, tmp_path / "data")

    assert (tmp_path / "data" / "cellranger-tiny-bcl-1.2.0" / "RunInfo.xml").read_bytes() == b"<RunInfo />"


def test_extract_tiny_bcl_rejects_path_traversal(tmp_path):
    fetcher = _load_fetcher()
    archive = tmp_path / "tiny.tar.gz"
    archive.write_bytes(_tar_payload({
        "cellranger-tiny-bcl-1.2.0/RunInfo.xml": b"must not be extracted",
        "cellranger-tiny-bcl-1.2.0/../../outside.txt": b"unsafe",
    }))

    try:
        fetcher.extract_tiny_bcl(archive, tmp_path / "data")
    except RuntimeError as exc:
        assert "unsafe archive member" in str(exc)
    else:
        raise AssertionError("unsafe archive should be rejected")
    assert not (tmp_path / "outside.txt").exists()
    assert not (tmp_path / "data" / "cellranger-tiny-bcl-1.2.0" / "RunInfo.xml").exists()


def test_extract_tiny_bcl_rejects_links(tmp_path):
    fetcher = _load_fetcher()
    archive = tmp_path / "tiny.tar.gz"
    with tarfile.open(archive, "w:gz") as tf:
        link = tarfile.TarInfo("cellranger-tiny-bcl-1.2.0/link")
        link.type = tarfile.SYMTYPE
        link.linkname = "../../outside.txt"
        tf.addfile(link)

    with pytest.raises(RuntimeError, match="unsupported archive member"):
        fetcher.extract_tiny_bcl(archive, tmp_path / "data")
