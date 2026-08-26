import gzip
import io
import importlib.util
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FETCHER = ROOT / "scripts" / "fetch_mageck_demo.py"


def _load_fetcher():
    spec = importlib.util.spec_from_file_location("fetch_mageck_demo", FETCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fastq_gz(records: int) -> bytes:
    import io

    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as fh:
        for i in range(records):
            fh.write(f"@r{i}\nACGT\n+\nIIII\n".encode("utf-8"))
    return buf.getvalue()


class FakeResponse:
    def __init__(self, payload: bytes):
        self.payload = payload
        self.headers = {"Content-Length": str(len(payload))}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def read(self, n: int = -1) -> bytes:
        if n == -1:
            out = self.payload
            self.payload = b""
            return out
        out = self.payload[:n]
        self.payload = self.payload[n:]
        return out


def test_subsample_download_rewrites_cache_with_wrong_record_count(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    dest = tmp_path / "ERR376998.fastq.gz"
    dest.write_bytes(_fastq_gz(3))

    payload = _fastq_gz(2)
    monkeypatch.setattr(fetcher.urllib.request, "urlopen", lambda url: FakeResponse(payload))

    fetcher.download_fastq_subsample("https://example.test/sample.fastq.gz", dest, 2)

    assert fetcher.count_fastq_records(dest) == 2


def test_full_download_replaces_smaller_cached_subsample(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    dest = tmp_path / "ERR376998.fastq.gz"
    dest.write_bytes(_fastq_gz(1))
    payload = _fastq_gz(4)

    monkeypatch.setattr(fetcher.urllib.request, "urlopen", lambda url: FakeResponse(payload))

    fetcher.download("https://example.test/full.fastq.gz", dest)

    assert dest.read_bytes() == payload


def _zip_payload(files: dict[str, bytes]) -> bytes:
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, "w") as archive:
        for name, content in files.items():
            archive.writestr(name, content)
    return payload.getvalue()


def test_fetch_library_extracts_only_expected_member(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    payload = _zip_payload({
        "yusa_library.csv": b"id,gRNA.sequence,Gene\n" + b"guide,ACGT,GENE\n" * 101,
        "../outside.txt": b"must not be extracted",
    })
    monkeypatch.setattr(fetcher.urllib.request, "urlopen", lambda url: FakeResponse(payload))

    fetcher.fetch_library(tmp_path / "data")

    assert (tmp_path / "data" / "yusa_library.csv").exists()
    assert not (tmp_path / "outside.txt").exists()


def test_fetch_library_rejects_archive_without_expected_member(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    payload = _zip_payload({"../yusa_library.csv": b"unexpected"})
    monkeypatch.setattr(fetcher.urllib.request, "urlopen", lambda url: FakeResponse(payload))

    try:
        fetcher.fetch_library(tmp_path / "data")
    except RuntimeError as exc:
        assert "missing yusa_library.csv" in str(exc)
    else:
        raise AssertionError("unsafe archive should be rejected")
