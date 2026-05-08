import importlib.util
import gzip
import http.client
import io
import json
import sys
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FETCHER = ROOT / "scripts" / "fetch_srp009896_barcode_demo.py"


def _load_fetcher():
    spec = importlib.util.spec_from_file_location("fetch_srp009896_barcode_demo", FETCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _zip_bytes() -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("BarcodesPerSample.csv", "sample_a\tACGT\tSRR391079\nsample_b\tTGCA\tSRR391079\n")
        zf.writestr("large.fastq.gz", b"ignored")
    return buf.getvalue()


def test_extracts_first_zip_member_from_prefix_without_central_directory():
    fetcher = _load_fetcher()
    payload = _zip_bytes()
    member, content = fetcher.extract_first_zip_member(payload[:180])

    assert member == "BarcodesPerSample.csv"
    assert content.decode("utf-8").splitlines() == [
        "sample_a\tACGT\tSRR391079",
        "sample_b\tTGCA\tSRR391079",
    ]


def test_fetches_public_example_barcodes_through_google_drive_warning(monkeypatch):
    fetcher = _load_fetcher()
    warning = (
        '<form id="download-form" action="https://drive.usercontent.google.com/download" method="get">'
        '<input type="hidden" name="id" value="FILEID">'
        '<input type="hidden" name="export" value="download">'
        '<input type="hidden" name="confirm" value="t">'
        '<input type="hidden" name="uuid" value="UUID">'
        "</form>"
    ).encode("utf-8")
    calls = []

    def fake_fetch_range(url, start, end, timeout=60):
        calls.append(url)
        if "confirm=t" not in url:
            return warning
        return _zip_bytes()[:180]

    monkeypatch.setattr(fetcher, "fetch_range", fake_fetch_range)

    member, content = fetcher.fetch_first_zip_member(
        "https://drive.google.com/file/d/FILEID/view?usp=sharing"
    )

    assert member == "BarcodesPerSample.csv"
    assert b"sample_a\tACGT" in content
    assert any("confirm=t" in url and "uuid=UUID" in url for url in calls)


def test_fetch_range_returns_nonempty_incomplete_read_partial(monkeypatch):
    fetcher = _load_fetcher()

    class FakeResponse:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def read(self):
            raise http.client.IncompleteRead(b"ACGT", 100)

    monkeypatch.setattr(fetcher, "urlopen_with_retries", lambda request, timeout=60: FakeResponse())

    assert fetcher.fetch_range("https://example.test/r.fastq.gz", 0, 100) == b"ACGT"


def test_main_can_install_public_example_barcodes_without_fastq_download(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    fake_root = tmp_path / "repo"
    out_dir = fake_root / "examples" / "barcode_demux" / "data"
    monkeypatch.setattr(fetcher, "ROOT", fake_root)

    def fake_ena_metadata(accession):
        return {
            "run_accession": accession,
            "fastq_ftp": f"ftp.sra.ebi.ac.uk/{accession}.fastq.gz",
            "fastq_bytes": "1000",
            "fastq_md5": "remote-md5",
            "read_count": "100",
            "base_count": "8000",
            "sample_alias": "sample",
            "experiment_title": "experiment",
            "study_accession": "SRP009896",
        }

    monkeypatch.setattr(fetcher, "ena_metadata", fake_ena_metadata)
    monkeypatch.setattr(
        fetcher,
        "fetch_first_zip_member",
        lambda url, member_name="BarcodesPerSample.csv": (
            "BarcodesPerSample.csv",
            b"sample_a\tACGT\tSRR391079\nsample_b\tACGT\tSRR391080\n",
        ),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "fetch_srp009896",
            "--out",
            str(out_dir),
            "--metadata-only",
            "--use-public-example-barcodes",
            "--require-barcodes",
        ],
    )

    fetcher.main()

    assert (out_dir / "barcodes.tsv").read_text(encoding="utf-8") == "sample_a\tACGT\tSRR391079\n"
    metadata = json.loads((out_dir / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["evidence_ready"] is True
    assert [run["accession"] for run in metadata["runs"]] == ["SRR391079"]
    assert metadata["barcodes"] == "examples/barcode_demux/data/barcodes.tsv"
    assert metadata["runs"][0]["local_fastq"] == "examples/barcode_demux/data/SRR391079.subsample100000.fastq.gz"
    assert not Path(metadata["barcodes"]).is_absolute()
    assert not Path(metadata["runs"][0]["local_fastq"]).is_absolute()


def test_explicit_default_accession_does_not_duplicate_metadata_run(tmp_path, monkeypatch):
    fetcher = _load_fetcher()

    def fake_ena_metadata(accession):
        return {
            "run_accession": accession,
            "fastq_ftp": f"ftp.sra.ebi.ac.uk/{accession}.fastq.gz",
            "fastq_bytes": "1000",
            "fastq_md5": "remote-md5",
            "read_count": "100",
            "base_count": "8000",
            "sample_alias": "sample",
            "experiment_title": "experiment",
            "study_accession": "SRP009896",
        }

    monkeypatch.setattr(fetcher, "ena_metadata", fake_ena_metadata)
    monkeypatch.setattr(
        fetcher,
        "fetch_first_zip_member",
        lambda url, member_name="BarcodesPerSample.csv": (
            "BarcodesPerSample.csv",
            b"sample_a\tACGT\tSRR391079\nsample_b\tACGT\tSRR391080\n",
        ),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "fetch_srp009896",
            "--out",
            str(tmp_path),
            "--metadata-only",
            "--accession",
            "SRR391079",
            "--use-public-example-barcodes",
        ],
    )

    fetcher.main()

    metadata = json.loads((tmp_path / "metadata.json").read_text(encoding="utf-8"))
    assert [run["accession"] for run in metadata["runs"]] == ["SRR391079"]
    assert metadata["barcode_filter_applied"] is True
    assert metadata["barcode_rows_before_filter"] == 2
    assert metadata["barcode_rows_after_filter"] == 1
    assert (tmp_path / "barcodes.tsv").read_text(encoding="utf-8") == "sample_a\tACGT\tSRR391079\n"


def test_subsample_fastq_retries_transient_stream_open_failure(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    payload = io.BytesIO()
    with gzip.GzipFile(fileobj=payload, mode="wb") as gz:
        gz.write(
            b"@r0\nACGT\n+\nIIII\n"
            b"@r1\nTGCA\n+\nJJJJ\n"
        )
    attempts = {"count": 0}

    def fake_fetch_range(url, start, end, timeout=60):
        attempts["count"] += 1
        if attempts["count"] == 1:
            raise http.client.RemoteDisconnected("closed")
        return payload.getvalue()

    monkeypatch.setattr(fetcher, "fetch_range", fake_fetch_range)
    monkeypatch.setattr(fetcher.time, "sleep", lambda seconds: None)

    written = fetcher.copy_first_fastq_records("https://example.test/run.fastq.gz", tmp_path / "sub.fastq.gz", 2)

    assert written == 2
    assert attempts["count"] == 2
    with gzip.open(tmp_path / "sub.fastq.gz", "rt") as fh:
        assert fh.read().splitlines() == ["@r0", "ACGT", "+", "IIII", "@r1", "TGCA", "+", "JJJJ"]


def test_subsample_fastq_uses_gzip_prefix_without_remote_trailer(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    payload = io.BytesIO()
    with gzip.GzipFile(fileobj=payload, mode="wb") as gz:
        gz.write(
            b"@r0\nACGT\n+\nIIII\n"
            b"@r1\nTGCA\n+\nJJJJ\n"
            b"@r2\nGGGG\n+\nHHHH\n"
        )
    truncated = payload.getvalue()[:-8]
    calls = []

    def fake_fetch_range(url, start, end, timeout=60):
        calls.append((start, end))
        return truncated if start == 0 else b""

    monkeypatch.setattr(fetcher, "fetch_range", fake_fetch_range)

    written = fetcher.copy_first_fastq_records("https://example.test/run.fastq.gz", tmp_path / "sub.fastq.gz", 2)

    assert written == 2
    assert calls == [(0, fetcher.FASTQ_RANGE_CHUNK_BYTES - 1)]
    with gzip.open(tmp_path / "sub.fastq.gz", "rt") as fh:
        assert fh.read().splitlines() == ["@r0", "ACGT", "+", "IIII", "@r1", "TGCA", "+", "JJJJ"]
