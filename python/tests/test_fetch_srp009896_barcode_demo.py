import importlib.util
import io
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


def test_main_can_install_public_example_barcodes_without_fastq_download(tmp_path, monkeypatch):
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
        lambda url, member_name="BarcodesPerSample.csv": ("BarcodesPerSample.csv", b"sample_a\tACGT\tSRR391079\n"),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "fetch_srp009896",
            "--out",
            str(tmp_path),
            "--metadata-only",
            "--use-public-example-barcodes",
            "--require-barcodes",
        ],
    )

    fetcher.main()

    assert (tmp_path / "barcodes.tsv").read_text(encoding="utf-8") == "sample_a\tACGT\tSRR391079\n"
    assert '"claim_grade_ready": true' in (tmp_path / "metadata.json").read_text(encoding="utf-8")
