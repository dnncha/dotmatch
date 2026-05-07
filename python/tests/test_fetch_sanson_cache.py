import gzip
import hashlib
import importlib.util
import io
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FETCHER = ROOT / "scripts" / "fetch_sanson_brunello_demo.py"


def _load_fetcher():
    spec = importlib.util.spec_from_file_location("fetch_sanson_brunello_demo", FETCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_fastq_gz(path: Path, records: int) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for i in range(records):
            fh.write(f"@read{i}\n")
            fh.write("ACGT\n")
            fh.write("+\n")
            fh.write("IIII\n")


def _remote_fastq_bytes(accession: str) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
        gz.write(f"@{accession}\nACGT\n+\nIIII\n".encode("ascii"))
    return buf.getvalue()


def test_sanson_fetcher_reuses_existing_subsample_fastqs(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    out = tmp_path / "sanson"
    out.mkdir()
    (out / "broadgpp-brunello-library-corrected.txt").write_text(
        "id,sequence,gene\nsg1,ACGT,GENE\n",
        encoding="utf-8",
    )
    for sample in fetcher.SAMPLES:
        _write_fastq_gz(out / f"{sample}.subsample2.fastq.gz", 2)

    def fake_ena_metadata(accession):
        return {
            "run_accession": accession,
            "fastq_ftp": f"ftp.sra.ebi.ac.uk/{accession}.fastq.gz",
            "read_count": "100",
            "base_count": "400",
            "fastq_bytes": "1000",
            "fastq_md5": "remote-md5",
            "sample_alias": "sample",
            "experiment_title": "experiment",
        }

    def fail_remote_copy(*args, **kwargs):
        raise AssertionError("cached subsamples should avoid remote FASTQ download")

    monkeypatch.setattr(fetcher, "ena_metadata", fake_ena_metadata)
    monkeypatch.setattr(fetcher, "copy_first_records_from_urls", fail_remote_copy)
    monkeypatch.setattr(sys, "argv", ["fetch_sanson", "--out", str(out), "--subsample", "2"])

    fetcher.main()

    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert [sample["written_records"] for sample in manifest["samples"]] == [2, 2, 2, 2]


def test_sanson_full_fetch_streams_without_source_fastq_cache(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    out = tmp_path / "sanson"
    out.mkdir()
    (out / "broadgpp-brunello-library-corrected.txt").write_text(
        "id,sequence,gene\nsg1,ACGT,GENE\n",
        encoding="utf-8",
    )

    def fake_ena_metadata(accession):
        return {
            "run_accession": accession,
            "fastq_ftp": f"ftp.sra.ebi.ac.uk/{accession}.fastq.gz",
            "read_count": "1",
            "base_count": "4",
            "fastq_bytes": str(len(_remote_fastq_bytes(accession))),
            "fastq_md5": hashlib.md5(_remote_fastq_bytes(accession)).hexdigest(),
            "sample_alias": "sample",
            "experiment_title": "experiment",
        }

    class Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            self.close()
            return False

    def fake_urlopen(url, timeout=0):
        accession = Path(url).name.split(".", 1)[0]
        return Response(_remote_fastq_bytes(accession))

    monkeypatch.setattr(fetcher, "ena_metadata", fake_ena_metadata)
    monkeypatch.setattr(fetcher.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(sys, "argv", ["fetch_sanson", "--out", str(out), "--subsample", "0"])

    fetcher.main()

    assert sorted(path.name for path in out.glob("*.source.fastq.gz")) == []
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert [sample["written_records"] for sample in manifest["samples"]] == [1, 2, 2, 2]
    assert fetcher.count_fastq_gz(out / "RepA.fastq.gz") == 2


def test_sanson_full_fetch_replaces_incomplete_cached_fastq(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    out = tmp_path / "sanson"
    out.mkdir()
    (out / "broadgpp-brunello-library-corrected.txt").write_text(
        "id,sequence,gene\nsg1,ACGT,GENE\n",
        encoding="utf-8",
    )
    _write_fastq_gz(out / "RepC.fastq.gz", 1)

    def fake_ena_metadata(accession):
        return {
            "run_accession": accession,
            "fastq_ftp": f"ftp.sra.ebi.ac.uk/{accession}.fastq.gz",
            "read_count": "1",
            "base_count": "4",
            "fastq_bytes": str(len(_remote_fastq_bytes(accession))),
            "fastq_md5": hashlib.md5(_remote_fastq_bytes(accession)).hexdigest(),
            "sample_alias": "sample",
            "experiment_title": "experiment",
        }

    class Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            self.close()
            return False

    def fake_urlopen(url, timeout=0):
        accession = Path(url).name.split(".", 1)[0]
        return Response(_remote_fastq_bytes(accession))

    monkeypatch.setattr(fetcher, "ena_metadata", fake_ena_metadata)
    monkeypatch.setattr(fetcher.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(sys, "argv", ["fetch_sanson", "--out", str(out), "--subsample", "0", "--samples", "RepC"])

    fetcher.main()

    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["samples"][0]["written_records"] == 2
    assert fetcher.count_fastq_gz(out / "RepC.fastq.gz") == 2


def test_sanson_fetcher_can_limit_to_selected_samples(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    out = tmp_path / "sanson"
    out.mkdir()
    (out / "broadgpp-brunello-library-corrected.txt").write_text(
        "id,sequence,gene\nsg1,ACGT,GENE\n",
        encoding="utf-8",
    )

    def fake_ena_metadata(accession):
        return {
            "run_accession": accession,
            "fastq_ftp": f"ftp.sra.ebi.ac.uk/{accession}.fastq.gz",
            "read_count": "1",
            "base_count": "4",
            "fastq_bytes": "40",
            "fastq_md5": "remote-md5",
            "sample_alias": "sample",
            "experiment_title": "experiment",
        }

    def fake_copy_first_records(urls, out_path, records):
        _write_fastq_gz(out_path, records)
        return records

    monkeypatch.setattr(fetcher, "ena_metadata", fake_ena_metadata)
    monkeypatch.setattr(fetcher, "copy_first_records_from_urls", fake_copy_first_records)
    monkeypatch.setattr(sys, "argv", ["fetch_sanson", "--out", str(out), "--subsample", "2", "--samples", "RepB"])

    fetcher.main()

    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert [sample["sample_id"] for sample in manifest["samples"]] == ["RepB"]
    assert (out / "RepB.subsample2.fastq.gz").exists()
    assert not (out / "RepA.subsample2.fastq.gz").exists()
