import gzip
import importlib.util
import io
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FETCHER = ROOT / "scripts" / "fetch_public_oligo_adapter_demo.py"


def _load_fetcher():
    spec = importlib.util.spec_from_file_location("fetch_public_oligo_adapter_demo", FETCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _gzip_member(payload: bytes) -> bytes:
    out = io.BytesIO()
    with gzip.GzipFile(fileobj=out, mode="wb") as gz:
        gz.write(payload)
    return out.getvalue()


def test_adapter_fasta_parser_deduplicates_prefix_targets(tmp_path):
    fetcher = _load_fetcher()
    fasta = tmp_path / "truseq.fa"
    targets = tmp_path / "adapter_oligos.tsv"
    fasta.write_text(
        ">TruSeq_I7_full_length\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n"
        ">TruSeq_I7_Recommended\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n"
        ">TruSeq_I7_threeprime\nATCTCGTATGCCGTCTTCTGCTTG\n",
        encoding="utf-8",
    )

    metadata = fetcher.write_targets_from_adapter_fasta(fasta, targets, target_length=20)

    assert targets.read_text(encoding="utf-8").splitlines() == [
        "target_id\ttarget_seq\tgene",
        "TruSeq_I7_full_length\tAGATCGGAAGAGCACACGTC\tTruSeq_I7_full_length",
        "TruSeq_I7_threeprime\tATCTCGTATGCCGTCTTCTG\tTruSeq_I7_threeprime",
    ]
    assert metadata["adapter_count"] == 3
    assert metadata["target_count"] == 2
    assert metadata["target_length"] == 20
    assert metadata["duplicate_prefixes"] == {
        "AGATCGGAAGAGCACACGTC": ["TruSeq_I7_full_length", "TruSeq_I7_Recommended"]
    }


def test_copy_first_public_fastq_records_writes_valid_gzip(tmp_path):
    fetcher = _load_fetcher()
    source = _gzip_member(
        b"@r0\n" + b"A" * 260 + b"\n+\n" + b"I" * 260 + b"\n"
        b"@r1\n" + b"C" * 260 + b"\n+\n" + b"I" * 260 + b"\n"
        b"@r2\n" + b"G" * 260 + b"\n+\n" + b"I" * 260 + b"\n"
    )

    monkey_data = {"payload": source}

    def fake_fetch_bytes(url, timeout=60):
        return monkey_data["payload"]

    fetcher.fetch_bytes = fake_fetch_bytes

    written = fetcher.copy_first_fastq_records(
        "https://example.test/reads.fastq.gz",
        tmp_path / "reads.subsample2.fastq.gz",
        records=2,
    )

    assert written == 2
    with gzip.open(tmp_path / "reads.subsample2.fastq.gz", "rt", encoding="utf-8") as fh:
        assert fh.read().splitlines() == [
            "@r0",
            "A" * 260,
            "+",
            "I" * 260,
            "@r1",
            "C" * 260,
            "+",
            "I" * 260,
        ]


def test_main_records_repo_relative_public_adapter_metadata(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    fake_root = tmp_path / "repo"
    out = fake_root / "examples" / "oligo_adapter" / "data"

    monkeypatch.setattr(fetcher, "ROOT", fake_root)
    monkeypatch.setattr(
        fetcher,
        "download_file",
        lambda url, dest: dest.write_text(
            ">TruSeq_I7_full_length\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n",
            encoding="utf-8",
        )
        if dest.suffix == ".fa"
        else dest.write_bytes(
            _gzip_member(b"@r0\n" + b"A" * 260 + b"\n+\n" + b"I" * 260 + b"\n")
        ),
    )
    monkeypatch.setattr(fetcher, "md5_file", lambda path: "md5")
    monkeypatch.setattr(fetcher, "sha256_file", lambda path: "sha256")

    result = fetcher.main(["--out", str(out), "--records", "1"])

    assert result == 0
    metadata = json.loads((out / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["evidence_ready"] is True
    assert metadata["dataset_id"] == "fast_adapter_trimming_truseq_r1"
    assert metadata["target_start"] == 229
    assert metadata["target_length"] == 20
    assert metadata["written_records"] == 1
    assert metadata["targets"] == "examples/oligo_adapter/data/adapter_oligos.tsv"
    assert metadata["local_fastq"] == "examples/oligo_adapter/data/788707_20180313_S_R1.small.subsample1.fastq.gz"
    assert not Path(metadata["targets"]).is_absolute()
    assert not Path(metadata["local_fastq"]).is_absolute()
