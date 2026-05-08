import gzip
import importlib.util
import io
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FETCHER = ROOT / "scripts" / "fetch_10x_feature_barcode_demo.py"


def _load_fetcher():
    spec = importlib.util.spec_from_file_location("fetch_10x_feature_barcode_demo", FETCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _tar_header(name: str, size: int) -> bytes:
    header = bytearray(512)
    encoded = name.encode("utf-8")
    header[: min(len(encoded), 100)] = encoded[:100]
    header[100:108] = b"0000777\0"
    header[108:116] = b"0000000\0"
    header[116:124] = b"0000000\0"
    header[124:136] = f"{size:011o}\0".encode("ascii")
    header[136:148] = b"00000000000\0"
    header[156:157] = b"0"
    header[257:263] = b"ustar\0"
    header[263:265] = b"00"
    header[148:156] = b"        "
    checksum = sum(header)
    header[148:156] = f"{checksum:06o}\0 ".encode("ascii")
    return bytes(header)


def _pad(data: bytes) -> bytes:
    return data + (b"\0" * ((512 - len(data) % 512) % 512))


def _gzip_member(payload: bytes) -> bytes:
    out = io.BytesIO()
    with gzip.GzipFile(fileobj=out, mode="wb") as gz:
        gz.write(payload)
    return out.getvalue()


def test_feature_reference_parser_writes_targets_and_fixed_position_metadata(tmp_path):
    fetcher = _load_fetcher()
    feature_ref = tmp_path / "feature_ref.csv"
    targets = tmp_path / "feature_barcodes.tsv"
    feature_ref.write_text(
        "id,name,read,pattern,sequence,feature_type\n"
        "CD3,CD3,R2,^NNNNNNNNNN(BC)NNNNNNNNN,CTCATTGTAACTCCT,Antibody Capture\n"
        "CD4,CD4,R2,^NNNNNNNNNN(BC)NNNNNNNNN,TGTTCCCGCTCAACT,Antibody Capture\n",
        encoding="utf-8",
    )

    metadata = fetcher.write_targets_from_feature_ref(feature_ref, targets)

    assert targets.read_text(encoding="utf-8").splitlines() == [
        "target_id\ttarget_seq\tgene",
        "CD3\tCTCATTGTAACTCCT\tCD3",
        "CD4\tTGTTCCCGCTCAACT\tCD4",
    ]
    assert metadata == {
        "feature_count": 2,
        "read": "R2",
        "pattern": "^NNNNNNNNNN(BC)NNNNNNNNN",
        "target_start": 10,
        "target_length": 15,
        "feature_type": "Antibody Capture",
    }


def test_tar_member_lookup_resolves_gnu_longlink_names(monkeypatch):
    fetcher = _load_fetcher()
    long_name = (
        "1k_PBMCs_TotalSeq_B_3p_fastqs/antibody/"
        "1k_PBMCs_TotalSeq_B_3p_antibody_S5_L001_R2_001.fastq.gz"
    )
    payload = b"compressed-fastq"
    longlink = long_name.encode("utf-8") + b"\0"
    blob = (
        _tar_header("././@LongLink", len(longlink))
        + _pad(longlink)
        + _tar_header(long_name[:100], len(payload))
        + _pad(payload)
        + (b"\0" * 1024)
    )

    def fake_fetch_range(url, start, end, timeout=60):
        return blob[start : end + 1]

    monkeypatch.setattr(fetcher, "fetch_range", fake_fetch_range)

    member = fetcher.find_tar_member("https://example.test/archive.tar", "L001_R2_001.fastq.gz")

    assert member.name == long_name
    assert member.size == len(payload)
    assert blob[member.data_offset : member.data_offset + member.size] == payload


def test_copy_ranged_fastq_records_handles_concatenated_gzip_prefix(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    records = (
        b"@r0\nAAAAAAAAAACCCCCCCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIII\n"
        b"@r1\nAAAAAAAAAAGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIIIIIIIIIII\n"
        b"@r2\nAAAAAAAAAATTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIII\n"
    )
    compressed = _gzip_member(records[: len(records) // 2]) + _gzip_member(records[len(records) // 2 :])
    truncated = compressed[:-4]
    member = fetcher.TarMember("reads_R2.fastq.gz", 1000, len(compressed))

    monkeypatch.setattr(fetcher, "find_tar_member", lambda url, suffix: member)
    monkeypatch.setattr(fetcher, "fetch_range", lambda url, start, end, timeout=60: truncated)

    written = fetcher.copy_first_tarred_fastq_records(
        "https://example.test/archive.tar",
        "R2.fastq.gz",
        tmp_path / "reads.fastq.gz",
        2,
        prefix_bytes=len(truncated),
    )

    assert written == 2
    with gzip.open(tmp_path / "reads.fastq.gz", "rt", encoding="utf-8") as fh:
        assert fh.read().splitlines() == [
            "@r0",
            "AAAAAAAAAACCCCCCCCCCCCCCC",
            "+",
            "IIIIIIIIIIIIIIIIIIIIIIIII",
            "@r1",
            "AAAAAAAAAAGGGGGGGGGGGGGGG",
            "+",
            "IIIIIIIIIIIIIIIIIIIIIIIII",
        ]


def test_main_records_repo_relative_metadata(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    fake_root = tmp_path / "repo"
    out = fake_root / "examples" / "feature_barcode" / "data"

    monkeypatch.setattr(fetcher, "ROOT", fake_root)
    monkeypatch.setattr(
        fetcher,
        "download_feature_reference",
        lambda url, dest: dest.write_text(
            "id,name,read,pattern,sequence,feature_type\n"
            "CD3,CD3,R2,^NNNNNNNNNN(BC)NNNNNNNNN,CTCATTGTAACTCCT,Antibody Capture\n",
            encoding="utf-8",
        ),
    )
    monkeypatch.setattr(fetcher, "md5_file", lambda path: fetcher.FEATURE_REF_MD5)
    monkeypatch.setattr(fetcher, "find_tar_member", lambda url, suffix: fetcher.TarMember("source_R2.fastq.gz", 2048, 4096))

    def fake_copy(tar_url, suffix, dest, records, prefix_bytes=fetcher.FASTQ_PREFIX_BYTES):
        dest.write_bytes(
            _gzip_member(b"@r0\nAAAAAAAAAACTCATTGTAACTCCT\n+\nIIIIIIIIIIIIIIIIIIIIIIIII\n")
        )
        return 1

    monkeypatch.setattr(fetcher, "copy_first_tarred_fastq_records", fake_copy)

    result = fetcher.main(["--out", str(out), "--records", "1"])

    assert result == 0
    metadata = json.loads((out / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["evidence_ready"] is True
    assert metadata["written_records"] == 1
    assert metadata["targets"] == "examples/feature_barcode/data/feature_barcodes.tsv"
    assert metadata["local_fastq"] == "examples/feature_barcode/data/1k_PBMCs_TotalSeq_B_3p_antibody_S5_L001_R2.subsample1.fastq.gz"
    assert not Path(metadata["targets"]).is_absolute()
    assert not Path(metadata["local_fastq"]).is_absolute()
