import gzip
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FETCHER = ROOT / "scripts" / "fetch_nfcore_amplicon_panel_demo.py"


def _load_fetcher():
    spec = importlib.util.spec_from_file_location("fetch_nfcore_amplicon_panel_demo", FETCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def test_primer_fasta_selection_uses_observed_full_length_group(tmp_path):
    fetcher = _load_fetcher()
    primer_fasta = tmp_path / "primers.fasta"
    targets = tmp_path / "targets.tsv"
    reads = tmp_path / "reads.fastq.gz"
    primer_fasta.write_text(
        ">primer_1_LEFT\nACGT\n"
        ">primer_2_LEFT\nTTTT\n"
        ">primer_3_LEFT\nGGGGAA\n",
        encoding="utf-8",
    )
    _write_fastq(
        reads,
        [
            ("r0", "ACGTAAAA"),
            ("r1", "ACGTCCCC"),
            ("r2", "TTTTCCCC"),
            ("r3", "GGGGAATTTT"),
        ],
    )

    metadata = fetcher.write_observed_targets_from_primer_fasta(primer_fasta, reads, targets)

    assert targets.read_text(encoding="utf-8").splitlines() == [
        "target_id\ttarget_seq\tgene",
        "primer_1_LEFT\tACGT\tprimer_1",
        "primer_2_LEFT\tTTTT\tprimer_2",
    ]
    assert metadata["selected_target_length"] == 4
    assert metadata["feature_count"] == 2
    assert metadata["exact_assigned_in_subsample"] == 3
    assert metadata["total_reads_scanned"] == 4
    assert metadata["dropped_feature_count"] == 1


def test_main_records_repo_relative_amplicon_metadata(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    fake_root = tmp_path / "repo"
    out = fake_root / "examples" / "amplicon_panel" / "data"

    monkeypatch.setattr(fetcher, "ROOT", fake_root)

    def fake_download(url, dest):
        if str(dest).endswith(".fasta"):
            dest.write_text(">primer_1_LEFT\nACGT\n", encoding="utf-8")
        else:
            _write_fastq(dest, [("r0", "ACGTAAAA")])

    monkeypatch.setattr(fetcher, "download_url", fake_download)
    def fake_copy(src, dest, records):
        _write_fastq(dest, [("r0", "ACGTAAAA")])
        return 1

    monkeypatch.setattr(fetcher, "copy_first_fastq_records", fake_copy)

    result = fetcher.main(["--out", str(out), "--records", "1"])

    assert result == 0
    metadata = json.loads((out / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["evidence_ready"] is True
    assert metadata["targets"] == "examples/amplicon_panel/data/artic_v3_primers_len4.tsv"
    assert metadata["local_fastq"] == "examples/amplicon_panel/data/nfcore_viralrecon_sample1_R1.subsample1.fastq.gz"
    assert metadata["target_start"] == 0
    assert metadata["target_length"] == 4
    assert not Path(metadata["targets"]).is_absolute()
    assert not Path(metadata["local_fastq"]).is_absolute()
