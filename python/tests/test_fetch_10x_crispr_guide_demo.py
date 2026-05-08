import gzip
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FETCHER = ROOT / "scripts" / "fetch_10x_crispr_guide_demo.py"


def _load_fetcher():
    spec = importlib.util.spec_from_file_location("fetch_10x_crispr_guide_demo", FETCHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def test_crispr_feature_reference_selection_uses_observed_fixed_window(tmp_path):
    fetcher = _load_fetcher()
    feature_ref = tmp_path / "feature_ref.csv"
    targets = tmp_path / "guides.tsv"
    reads = tmp_path / "reads.fastq.gz"
    feature_ref.write_text(
        "id,name,read,pattern,sequence,feature_type,target_gene_id,target_gene_name\n"
        "RAB1A-2,RAB1A-2_MS,R2,TTCCAGCTTAGCTCTTAAAC(BC),TATTTCCTGGTTCGCCGGC,CRISPR Guide Capture,ENSG00000138069,RAB1A\n"
        "NT1,NT1,R2,TTCCAGCTTAGCTCTTAAAC(BC),GCCCGCATCGTCAGCACGTT,CRISPR Guide Capture,Non-Targeting,Non-Targeting\n"
        "ADT,ADT,R2,^NNNNNNNNNN(BC)NNNNNNNNN,ACGTACGTACGTACG,Antibody Capture,Ignore,Ignore\n",
        encoding="utf-8",
    )
    _write_fastq(
        reads,
        [
            ("r0", "CAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCTTAGCTCTTAAACTATTTCCTGGTTCGCCGGCC"),
            ("r1", "CAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCTTAGCTCTTAAACTATTTCCTGGTTCGCCGGCA"),
            ("r2", "CAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCTTAGCTCTTAAACGCCCGCATCGTCAGCACGTT"),
        ],
    )

    metadata = fetcher.write_observed_targets_from_feature_ref(feature_ref, reads, targets)

    assert targets.read_text(encoding="utf-8").splitlines() == [
        "target_id\ttarget_seq\tgene",
        "RAB1A-2\tTATTTCCTGGTTCGCCGGC\tRAB1A",
    ]
    assert metadata["feature_count"] == 1
    assert metadata["target_start"] == 63
    assert metadata["target_length"] == 19
    assert metadata["selected_pattern"] == "TTCCAGCTTAGCTCTTAAAC(BC)"
    assert metadata["exact_assigned_in_subsample"] == 2
    assert metadata["total_reads_scanned"] == 3
    assert metadata["dropped_feature_count"] == 1


def test_main_records_repo_relative_crispr_metadata(tmp_path, monkeypatch):
    fetcher = _load_fetcher()
    fake_root = tmp_path / "repo"
    out = fake_root / "examples" / "perturb_seq" / "data"

    monkeypatch.setattr(fetcher, "ROOT", fake_root)
    monkeypatch.setattr(
        fetcher,
        "download_feature_reference",
        lambda url, dest: dest.write_text(
            "id,name,read,pattern,sequence,feature_type,target_gene_id,target_gene_name\n"
            "RAB1A-2,RAB1A-2_MS,R2,TTCCAGCTTAGCTCTTAAAC(BC),TATTTCCTGGTTCGCCGGC,CRISPR Guide Capture,ENSG00000138069,RAB1A\n",
            encoding="utf-8",
        ),
    )
    monkeypatch.setattr(fetcher, "md5_file", lambda path: fetcher.FEATURE_REF_MD5)

    def fake_write_fastq(tar_url, suffix, dest, records, prefix_bytes=fetcher.FASTQ_PREFIX_BYTES, tar_bytes=None):
        _write_fastq(
            dest,
            [("r0", "CAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCTTAGCTCTTAAACTATTTCCTGGTTCGCCGGCC")],
        )
        return 1

    monkeypatch.setattr(fetcher, "copy_first_tarred_fastq_records", fake_write_fastq)

    result = fetcher.main(["--out", str(out), "--records", "1"])

    assert result == 0
    metadata = json.loads((out / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["evidence_ready"] is True
    assert metadata["targets"] == "examples/perturb_seq/data/crispr_guides.tsv"
    assert metadata["local_fastq"] == "examples/perturb_seq/data/1k_CRISPR_5p_gemx_crispr_S1_L001_R2.subsample1.fastq.gz"
    assert metadata["target_start"] == 63
    assert metadata["target_length"] == 19
    assert not Path(metadata["targets"]).is_absolute()
    assert not Path(metadata["local_fastq"]).is_absolute()
