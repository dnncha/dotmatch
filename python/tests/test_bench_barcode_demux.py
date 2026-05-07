import gzip
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "scripts" / "bench_barcode_demux.py"


def _load_bench():
    spec = importlib.util.spec_from_file_location("bench_barcode_demux", BENCH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_fastq(path: Path, reads: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for name, seq in reads:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def test_hash_splitter_exact_demuxes_longest_unique_prefix(tmp_path):
    bench = _load_bench()
    reads = tmp_path / "reads.fastq.gz"
    barcodes = tmp_path / "barcodes.tsv"
    out_dir = tmp_path / "split"
    _write_fastq(reads, [
        ("r0", "ACGTAAAA"),
        ("r1", "ACAAAAAA"),
        ("r2", "TTTTAAAA"),
        ("r3", "GGGGAAAA"),
    ])
    barcodes.write_text("short\tAC\nlong\tACGT\nother\tTTTT\n", encoding="utf-8")

    stats = bench.hash_splitter_exact(reads, bench.read_barcodes(barcodes), out_dir)

    assert stats == {"assigned_reads": "3", "unmatched_reads": "1"}
    assert (out_dir / "long.fastq").read_text(encoding="utf-8").startswith("@r0\nACGTAAAA\n")
    assert (out_dir / "short.fastq").read_text(encoding="utf-8").startswith("@r1\nACAAAAAA\n")
    assert "@r3" in (out_dir / "unknown.fastq").read_text(encoding="utf-8")
