import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "scripts" / "run_crispr_dataset_benchmark.py"


def _load_bench():
    spec = importlib.util.spec_from_file_location("run_crispr_dataset_benchmark", BENCH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_single_sample_scope_sets_sample_id_for_full_chunk_rows():
    bench = _load_bench()

    assert bench.sample_scope(["RepB"]) == "RepB"
    assert bench.sample_scope(["plasmid", "RepA"]) == ""
