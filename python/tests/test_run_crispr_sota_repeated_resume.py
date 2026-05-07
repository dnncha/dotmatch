import argparse
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "scripts" / "run_crispr_sota_repeated.py"


def _load_runner():
    spec = importlib.util.spec_from_file_location("run_crispr_sota_repeated", RUNNER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_resume_completeness_requires_all_requested_successful_tools():
    runner = _load_runner()
    args = argparse.Namespace(run_mageck=True, run_guide_counter=True, skip_levenshtein=False)
    rows = [
        {"dataset_id": "sanson_brunello", "requested_records_per_sample": "100000", "repeat": "1", "run_level": "subsample", "tool": "dotmatch_exact_k0", "exit_code": "0"},
        {"dataset_id": "sanson_brunello", "requested_records_per_sample": "100000", "repeat": "1", "run_level": "subsample", "tool": "dotmatch_hamming_k1", "exit_code": "0"},
        {"dataset_id": "sanson_brunello", "requested_records_per_sample": "100000", "repeat": "1", "run_level": "subsample", "tool": "dotmatch_levenshtein_k1", "exit_code": "0"},
        {"dataset_id": "sanson_brunello", "requested_records_per_sample": "100000", "repeat": "1", "run_level": "subsample", "tool": "mageck_count_exact", "exit_code": "0"},
    ]

    assert not runner.run_complete(rows, "sanson_brunello", "100000", 1, "subsample", runner.expected_tools(args))

    rows.append({"dataset_id": "sanson_brunello", "requested_records_per_sample": "100000", "repeat": "1", "run_level": "subsample", "tool": "guide_counter_one_mismatch", "exit_code": "0"})

    assert runner.run_complete(rows, "sanson_brunello", "100000", 1, "subsample", runner.expected_tools(args))


def test_resume_ignores_failed_existing_rows():
    runner = _load_runner()
    args = argparse.Namespace(run_mageck=False, run_guide_counter=True, skip_levenshtein=True)
    rows = [
        {"dataset_id": "mageck_yusa", "requested_records_per_sample": "100000", "repeat": "2", "run_level": "subsample", "tool": "dotmatch_exact_k0", "exit_code": "0"},
        {"dataset_id": "mageck_yusa", "requested_records_per_sample": "100000", "repeat": "2", "run_level": "subsample", "tool": "dotmatch_hamming_k1", "exit_code": "0"},
        {"dataset_id": "mageck_yusa", "requested_records_per_sample": "100000", "repeat": "2", "run_level": "subsample", "tool": "guide_counter_one_mismatch", "exit_code": "127"},
    ]

    assert not runner.run_complete(rows, "mageck_yusa", "100000", 2, "subsample", runner.expected_tools(args))


def test_full_only_mode_skips_subsample_phase_and_runs_full_phase():
    runner = _load_runner()

    args = argparse.Namespace(full=False, full_only=True)
    assert not runner.should_run_subsamples(args)
    assert runner.should_run_full(args)

    args = argparse.Namespace(full=False, full_only=False)
    assert runner.should_run_subsamples(args)
    assert not runner.should_run_full(args)


def test_dataset_filter_selects_named_runner_and_rejects_unknown_dataset():
    runner = _load_runner()

    args = argparse.Namespace(datasets="mageck_yusa")
    selected = runner.selected_dataset_runners(args)
    assert [name for name, _ in selected] == ["mageck_yusa"]

    args = argparse.Namespace(datasets="mageck_yusa,unknown")
    try:
        runner.selected_dataset_runners(args)
    except SystemExit as exc:
        assert "unknown dataset" in str(exc)
    else:
        raise AssertionError("expected unknown dataset to fail")


def test_full_sample_resume_is_scoped_to_selected_sanson_sample():
    runner = _load_runner()
    args = argparse.Namespace(run_mageck=False, run_guide_counter=False, skip_levenshtein=True,
                              sanson_samples="RepB")
    rows = [
        {"dataset_id": "sanson_brunello", "requested_records_per_sample": "full", "repeat": "1",
         "run_level": "full_sample", "sample_id": "RepA", "tool": "dotmatch_exact_k0", "exit_code": "0"},
        {"dataset_id": "sanson_brunello", "requested_records_per_sample": "full", "repeat": "1",
         "run_level": "full_sample", "sample_id": "RepA", "tool": "dotmatch_hamming_k1", "exit_code": "0"},
    ]

    assert runner.full_sample_scope("sanson_brunello", args) == "RepB"
    assert runner.full_run_level("sanson_brunello", args) == "full_sample"
    assert not runner.run_complete(rows, "sanson_brunello", "full", 1, "full_sample",
                                   runner.expected_tools(args), sample_id="RepB")
    assert runner.run_complete(rows, "sanson_brunello", "full", 1, "full_sample",
                               runner.expected_tools(args), sample_id="RepA")


def test_fetch_sanson_passes_selected_samples(monkeypatch):
    runner = _load_runner()
    captured = []

    def fake_run_cmd(cmd, env=None):
        captured.append(cmd)

    monkeypatch.setattr(runner, "run_cmd", fake_run_cmd)

    manifest = runner.fetch_sanson(None, argparse.Namespace(sanson_samples="RepC"))

    assert manifest == runner.SANSON_DATA / "manifest.json"
    assert "--samples" in captured[0]
    assert captured[0][captured[0].index("--samples") + 1] == "RepC"


def test_scoped_manifest_path_is_sample_specific_for_full_sample_runs():
    runner = _load_runner()
    manifest = Path("/tmp/sanson/manifest.json")

    assert runner.scoped_manifest_path(manifest, "RepB") == Path("/tmp/sanson/manifest.RepB.json")
    assert runner.scoped_manifest_path(manifest, "") == manifest
