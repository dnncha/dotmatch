import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from dotmatch.feature_matrix import build_feature_matrix


ROOT = Path(__file__).resolve().parents[2]
LEGACY_ENV = {**os.environ, "DOTMATCH_PYTHON_NO_DELEGATE": "1"}


def _write_inputs(tmp_path: Path) -> tuple[Path, Path]:
    targets = tmp_path / "features.tsv"
    targets.write_text(
        "target_id\ttarget_seq\n"
        "feature_a\tAAAA\n"
        "feature_b\tAACC\n"
        "feature_c\tTTTT\n",
        encoding="utf-8",
    )
    observations = tmp_path / "observations.tsv"
    observations.write_text(
        "observation_id\tcell_barcode\tfeature_seq\n"
        "obs_001\tcell_z\tAAAA\n"
        "obs_002\tcell_z\tAAAT\n"
        "obs_003\tcell_a\tTTTT\n"
        "obs_004\tcell_a\tAAAC\n"
        "obs_005\tcell_a\tCCCC\n"
        "obs_006\tcell_a\t\n",
        encoding="utf-8",
    )
    return observations, targets


def test_build_feature_matrix_writes_deterministic_sparse_artifacts(tmp_path: Path) -> None:
    observations, targets = _write_inputs(tmp_path)
    result = build_feature_matrix(
        observations,
        targets,
        tmp_path / "matrix",
        cell_column="cell_barcode",
        sequence_column="feature_seq",
        id_column="observation_id",
        k=1,
        metric="hamming",
    )

    assert result.summary["matrix_orientation"] == "cells_by_features"
    assert result.summary["total_observations"] == 6
    assert result.summary["valid_observations"] == 5
    assert result.summary["assigned_unique"] == 3
    assert result.summary["assigned_exact"] == 2
    assert result.summary["assigned_corrected"] == 1
    assert result.summary["ambiguous"] == 1
    assert result.summary["unmatched"] == 1
    assert result.summary["invalid"] == 1
    assert result.summary["cells"] == 2
    assert result.summary["features"] == 3
    assert result.summary["nonzero_entries"] == 2
    assert result.summary["scope"]["umi_deduplication"] == "not_performed"

    output = result.output_dir
    assert (output / "barcodes.tsv").read_text(encoding="utf-8") == "cell_barcode\ncell_a\ncell_z\n"
    assert (output / "features.tsv").read_text(encoding="utf-8") == (
        "target_id\ttarget_seq\nfeature_a\tAAAA\nfeature_b\tAACC\nfeature_c\tTTTT\n"
    )
    assert (output / "cell_feature_counts.tsv").read_text(encoding="utf-8") == (
        "cell_barcode\ttarget_id\tcount\ncell_a\tfeature_c\t1\ncell_z\tfeature_a\t2\n"
    )
    assert (output / "matrix.mtx").read_text(encoding="utf-8") == (
        "%%MatrixMarket matrix coordinate integer general\n"
        "% DotMatch cell-by-feature unique-assignment counts\n"
        "2 3 2\n"
        "1 3 1\n"
        "2 1 2\n"
    )

    assignments = (output / "assignments.tsv").read_text(encoding="utf-8")
    assert "obs_004\tcell_a\tAAAC\t\t\t1\tambiguous\t2\t-1" in assignments
    assert "obs_005\tcell_a\tCCCC\t\t\t-1\tnone\t0\t-1" in assignments
    assert "obs_006\tcell_a\t\t\t\t-1\tinvalid\t0\t-1" in assignments

    qc = (output / "cell_qc.tsv").read_text(encoding="utf-8").splitlines()
    assert qc[1] == "cell_a\t4\t3\t1\t1\t1\t1\t1\t0.3333333333333333"
    assert qc[2] == "cell_z\t2\t2\t2\t0\t0\t0\t1\t1.0"

    written_summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert written_summary == result.summary


def test_build_feature_matrix_requires_explicit_existing_columns(tmp_path: Path) -> None:
    observations, targets = _write_inputs(tmp_path)

    with pytest.raises(ValueError, match="missing required column 'wrong_cell'"):
        build_feature_matrix(
            observations,
            targets,
            tmp_path / "matrix",
            cell_column="wrong_cell",
            sequence_column="feature_seq",
        )

    assert not (tmp_path / "matrix").exists()


def test_feature_matrix_cli_writes_full_artifact_set(tmp_path: Path) -> None:
    observations, targets = _write_inputs(tmp_path)
    output = tmp_path / "matrix"
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "dotmatch.cli",
            "feature",
            "matrix",
            "--observations",
            str(observations),
            "--targets",
            str(targets),
            "--id-column",
            "observation_id",
            "--cell-column",
            "cell_barcode",
            "--sequence-column",
            "feature_seq",
            "--k",
            "1",
            "--metric",
            "hamming",
            "--out-dir",
            str(output),
        ],
        check=False,
        env=LEGACY_ENV,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout)["assigned_unique"] == 3
    assert {path.name for path in output.iterdir()} == {
        "assignments.tsv",
        "barcodes.tsv",
        "cell_feature_counts.tsv",
        "cell_qc.tsv",
        "features.tsv",
        "matrix.mtx",
        "summary.json",
    }


def test_assignments_to_anndata_accepts_text_status_when_optional_deps_are_available() -> None:
    pd = pytest.importorskip("pandas")
    pytest.importorskip("anndata")
    import dotmatch

    assignments = pd.DataFrame(
        {
            "cell_barcode": ["cell_1", "cell_1", "cell_2"],
            "target_id": ["feature_a", "feature_b", "feature_a"],
            "status": ["unique", "ambiguous", "unique"],
        }
    )
    adata = dotmatch.assignments_to_anndata(assignments, feature_col="target_id")

    assert list(adata.obs_names) == ["cell_1", "cell_2"]
    assert list(adata.var_names) == ["feature_a"]
    from scipy import sparse
    assert sparse.isspmatrix_csr(adata.X)
    assert adata.X.toarray().tolist() == [[1], [1]]
