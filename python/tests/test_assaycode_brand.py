from __future__ import annotations

import assaycode
from assaycode import cli


def test_platform_identity_uses_dotmatch_engine() -> None:
    assert assaycode.PLATFORM_NAME == "AssayCode"
    assert assaycode.SPEC_NAME == "AssayScript"
    assert assaycode.ENGINE_NAME == "DotMatch"
    assert assaycode.__version__ == assaycode.engine.__version__


def test_assaycode_help_states_compatibility_boundary(capsys) -> None:
    assert cli.main(["--help"]) == 0
    output = capsys.readouterr().out
    assert "powered by the DotMatch engine" in output
    assert "additive assay-level identity" in output
    assert "assaycode check assay.toml" in output


def test_assaycode_version_names_platform_and_engine(capsys) -> None:
    assert cli.main(["--version"]) == 0
    output = capsys.readouterr().out
    assert output.startswith(f"assaycode {assaycode.__version__}")
    assert f"DotMatch engine {assaycode.__version__}" in output


def test_assaycode_shortcut_delegates_to_assay_namespace(monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(cli._engine_cli, "main", lambda argv: calls.append(list(argv)) or 0)
    assert cli.main(["check", "assay.toml"]) == 0
    assert calls == [["assay", "check", "assay.toml"]]


def test_assaycode_engine_escape_hatch_delegates_unchanged(monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(cli._engine_cli, "main", lambda argv: calls.append(list(argv)) or 0)
    assert cli.main(["engine", "dist", "ACGT", "AGGT"]) == 0
    assert calls == [["dist", "ACGT", "AGGT"]]


def test_assaycode_specialized_namespace_delegates_unchanged(monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(cli._engine_cli, "main", lambda argv: calls.append(list(argv)) or 0)
    assert cli.main(["barcode", "infer", "--help"]) == 0
    assert calls == [["barcode", "infer", "--help"]]
