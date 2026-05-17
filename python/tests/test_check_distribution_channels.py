import importlib.util
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "scripts" / "check_distribution_channels.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_distribution_channels", CHECKER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_repo(root: Path, doi: str = "10.5281/zenodo.1234567") -> None:
    files = {
        "pyproject.toml": '[project]\nname = "dotmatch"\nversion = "0.1.0"\n',
        "CITATION.cff": f'cff-version: 1.2.0\ntitle: "DotMatch"\nversion: "0.1.0"\ndoi: {doi}\n',
    }
    for path, text in files.items():
        full = root / path
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text(text, encoding="utf-8")


def test_distribution_channels_accepts_mocked_public_release(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_repo(tmp_path)

    def fake_fetch_json(url: str):
        if url == "https://pypi.org/pypi/dotmatch/0.1.0/json":
            return {
                "info": {"version": "0.1.0"},
                "urls": [
                    {"packagetype": "sdist", "filename": "dotmatch-0.1.0.tar.gz"},
                    {
                        "packagetype": "bdist_wheel",
                        "filename": "dotmatch-0.1.0-py3-none-macosx_11_0_universal2.whl",
                    },
                    {
                        "packagetype": "bdist_wheel",
                        "filename": "dotmatch-0.1.0-py3-none-manylinux_2_28_x86_64.whl",
                    },
                    {
                        "packagetype": "bdist_wheel",
                        "filename": "dotmatch-0.1.0-py3-none-musllinux_1_2_x86_64.whl",
                    },
                ],
            }
        if url == "https://api.anaconda.org/package/bioconda/dotmatch":
            return {"files": [{"version": "0.1.0", "basename": "dotmatch-0.1.0-0.tar.bz2"}]}
        if url == "https://quay.io/api/v1/repository/biocontainers/dotmatch/tag/?onlyActiveTags=true&page=1&limit=100":
            return {"tags": [{"name": "0.1.0--h123_0"}], "has_additional": False}
        raise AssertionError(url)

    completed = subprocess.CompletedProcess(["docker"], 0, stdout="ok", stderr="")
    monkeypatch.setattr(checker, "fetch_json", fake_fetch_json)
    monkeypatch.setattr(checker.subprocess, "run", lambda *args, **kwargs: completed)
    monkeypatch.setattr(checker, "verify_pypi_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker, "verify_bioconda_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker, "verify_ghcr_run", lambda image, version: None, raising=False)
    monkeypatch.setattr(checker, "verify_biocontainers_run", lambda image, version: None, raising=False)
    monkeypatch.setattr(checker, "url_ok", lambda url: url == "https://doi.org/10.5281/zenodo.1234567")

    result = checker.audit(tmp_path)

    assert result.failures == []
    assert {item.channel for item in result.passed} == {
        "pypi",
        "pypi-install",
        "bioconda",
        "bioconda-install",
        "biocontainers",
        "biocontainers-run",
        "ghcr",
        "ghcr-run",
        "zenodo",
    }


def test_distribution_channels_reports_failed_pypi_install(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_repo(tmp_path)

    def fake_fetch_json(url: str):
        if "pypi.org" in url:
            return {
                "info": {"version": "0.1.0"},
                "urls": [
                    {"packagetype": "sdist", "filename": "dotmatch-0.1.0.tar.gz"},
                    {
                        "packagetype": "bdist_wheel",
                        "filename": "dotmatch-0.1.0-py3-none-macosx_11_0_universal2.whl",
                    },
                    {
                        "packagetype": "bdist_wheel",
                        "filename": "dotmatch-0.1.0-py3-none-manylinux_2_28_x86_64.whl",
                    },
                    {
                        "packagetype": "bdist_wheel",
                        "filename": "dotmatch-0.1.0-py3-none-musllinux_1_2_x86_64.whl",
                    },
                ],
            }
        return {"files": [{"version": "0.1.0"}]}

    monkeypatch.setattr(checker, "fetch_json", fake_fetch_json)
    monkeypatch.setattr(checker, "verify_pypi_install", lambda version: (_ for _ in ()).throw(RuntimeError("pip failed")), raising=False)
    monkeypatch.setattr(checker, "verify_bioconda_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker.subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(["docker"], 0))
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any("PyPI one-command install failed for 0.1.0" in failure.message for failure in result.failures)


def test_distribution_channels_reports_missing_zenodo_doi(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_repo(tmp_path, doi="")

    monkeypatch.setattr(checker, "fetch_json", lambda url: {})
    monkeypatch.setattr(checker.subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(["docker"], 0))

    result = checker.audit(tmp_path)

    assert any("CITATION.cff must include a DOI" in failure.message for failure in result.failures)


def test_distribution_channels_reports_pypi_missing_version(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_repo(tmp_path)

    def fake_fetch_json(url: str):
        if "pypi.org" in url:
            return {"info": {"version": "0.0.9"}, "urls": [{"packagetype": "bdist_wheel"}]}
        return {"files": [{"version": "0.1.0"}]}

    monkeypatch.setattr(checker, "fetch_json", fake_fetch_json)
    monkeypatch.setattr(checker, "verify_pypi_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker, "verify_bioconda_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker.subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(["docker"], 0))
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any("PyPI version 0.1.0 is not available as an sdist" in failure.message for failure in result.failures)


def test_distribution_channels_reports_pypi_missing_repaired_wheels(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_repo(tmp_path)

    def fake_fetch_json(url: str):
        if "pypi.org" in url:
            return {
                "info": {"version": "0.1.0"},
                "urls": [
                    {"packagetype": "sdist", "filename": "dotmatch-0.1.0.tar.gz"},
                    {
                        "packagetype": "bdist_wheel",
                        "filename": "dotmatch-0.1.0-py3-none-macosx_11_0_universal2.whl",
                    },
                ],
            }
        return {"files": [{"version": "0.1.0"}]}

    monkeypatch.setattr(checker, "fetch_json", fake_fetch_json)
    monkeypatch.setattr(checker, "verify_pypi_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker, "verify_bioconda_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker.subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(["docker"], 0))
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any("must include repaired manylinux and musllinux wheels" in failure.message for failure in result.failures)


def test_distribution_channels_rejects_raw_pypi_linux_wheel(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_repo(tmp_path)

    def fake_fetch_json(url: str):
        if "pypi.org" in url:
            return {
                "info": {"version": "0.1.0"},
                "urls": [
                    {"packagetype": "sdist", "filename": "dotmatch-0.1.0.tar.gz"},
                    {
                        "packagetype": "bdist_wheel",
                        "filename": "dotmatch-0.1.0-py3-none-macosx_11_0_universal2.whl",
                    },
                    {
                        "packagetype": "bdist_wheel",
                        "filename": "dotmatch-0.1.0-py3-none-manylinux_2_28_x86_64.whl",
                    },
                    {
                        "packagetype": "bdist_wheel",
                        "filename": "dotmatch-0.1.0-py3-none-musllinux_1_2_x86_64.whl",
                    },
                    {"packagetype": "bdist_wheel", "filename": "dotmatch-0.1.0-py3-none-linux_x86_64.whl"},
                ],
            }
        return {"files": [{"version": "0.1.0"}]}

    monkeypatch.setattr(checker, "fetch_json", fake_fetch_json)
    monkeypatch.setattr(checker, "verify_pypi_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker, "verify_bioconda_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker.subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(["docker"], 0))
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any("must not include raw linux_x86_64 wheels" in failure.message for failure in result.failures)


def test_distribution_channels_reports_bioconda_missing_version(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_repo(tmp_path)

    def fake_fetch_json(url: str):
        if "pypi.org" in url:
            return {"info": {"version": "0.1.0"}, "urls": [{"packagetype": "sdist"}]}
        return {"files": [{"version": "0.0.9"}]}

    monkeypatch.setattr(checker, "fetch_json", fake_fetch_json)
    monkeypatch.setattr(checker, "verify_pypi_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker.subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(["docker"], 0))
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any("Bioconda version 0.1.0 is not available" in failure.message for failure in result.failures)


def test_distribution_channels_reports_failed_bioconda_install(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_repo(tmp_path)

    def fake_fetch_json(url: str):
        if "pypi.org" in url:
            return {"info": {"version": "0.1.0"}, "urls": [{"packagetype": "sdist"}]}
        return {"files": [{"version": "0.1.0"}]}

    monkeypatch.setattr(checker, "fetch_json", fake_fetch_json)
    monkeypatch.setattr(checker, "verify_pypi_install", lambda version: None, raising=False)
    monkeypatch.setattr(
        checker,
        "verify_bioconda_install",
        lambda version: (_ for _ in ()).throw(RuntimeError("conda failed")),
        raising=False,
    )
    monkeypatch.setattr(checker.subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(["docker"], 0))
    monkeypatch.setattr(checker, "verify_ghcr_run", lambda image, version: None, raising=False)
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any("Bioconda one-command install failed for 0.1.0" in failure.message for failure in result.failures)


def test_distribution_channels_reports_missing_biocontainers_tag(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_repo(tmp_path)

    def fake_fetch_json(url: str):
        if "pypi.org" in url:
            return {"info": {"version": "0.1.0"}, "urls": [{"packagetype": "sdist"}]}
        if "api.anaconda.org" in url:
            return {"files": [{"version": "0.1.0"}]}
        return {"tags": [{"name": "0.0.9--h123_0"}], "has_additional": False}

    monkeypatch.setattr(checker, "fetch_json", fake_fetch_json)
    monkeypatch.setattr(checker, "verify_pypi_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker, "verify_bioconda_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker.subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(["docker"], 0))
    monkeypatch.setattr(checker, "verify_ghcr_run", lambda image, version: None, raising=False)
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any("BioContainers image tag for version 0.1.0 is not available" in failure.message for failure in result.failures)


def test_distribution_channels_reports_failed_biocontainers_runtime_smoke(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_repo(tmp_path)

    def fake_fetch_json(url: str):
        if "pypi.org" in url:
            return {"info": {"version": "0.1.0"}, "urls": [{"packagetype": "sdist"}]}
        if "api.anaconda.org" in url:
            return {"files": [{"version": "0.1.0"}]}
        return {"tags": [{"name": "0.1.0--h123_0"}], "has_additional": False}

    monkeypatch.setattr(checker, "fetch_json", fake_fetch_json)
    monkeypatch.setattr(checker, "verify_pypi_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker, "verify_bioconda_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker.subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(["docker"], 0))
    monkeypatch.setattr(
        checker,
        "verify_biocontainers_run",
        lambda image, version: (_ for _ in ()).throw(RuntimeError("biocontainer failed")),
        raising=False,
    )
    monkeypatch.setattr(checker, "verify_ghcr_run", lambda image, version: None, raising=False)
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any(
        "BioContainers image runtime smoke test failed for quay.io/biocontainers/dotmatch:0.1.0--h123_0"
        in failure.message
        for failure in result.failures
    )


def test_distribution_channels_reports_missing_ghcr_manifest(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_repo(tmp_path)

    def fake_fetch_json(url: str):
        if "pypi.org" in url:
            return {"info": {"version": "0.1.0"}, "urls": [{"packagetype": "sdist"}]}
        return {"files": [{"version": "0.1.0"}]}

    failed = subprocess.CompletedProcess(["docker"], 1, stdout="", stderr="manifest unknown")
    monkeypatch.setattr(checker, "fetch_json", fake_fetch_json)
    monkeypatch.setattr(checker, "verify_pypi_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker, "verify_bioconda_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker.subprocess, "run", lambda *args, **kwargs: failed)
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any("GHCR image tag ghcr.io/dnncha/dotmatch:v0.1.0 is not available" in failure.message for failure in result.failures)


def test_distribution_channels_reports_failed_ghcr_runtime_smoke(tmp_path, monkeypatch):
    checker = _load_checker()
    _write_repo(tmp_path)

    def fake_fetch_json(url: str):
        if "pypi.org" in url:
            return {"info": {"version": "0.1.0"}, "urls": [{"packagetype": "sdist"}]}
        return {"files": [{"version": "0.1.0"}]}

    monkeypatch.setattr(checker, "fetch_json", fake_fetch_json)
    monkeypatch.setattr(checker, "verify_pypi_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker, "verify_bioconda_install", lambda version: None, raising=False)
    monkeypatch.setattr(checker.subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(["docker"], 0))
    monkeypatch.setattr(
        checker,
        "verify_ghcr_run",
        lambda image, version: (_ for _ in ()).throw(RuntimeError("container failed")),
        raising=False,
    )
    monkeypatch.setattr(checker, "url_ok", lambda url: True)

    result = checker.audit(tmp_path)

    assert any(
        "GHCR image runtime smoke test failed for ghcr.io/dnncha/dotmatch:v0.1.0" in failure.message
        for failure in result.failures
    )
