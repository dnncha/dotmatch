import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "check_download_metrics.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("check_download_metrics", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_artifact_breakdown_preserves_version_platform_and_python() -> None:
    script = _load_script()
    result = script.aggregate_files(
        [
            {
                "version": "0.2.2",
                "basename": "linux-64/dotmatch-0.2.2-py311h123_0.conda",
                "ndownloads": 12,
            },
            {
                "version": "0.2.2",
                "basename": "osx-arm64/dotmatch-0.2.2-py312h456_0.conda",
                "ndownloads": 8,
            },
        ]
    )

    assert result["total_downloads"] == 20
    assert result["artifact_files"] == 2
    assert result["by_version"]["0.2.2"] == {"downloads": 20, "artifact_files": 2}
    assert result["by_platform"]["linux-64"] == {"downloads": 12, "artifact_files": 1}
    assert result["by_python"]["3.11"] == {"downloads": 12, "artifact_files": 1}
    assert result["by_python"]["3.12"] == {"downloads": 8, "artifact_files": 1}


def test_snapshot_keeps_downloads_separate_from_discovery_signals(monkeypatch) -> None:
    script = _load_script()
    responses = {
        script.ANACONDA_PACKAGE_URL: {
            "files": [
                {
                    "version": "0.2.2",
                    "basename": "linux-64/dotmatch-0.2.2-py311h123_0.conda",
                    "ndownloads": 12,
                }
            ]
        },
        script.ANACONDA_REPOCORE_URL: {"latest_version": "0.2.2", "download_count": 12},
        script.PYPI_JSON_URL: {"info": {"version": "0.2.2"}},
        script.PYPISTATS_RECENT_URL: {"data": {"last_month": 4}},
        script.PYPISTATS_OVERALL_URL: {
            "data": [
                {"category": "without_mirrors", "downloads": 7},
                {"category": "without_mirrors", "downloads": 5},
                {"category": "with_mirrors", "downloads": 100},
            ]
        },
        script.BIOCONTAINERS_STATS_URL: {
            "stats": [{"date": "2026-08-01", "count": 2}],
            "tags": {"0.2.2--py311h123_1": {}},
        },
        script.ZENODO_RECORD_URL: {
            "id": 21511337,
            "doi": "10.5281/zenodo.21511337",
            "metadata": {"version": "0.2.2"},
            "stats": {"downloads": 1, "unique_downloads": 1, "views": 3},
            "files": [{"key": "dnncha/dotmatch-v0.2.2.zip", "size": 12}],
        },
        script.GITHUB_REPO_URL: {"stargazers_count": 2, "forks_count": 0},
        script.GITHUB_RELEASES_URL: [
            {
                "tag_name": "v0.2.2",
                "draft": False,
                "published_at": "2026-07-23T13:45:21Z",
                "html_url": "https://github.com/dnncha/dotmatch/releases/tag/v0.2.2",
                "assets": [{"name": "dotmatch-0.2.2.tar.gz", "download_count": 3}],
            }
        ],
    }
    monkeypatch.setattr(script, "fetch_json", lambda url: responses[url])

    snapshot = script.collect_snapshot()

    assert snapshot["aggregate_downloads"]["reported_downloads"] == 30
    assert snapshot["channels"]["anaconda_bioconda"]["breakdown"]["total_downloads"] == 12
    assert snapshot["channels"]["pypi"]["cumulative_downloads_excluding_mirrors"] == 12
    assert snapshot["channels"]["pypi"]["daily_downloads_excluding_mirrors"] == [
        {"date": "unknown", "downloads": 7},
        {"date": "unknown", "downloads": 5},
    ]
    assert snapshot["channels"]["ghcr"]["included_in_aggregate_total"] is False
    assert snapshot["channels"]["github_releases"]["breakdown"]["by_version"]["0.2.2"] == {
        "downloads": 3,
        "asset_files": 1,
    }
    assert snapshot["channels"]["biocontainers"]["breakdown"]["total_downloads"] == 2
    assert snapshot["channels"]["biocontainers"]["breakdown"]["by_version"] == {}
    assert snapshot["channels"]["biocontainers"]["included_in_aggregate_total"] is True
    assert snapshot["channels"]["zenodo"]["breakdown"]["total_downloads"] == 1
    assert snapshot["channels"]["zenodo"]["breakdown"]["by_version"]["0.2.2"]["downloads"] == 1
    assert snapshot["channels"]["zenodo"]["included_in_aggregate_total"] is True
    assert snapshot["channels"]["pypi"]["recent_month_downloads"] == 4
    assert snapshot["discovery_signals"]["github"]["stars"] == 2
    assert snapshot["interpretation"]["human_use_requires_separate_evidence"] is True


def test_pypi_series_retains_rolled_off_days_and_refreshes_overlap() -> None:
    script = _load_script()

    result = script.merge_pypi_daily_downloads(
        [
            {"date": "2026-08-02", "downloads": 7},
            {"date": "2026-08-03", "downloads": 5},
        ],
        {
            "channels": {
                "pypi": {
                    "daily_downloads_excluding_mirrors": [
                        {"date": "2026-07-01", "downloads": 11},
                        {"date": "2026-08-02", "downloads": 2},
                    ]
                }
            }
        },
    )

    assert result == [
        {"date": "2026-07-01", "downloads": 11},
        {"date": "2026-08-02", "downloads": 7},
        {"date": "2026-08-03", "downloads": 5},
    ]


def test_biocontainers_stats_preserves_window_and_missing_dimensions() -> None:
    script = _load_script()
    result = script.aggregate_biocontainers_stats(
        {
            "stats": [
                {"date": "2026-08-01", "count": 2},
                {"date": "2026-08-02", "count": 5},
            ],
            "tags": {"0.2.2--py311h123_1": {}},
        }
    )

    assert result["availability"] == "available"
    assert result["total_downloads"] == 7
    assert result["observed_start"] == "2026-08-01"
    assert result["observed_end"] == "2026-08-02"
    assert result["tag_names"] == ["0.2.2--py311h123_1"]
    assert result["by_version"] == {}
    assert result["by_platform"] == {}


def test_zenodo_stats_keeps_record_version_and_archive_platform() -> None:
    script = _load_script()
    result = script.aggregate_zenodo_stats(
        {
            "id": 21511337,
            "doi": "10.5281/zenodo.21511337",
            "metadata": {"version": "0.2.2"},
            "stats": {"downloads": 8, "unique_downloads": 6, "views": 92},
            "files": [{"key": "dnncha/dotmatch-v0.2.2.zip", "size": 12}],
        }
    )

    assert result["availability"] == "available"
    assert result["total_downloads"] == 8
    assert result["provider_unique_downloads"] == 6
    assert result["by_version"]["0.2.2"] == {"downloads": 8, "artifact_files": 1}
    assert result["by_platform"]["source_or_archive"] == {
        "downloads": 8,
        "artifact_files": 1,
    }


def test_release_asset_aggregation_ignores_drafts_and_breaks_down_platform() -> None:
    script = _load_script()
    result = script.aggregate_github_release_assets(
        [
            {
                "tag_name": "v0.2.2",
                "draft": False,
                "published_at": "2026-07-23T13:45:21Z",
                "assets": [
                    {"name": "dotmatch-0.2.2-py3-none-linux_x86_64.whl", "download_count": 4},
                    {"name": "dotmatch-0.2.2.tar.gz", "download_count": 2},
                ],
            },
            {
                "tag_name": "v0.3.0",
                "draft": True,
                "published_at": None,
                "assets": [{"name": "draft.tar.gz", "download_count": 100}],
            },
        ]
    )

    assert result["total_downloads"] == 6
    assert result["release_count"] == 1
    assert result["by_platform"]["linux"] == {"downloads": 4, "asset_files": 1}
    assert result["by_platform"]["source_or_archive"] == {"downloads": 2, "asset_files": 1}


def test_fetch_json_retries_transient_http_errors(monkeypatch) -> None:
    script = _load_script()
    attempts = {"count": 0}

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self):
            return b'{"ok": true}'

    def flaky_urlopen(request, timeout):
        attempts["count"] += 1
        if attempts["count"] < 3:
            raise script.HTTPError(request.full_url, 504, "gateway timeout", {}, None)
        return Response()

    monkeypatch.setattr(script.urllib.request, "urlopen", flaky_urlopen)
    monkeypatch.setattr(script.time, "sleep", lambda _: None)

    assert script.fetch_json("https://example.test/metrics") == {"ok": True}
    assert attempts["count"] == 3


def test_optional_fetch_records_provider_failure(monkeypatch) -> None:
    script = _load_script()

    def unavailable(_url):
        raise RuntimeError("provider unavailable")

    monkeypatch.setattr(script, "fetch_json", unavailable)

    result = script.fetch_optional_json("https://example.test/metrics")

    assert result["_fetch_error"] == "RuntimeError: provider unavailable"
